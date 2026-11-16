from __future__ import annotations

import hashlib
import itertools
from pathlib import Path

import pytest

from cellulose_forge.energy_drivers.openmm import (
    OpenMMForceClassification,
    openmm_energy_report,
    run_openmm_single_point,
)
from cellulose_forge.energy_drivers.report import validate_potential_sum


def test_openmm_components_map_to_canonical_terms() -> None:
    report = openmm_energy_report(
        {
            "bond": 1.0,
            "angle": 2.0,
            "urey_bradley": 3.0,
            "proper_dihedral": 4.0,
            "periodic_improper": 5.0,
            "harmonic_improper": 6.0,
            "vdw": 7.0,
            "vdw_14": 8.0,
            "electrostatics": 9.0,
            "electrostatics_14": 10.0,
            "restraint": 11.0,
            "dispersion_correction": 12.0,
            "potential": 78.0,
        },
        coordinate_sha256="abc",
    )

    assert report.terms == {
        "bond": 1.0,
        "angle": 5.0,
        "dihedral": 4.0,
        "improper": 11.0,
        "vdw": 15.0,
        "electrostatics": 19.0,
        "restraint": 11.0,
        "dispersion_correction": 12.0,
        "potential": 78.0,
    }
    validate_potential_sum(report, absolute_kj_mol=0.0, relative=0.0)


def test_openmm_missing_components_are_explicit_zeros() -> None:
    report = openmm_energy_report({"potential": 0}, coordinate_sha256="abc")

    assert set(report.terms.values()) == {0}


@pytest.mark.parametrize("components", [{"bond": 1.0}, {"potential": 0, "cmap": 1.0}])
def test_openmm_refuses_incomplete_or_unclassified_components(
    components: dict[str, float],
) -> None:
    with pytest.raises(ValueError):
        openmm_energy_report(components, coordinate_sha256="abc")


def test_openmm_extra_metadata_is_copied() -> None:
    metadata: dict[str, object] = {"platform": "CPU"}
    report = openmm_energy_report({"potential": 0}, coordinate_sha256="abc", extra=metadata)
    metadata["platform"] = "Reference"

    assert report.extra["platform"] == "CPU"


def _write_openmm_fixture(directory: Path) -> tuple[Path, Path, str]:
    openmm = pytest.importorskip("openmm")
    app = pytest.importorskip("openmm.app")

    system = openmm.System()
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("MOL", chain)
    atoms = []
    for name in ("C1", "C2", "C3", "C4"):
        system.addParticle(12.0)
        atoms.append(topology.addAtom(name, app.Element.getBySymbol("C"), residue))
    for first, second in itertools.pairwise(atoms):
        topology.addBond(first, second)

    bonds = openmm.HarmonicBondForce()
    bonds.addBond(0, 1, 0.14, 1000.0)
    bonds.addBond(0, 2, 0.25, 500.0)  # Urey--Bradley.
    system.addForce(bonds)
    angles = openmm.HarmonicAngleForce()
    angles.addAngle(0, 1, 2, 1.9, 100.0)
    system.addForce(angles)
    torsions = openmm.PeriodicTorsionForce()
    torsions.addTorsion(0, 1, 2, 3, 1, 0.0, 2.0)
    torsions.addTorsion(1, 0, 2, 3, 2, 3.14, 1.0)  # Improper.
    system.addForce(torsions)
    nonbonded = openmm.NonbondedForce()
    nonbonded.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
    for charge in (0.2, -0.2, 0.1, -0.1):
        nonbonded.addParticle(charge, 0.3, 0.2)
    nonbonded.addException(0, 1, 0.0, 0.3, 0.0)
    nonbonded.addException(1, 2, 0.0, 0.3, 0.0)
    nonbonded.addException(2, 3, 0.0, 0.3, 0.0)
    system.addForce(nonbonded)
    restraint = openmm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", 10.0)
    for parameter in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(parameter)
    restraint.addParticle(0, [0.0, 0.0, 0.0])
    system.addForce(restraint)
    system.addForce(openmm.CMMotionRemover())

    system_path = directory / "system.xml"
    system_path.write_text(openmm.XmlSerializer.serialize(system), encoding="utf-8")
    coordinate_path = directory / "coordinates.gro"
    coordinate_path.write_text(
        "extended-precision fixture\n"
        "4\n"
        "    1MOL     C1    1   0.01000000   0.00000000   0.00000000\n"
        "    1MOL     C2    2   0.16000000   0.00000000   0.00000000\n"
        "    1MOL     C3    3   0.27000000   0.08000000   0.00000000\n"
        "    1MOL     C4    4   0.39000000   0.10000000   0.08000000\n"
        "  2.0000000  2.1000000  2.2000000  0.0000000  0.0000000"
        " -0.1000000  0.0000000  0.0000000  0.0000000\n",
        encoding="utf-8",
    )
    checksum = hashlib.sha256(coordinate_path.read_bytes()).hexdigest()
    return system_path, coordinate_path, checksum


def test_run_openmm_single_point_splits_mixed_forces(tmp_path: Path) -> None:
    system, coordinates, checksum = _write_openmm_fixture(tmp_path)
    report = run_openmm_single_point(
        system_file=system,
        coordinate_file=coordinates,
        coordinate_sha256=checksum,
        classification=OpenMMForceClassification(
            urey_bradley_bonds={0: (1,)},
            periodic_impropers={2: (1,)},
            force_roles={4: "restraint"},
        ),
    )

    assert report.terms["bond"] != 0
    assert report.terms["angle"] != 0
    assert report.terms["dihedral"] != 0
    assert report.terms["improper"] != 0
    assert report.terms["vdw"] != 0
    assert report.terms["electrostatics"] != 0
    assert report.terms["restraint"] != 0
    assert report.terms["dispersion_correction"] == pytest.approx(0.0)
    validate_potential_sum(report, absolute_kj_mol=1.0e-6, relative=1.0e-12)


def test_run_openmm_single_point_refuses_unclassified_custom_force(tmp_path: Path) -> None:
    system, coordinates, checksum = _write_openmm_fixture(tmp_path)

    with pytest.raises(ValueError, match="requires an explicit role"):
        run_openmm_single_point(
            system_file=system,
            coordinate_file=coordinates,
            coordinate_sha256=checksum,
            classification=OpenMMForceClassification(
                urey_bradley_bonds={0: (1,)}, periodic_impropers={2: (1,)}
            ),
        )


def test_run_openmm_single_point_checks_coordinate_hash(tmp_path: Path) -> None:
    system, coordinates, _ = _write_openmm_fixture(tmp_path)

    with pytest.raises(ValueError, match="shared checksum"):
        run_openmm_single_point(
            system_file=system,
            coordinate_file=coordinates,
            coordinate_sha256="wrong",
            classification=OpenMMForceClassification(),
        )
