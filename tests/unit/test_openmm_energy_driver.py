from __future__ import annotations

import hashlib
import itertools
from pathlib import Path

import pytest

from cellulose_forge.energy_drivers.coordinates import write_coordinate_reference
from cellulose_forge.energy_drivers.openmm import (
    OpenMMForceClassification,
    openmm_classification_from_topology,
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


def _write_openmm_fixture(directory: Path) -> tuple[Path, Path, Path, str]:
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
    bonds.addBond(1, 2, 0.14, 1000.0)
    bonds.addBond(2, 3, 0.14, 1000.0)
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
    rb = openmm.RBTorsionForce()
    rb.addTorsion(0, 1, 2, 3, 1.0, 0.5, 0.25, 0.125, 0.0, 0.0)
    system.addForce(rb)
    system.addForce(openmm.CMMotionRemover())

    system_path = directory / "system.xml"
    system_path.write_text(openmm.XmlSerializer.serialize(system), encoding="utf-8")
    coordinate_path = directory / "coordinates.gro"
    coordinate_path.write_text(
        "extended-precision fixture\n"
        "4\n"
        "    1MOL     C1    1   0.01000000   0.00000000   0.00000000  9.0  8.0  7.0\n"
        "    1MOL     C2    2   0.16000000   0.00000000   0.00000000  9.0  8.0  7.0\n"
        "    1MOL     C3    3   0.27000000   0.08000000   0.00000000  9.0  8.0  7.0\n"
        "    1MOL     C4    4   0.39000000   0.10000000   0.08000000  9.0  8.0  7.0\n"
        "  2.0000000  2.1000000  2.2000000  0.0000000  0.0000000"
        " -0.1000000  0.0000000  0.0000000  0.0000000\n",
        encoding="utf-8",
    )
    reference = directory / "coordinate_reference.json"
    checksum = write_coordinate_reference(
        reference,
        [[0.1, 0.0, 0.0], [1.6, 0.0, 0.0], [2.7, 0.8, 0.0], [3.9, 1.0, 0.8]],
        [[20.0, 0.0, 0.0], [-1.0, 21.0, 0.0], [0.0, 0.0, 22.0]],
    )
    return system_path, coordinate_path, reference, checksum


def test_run_openmm_single_point_splits_mixed_forces(tmp_path: Path) -> None:
    system, coordinates, reference, checksum = _write_openmm_fixture(tmp_path)
    report = run_openmm_single_point(
        system_file=system,
        coordinate_file=coordinates,
        coordinate_reference_file=reference,
        coordinate_sha256=checksum,
        classification=OpenMMForceClassification(
            urey_bradley_bonds={0: (1,)},
            periodic_impropers={2: (1,)},
            rb_impropers={5: ()},
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
    system, coordinates, reference, checksum = _write_openmm_fixture(tmp_path)

    with pytest.raises(ValueError, match="requires an explicit role"):
        run_openmm_single_point(
            system_file=system,
            coordinate_file=coordinates,
            coordinate_reference_file=reference,
            coordinate_sha256=checksum,
            classification=OpenMMForceClassification(
                urey_bradley_bonds={0: (1,)},
                periodic_impropers={2: (1,)},
                rb_impropers={5: ()},
            ),
        )


def test_run_openmm_single_point_checks_coordinate_hash(tmp_path: Path) -> None:
    system, coordinates, reference, _ = _write_openmm_fixture(tmp_path)

    with pytest.raises(ValueError, match="manifest checksum"):
        run_openmm_single_point(
            system_file=system,
            coordinate_file=coordinates,
            coordinate_reference_file=reference,
            coordinate_sha256="wrong",
            classification=OpenMMForceClassification(),
        )


def test_run_openmm_single_point_requires_every_mixed_force_classification(
    tmp_path: Path,
) -> None:
    system, coordinates, reference, checksum = _write_openmm_fixture(tmp_path)

    with pytest.raises(ValueError, match="every HarmonicBondForce"):
        run_openmm_single_point(
            system_file=system,
            coordinate_file=coordinates,
            coordinate_reference_file=reference,
            coordinate_sha256=checksum,
            classification=OpenMMForceClassification(
                periodic_impropers={2: ()}, rb_impropers={5: ()}
            ),
        )


def test_run_openmm_single_point_refuses_wrong_force_classification(tmp_path: Path) -> None:
    system, coordinates, reference, checksum = _write_openmm_fixture(tmp_path)

    with pytest.raises(ValueError, match="every HarmonicBondForce"):
        run_openmm_single_point(
            system_file=system,
            coordinate_file=coordinates,
            coordinate_reference_file=reference,
            coordinate_sha256=checksum,
            classification=OpenMMForceClassification(
                urey_bradley_bonds={2: ()},
                periodic_impropers={2: ()},
                rb_impropers={5: ()},
            ),
        )


def test_run_openmm_single_point_refuses_periodic_system_without_coordinate_box(
    tmp_path: Path,
) -> None:
    openmm = pytest.importorskip("openmm")
    system_path, _, reference, _ = _write_openmm_fixture(tmp_path)
    system = openmm.XmlSerializer.deserialize(system_path.read_text(encoding="utf-8"))
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(2, 0, 0), openmm.Vec3(0, 2, 0), openmm.Vec3(0, 0, 2)
    )
    system.getForce(3).setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    system_path.write_text(openmm.XmlSerializer.serialize(system), encoding="utf-8")
    coordinates = tmp_path / "no_box.pdb"
    coordinates.write_text(
        "ATOM      1  C1  MOL A   1       0.100   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  C2  MOL A   1       1.600   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C3  MOL A   1       2.700   0.800   0.000  1.00  0.00           C\n"
        "ATOM      4  C4  MOL A   1       3.900   1.000   0.800  1.00  0.00           C\n"
        "END\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="require box vectors"):
        run_openmm_single_point(
            system_file=system_path,
            coordinate_file=coordinates,
            coordinate_reference_file=reference,
            coordinate_sha256=hashlib.sha256(coordinates.read_bytes()).hexdigest(),
            classification=OpenMMForceClassification(
                urey_bradley_bonds={0: (1,)},
                periodic_impropers={2: (1,)},
                rb_impropers={5: ()},
                force_roles={4: "restraint"},
            ),
        )


def test_run_openmm_single_point_refuses_nonbonded_parameter_offsets(tmp_path: Path) -> None:
    openmm = pytest.importorskip("openmm")
    system_path, coordinates, reference, checksum = _write_openmm_fixture(tmp_path)
    system = openmm.XmlSerializer.deserialize(system_path.read_text(encoding="utf-8"))
    nonbonded = system.getForce(3)
    nonbonded.addGlobalParameter("lambda", 0.0)
    nonbonded.addParticleParameterOffset("lambda", 0, 0.1, 0.0, 0.0)
    system_path.write_text(openmm.XmlSerializer.serialize(system), encoding="utf-8")

    with pytest.raises(ValueError, match="parameter offsets"):
        run_openmm_single_point(
            system_file=system_path,
            coordinate_file=coordinates,
            coordinate_reference_file=reference,
            coordinate_sha256=checksum,
            classification=OpenMMForceClassification(
                urey_bradley_bonds={0: (1,)},
                periodic_impropers={2: (1,)},
                rb_impropers={5: ()},
                force_roles={4: "restraint"},
            ),
        )


def test_openmm_classification_from_topology_finds_urey_bradley(tmp_path: Path) -> None:
    system, _, _, _ = _write_openmm_fixture(tmp_path)

    classification = openmm_classification_from_topology(
        system,
        covalent_bonds={(0, 1), (1, 2), (2, 3)},
        periodic_impropers={2: (1,)},
        rb_impropers={5: ()},
        force_roles={4: "restraint"},
    )

    assert classification.urey_bradley_bonds == {0: (1,)}


def test_openmm_classification_refuses_non_one_three_pseudo_bond(tmp_path: Path) -> None:
    openmm = pytest.importorskip("openmm")
    system_path, _, _, _ = _write_openmm_fixture(tmp_path)
    system = openmm.XmlSerializer.deserialize(system_path.read_text(encoding="utf-8"))
    system.getForce(0).addBond(0, 3, 0.4, 10.0)
    system_path.write_text(openmm.XmlSerializer.serialize(system), encoding="utf-8")

    with pytest.raises(ValueError, match="not a 1-3 Urey"):
        openmm_classification_from_topology(system_path, covalent_bonds={(0, 1), (1, 2), (2, 3)})


def test_openmm_classification_refuses_graph_bond_missing_from_force(tmp_path: Path) -> None:
    system_path, _, _, _ = _write_openmm_fixture(tmp_path)

    with pytest.raises(ValueError, match="absent from HarmonicBondForce"):
        openmm_classification_from_topology(
            system_path,
            covalent_bonds={(0, 1), (1, 2), (2, 3), (0, 3)},
        )


def test_run_openmm_single_point_refuses_constrained_system(tmp_path: Path) -> None:
    openmm = pytest.importorskip("openmm")
    system_path, coordinates, reference, checksum = _write_openmm_fixture(tmp_path)
    system = openmm.XmlSerializer.deserialize(system_path.read_text(encoding="utf-8"))
    system.addConstraint(0, 1, 0.14)
    system_path.write_text(openmm.XmlSerializer.serialize(system), encoding="utf-8")

    with pytest.raises(ValueError, match="must be unconstrained"):
        run_openmm_single_point(
            system_file=system_path,
            coordinate_file=coordinates,
            coordinate_reference_file=reference,
            coordinate_sha256=checksum,
            classification=OpenMMForceClassification(
                urey_bradley_bonds={0: (1,)},
                periodic_impropers={2: (1,)},
                rb_impropers={5: ()},
                force_roles={4: "restraint"},
            ),
        )
