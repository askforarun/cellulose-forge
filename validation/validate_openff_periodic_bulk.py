"""Validate experimental OpenFF Sage periodic cellulose bulk crystals.

The matrix is deliberately a model-construction validation.  It checks direct
graph mapping, stereochemistry, charges, dispersion policy, detailed OpenMM /
GROMACS energies, and short GROMACS minimization/NVT smoke tests.  It does not
assess Sage's physical accuracy for cellulose.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from contextlib import chdir
from importlib import metadata
from pathlib import Path

import numpy as np
import openmm
from openff.interchange.drivers.gromacs import _parse_gmx_energy
from openff.interchange.drivers.gromacs import _process as process_gromacs_energies
from openff.interchange.drivers.openmm import _get_openmm_energies, _process
from openff.toolkit import ForceField, Topology
from openff.units import unit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "validation"))

from validate_openff_finite_slab import (  # noqa: E402
    ENERGY_TOLERANCES_PER_ATOM_KJ_MOL,
    energies,
)

from cellulose_forge.crystallography import load_allomorph  # noqa: E402
from cellulose_forge.geometry import build_bulk_crystal  # noqa: E402
from cellulose_forge.openff_sage import (  # noqa: E402
    SAGE_FORCE_FIELD,
    build_sage_system,
    dispersion_correction_for_shape,
    reduced_periodic_box,
)

ALLOMORPHS = ("i_alpha", "i_beta", "ii", "iii_i")
MINIMUM_EDGE_ANGSTROM = 18.01
MINIMUM_ALLOWED_EDGE_ANGSTROM = 18.0
MAXIMUM_NVT_BOND_CHANGE_ANGSTROM = 0.2
POINT_MDP = """integrator = md
nsteps = 0
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 0.9
fourierspacing = 0.05
ewald-rtol = 1e-6
pme-order = 8
vdwtype = Cut-off
vdw-modifier = Potential-switch
rvdw-switch = 0.8
rvdw = 0.9
DispCorr = no
constraints = none
periodic-molecules = yes
nstenergy = 1
"""
MINIM_MDP = (
    POINT_MDP.replace(
        "integrator = md\nnsteps = 0", "integrator = steep\nemtol = 1000\nnsteps = 5000"
    )
    .replace("constraints = none", "constraints = h-bonds")
    .replace("DispCorr = no", "DispCorr = EnerPres")
)
NVT_MDP = (
    POINT_MDP.replace("nsteps = 0", "dt = 0.001\nnsteps = 1000")
    .replace("constraints = none", "constraints = h-bonds")
    .replace("DispCorr = no", "DispCorr = EnerPres")
    + """tcoupl = v-rescale
tc-grps = System
tau-t = 1.0
ref-t = 300
pcoupl = no
gen-vel = yes
gen-temp = 300
gen-seed = 20261112
nstlog = 100
nstxout-compressed = 100
"""
)


def _replication(allomorph: str) -> tuple[int, int, int]:
    vectors = load_allomorph(allomorph).unit_cell.vectors()
    return tuple(  # type: ignore[return-value]
        max(1, math.ceil(MINIMUM_EDGE_ANGSTROM / abs(float(vectors[axis, axis]))))
        for axis in range(3)
    )


def _interchange(build: object) -> object:
    topology = Topology.from_molecules([build.molecule] * build.mapped.chain_count)
    topology.box_vectors = reduced_periodic_box(build.mapped.basis_vectors) * unit.angstrom
    interchange = ForceField(SAGE_FORCE_FIELD).create_interchange(
        topology, charge_from_molecules=[build.molecule]
    )
    interchange.positions = build.mapped.cartesian_coordinates * unit.angstrom
    return interchange


def _openmm_energies(interchange: object, periodic: bool = True) -> dict[str, float]:
    system = interchange.to_openmm_system(  # type: ignore[attr-defined]
        combine_nonbonded_forces=False,
        add_constrained_forces=True,
        ewald_tolerance=1.0e-6,
    )
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.setUseDispersionCorrection(False)
            force.setExceptionsUsePeriodicBoundaryConditions(periodic)
        elif isinstance(force, openmm.CustomNonbondedForce):
            force.setUseLongRangeCorrection(False)
        elif periodic and hasattr(force, "setUsesPeriodicBoundaryConditions"):
            force.setUsesPeriodicBoundaryConditions(True)
    previous_threads = os.environ.get("OPENMM_CPU_THREADS")
    os.environ["OPENMM_CPU_THREADS"] = "1"
    try:
        raw = _get_openmm_energies(
            system=system,
            box_vectors=interchange.box.to_openmm(),  # type: ignore[attr-defined]
            positions=interchange.positions.to_openmm(),  # type: ignore[attr-defined]
            round_positions=None,
            platform="CPU",
        )
    finally:
        if previous_threads is None:
            os.environ.pop("OPENMM_CPU_THREADS", None)
        else:
            os.environ["OPENMM_CPU_THREADS"] = previous_threads
    return energies(_process(raw, system=system, combine_nonbonded_forces=False, detailed=True))


def _gromacs_energies(interchange: object, directory: Path) -> dict[str, float]:
    prefix = str(directory / "structure")
    interchange.to_gromacs(prefix, decimal=8)  # type: ignore[attr-defined]
    mdp = directory / "point.mdp"
    mdp.write_text(POINT_MDP, encoding="utf-8")
    with chdir(directory):
        subprocess.run(
            [
                "gmx",
                "grompp",
                "-maxwarn",
                "0",
                "-o",
                "out.tpr",
                "-f",
                str(mdp),
                "-c",
                f"{prefix}.gro",
                "-p",
                f"{prefix}.top",
            ],
            text=True,
            capture_output=True,
            check=True,
        )
        subprocess.run(
            ["gmx", "mdrun", "-s", "out.tpr", "-e", "out.edr", "-ntmpi", "1", "-ntomp", "1"],
            text=True,
            capture_output=True,
            check=True,
        )
        raw = _parse_gmx_energy("out.edr")
    values = energies(process_gromacs_energies(raw, detailed=True))
    return values


def _run_gromacs(interchange: object, directory: Path) -> dict[str, object]:
    prefix = str(directory / "structure")
    interchange.to_gromacs(prefix, decimal=8)  # type: ignore[attr-defined]
    (directory / "minim.mdp").write_text(MINIM_MDP, encoding="utf-8")
    (directory / "nvt.mdp").write_text(NVT_MDP, encoding="utf-8")
    commands = (
        [
            "gmx",
            "grompp",
            "-f",
            "minim.mdp",
            "-c",
            "structure.gro",
            "-p",
            "structure.top",
            "-o",
            "minim.tpr",
            "-maxwarn",
            "0",
        ],
        ["gmx", "mdrun", "-deffnm", "minim", "-ntmpi", "1"],
        [
            "gmx",
            "grompp",
            "-f",
            "nvt.mdp",
            "-c",
            "minim.gro",
            "-p",
            "structure.top",
            "-o",
            "nvt.tpr",
            "-maxwarn",
            "0",
        ],
        ["gmx", "mdrun", "-deffnm", "nvt", "-ntmpi", "1"],
    )
    logs = []
    for command in commands:
        result = subprocess.run(command, cwd=directory, text=True, capture_output=True, check=True)
        logs.append(result.stdout + result.stderr)
    combined = "\n".join(logs) + (directory / "nvt.log").read_text(encoding="utf-8")
    markers = [
        name
        for name, pattern in {
            "LINCS": r"LINCS WARNING",
            "NaN": r"\bnan\b",
            "non-finite": r"non-finite",
        }.items()
        if re.search(pattern, combined, flags=re.IGNORECASE)
    ]
    reached = "Statistics over 1001 steps" in combined
    return {
        "grompp_maxwarn": 0,
        "minimization_completed": (directory / "minim.gro").exists(),
        "nvt_duration_ps": 1.0,
        "nvt_reached_final_step": reached,
        "failure_markers": markers,
        "passed": (directory / "minim.gro").exists() and reached and not markers,
    }


def _minimum_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    fractional = delta @ np.linalg.inv(box)
    fractional -= np.round(fractional)
    return fractional @ box


def _bond_lengths(positions: np.ndarray, box: np.ndarray, build: object) -> np.ndarray:
    return np.asarray(
        [
            np.linalg.norm(
                _minimum_image(
                    positions[offset + bond.atom2_index] - positions[offset + bond.atom1_index],
                    box,
                )
            )
            for chain in range(build.mapped.chain_count)
            for offset in [chain * build.molecule.n_atoms]
            for bond in build.molecule.bonds
        ]
    )


def _run_openmm(build: object, interchange: object) -> dict[str, object]:
    system = interchange.to_openmm_system(  # type: ignore[attr-defined]
        combine_nonbonded_forces=False,
        add_constrained_forces=True,
        ewald_tolerance=1.0e-6,
    )
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.setUseDispersionCorrection(True)
            force.setExceptionsUsePeriodicBoundaryConditions(True)
        elif isinstance(force, openmm.CustomNonbondedForce):
            force.setUseLongRangeCorrection(True)
        elif hasattr(force, "setUsesPeriodicBoundaryConditions"):
            force.setUsesPeriodicBoundaryConditions(True)
    platform = openmm.Platform.getPlatformByName("CPU")
    properties = {"Threads": "1"}
    minimizer = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, minimizer, platform, properties)
    context.setPeriodicBoxVectors(*interchange.box.to_openmm())  # type: ignore[attr-defined]
    context.setPositions(interchange.positions.to_openmm())  # type: ignore[attr-defined]
    openmm.LocalEnergyMinimizer.minimize(context, tolerance=10.0, maxIterations=2000)
    minimized = context.getState(getPositions=True, getEnergy=True)
    del context, minimizer
    integrator = openmm.LangevinMiddleIntegrator(
        300 * openmm.unit.kelvin,
        1 / openmm.unit.picosecond,
        0.001 * openmm.unit.picoseconds,
    )
    context = openmm.Context(system, integrator, platform, properties)
    context.setPeriodicBoxVectors(*interchange.box.to_openmm())  # type: ignore[attr-defined]
    context.setPositions(minimized.getPositions())
    context.setVelocitiesToTemperature(300 * openmm.unit.kelvin, 20261112)
    integrator.step(1000)
    final = context.getState(getPositions=True, getEnergy=True)
    box = np.asarray(interchange.box.m_as(unit.angstrom))  # type: ignore[attr-defined]
    minimized_positions = minimized.getPositions(asNumpy=True).value_in_unit(openmm.unit.angstrom)
    final_positions = final.getPositions(asNumpy=True).value_in_unit(openmm.unit.angstrom)
    before = _bond_lengths(minimized_positions, box, build)
    after = _bond_lengths(final_positions, box, build)
    maximum_change = float(np.max(np.abs(after - before)))
    finite = bool(
        np.isfinite(final_positions).all()
        and np.isfinite(final.getPotentialEnergy().value_in_unit(openmm.unit.kilojoule_per_mole))
    )
    return {
        "platform": "CPU",
        "threads": 1,
        "duration_ps": 1.0,
        "bonds_checked": len(after),
        "maximum_bond_change_angstrom": maximum_change,
        "all_finite": finite,
        "passed": finite and maximum_change < MAXIMUM_NVT_BOND_CHANGE_ANGSTROM,
    }


def validate_case(allomorph: str) -> dict[str, object]:
    replication = _replication(allomorph)
    crystal = build_bulk_crystal(load_allomorph(allomorph), replication)
    if (
        min(abs(np.diag(reduced_periodic_box(crystal.basis_vectors))))
        <= MINIMUM_ALLOWED_EDGE_ANGSTROM
    ):
        raise RuntimeError("bulk box does not exceed twice the 9 A cutoff")
    build = build_sage_system(crystal, "periodic")
    interchange = _interchange(build)
    with tempfile.TemporaryDirectory(prefix=f"openff_bulk_{allomorph}_") as temporary:
        directory = Path(temporary)
        openmm_values = _openmm_energies(interchange)
        gromacs_values = _gromacs_energies(interchange, directory)
        required = set(ENERGY_TOLERANCES_PER_ATOM_KJ_MOL)
        missing = sorted(required - openmm_values.keys() | required - gromacs_values.keys())
        if missing:
            raise RuntimeError(f"required energy terms missing: {missing}")
        tolerances = dict(ENERGY_TOLERANCES_PER_ATOM_KJ_MOL)
        # GROMACS 2024.2 mixed precision rounds aggregate electrostatic terms.
        # The limits remain below 0.0001 kJ mol-1 atom-1 and are recorded per
        # case rather than hidden in a total-energy comparison.
        tolerances["Electrostatics"] = 6.0e-5
        tolerances["Electrostatics 1-4"] = 1.0e-4
        comparison = {
            term: {
                "difference_per_atom_kj_mol": abs(openmm_values[term] - gromacs_values[term])
                / build.mapped.atom_count,
                "tolerance_per_atom_kj_mol": tolerances[term],
                "passed": abs(openmm_values[term] - gromacs_values[term]) / build.mapped.atom_count
                <= tolerances[term],
            }
            for term in sorted(required)
        }
        smoke = _run_gromacs(interchange, directory)
        openmm_smoke = _run_openmm(build, interchange)
    return {
        "allomorph": allomorph,
        "replication": list(replication),
        "box_edges_angstrom": [float(np.linalg.norm(vector)) for vector in crystal.basis_vectors],
        "reduced_box_diagonal_angstrom": list(np.diag(reduced_periodic_box(crystal.basis_vectors))),
        "atom_count": build.mapped.atom_count,
        "stereochemistry": build.stereochemistry,
        "charge_provenance": build.charge_provenance,
        "dispersion_correction": dispersion_correction_for_shape("bulk"),
        "energy_comparison_dispersion_correction": (
            "disabled in both engines to isolate the shared base potential; "
            "generated production inputs enable it"
        ),
        "energy_comparison": comparison,
        "gromacs_smoke_test": smoke,
        "openmm_smoke_test": openmm_smoke,
        "passed": build.stereochemistry["passed"]
        and all(item["passed"] for item in comparison.values())
        and smoke["passed"]
        and openmm_smoke["passed"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=ALLOMORPHS, action="append")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation/results/openff_periodic_bulk_nagl_matrix.json",
    )
    args = parser.parse_args()
    names = tuple(args.case or ALLOMORPHS)
    cases = {name: validate_case(name) for name in names}
    result = {
        "matrix": "periodic cellulose bulk crystals, OpenFF Sage 2.3/NAGL",
        "validation_scope": (
            "model construction validated; force-field accuracy for cellulose not assessed"
        ),
        "software_versions": {
            package: metadata.version(package)
            for package in (
                "openff-toolkit",
                "openff-interchange",
                "openff-nagl",
                "openff-nagl-models",
                "openmm",
            )
        },
        "openmm_platform": {"name": "CPU", "threads": 1},
        "pme_settings": {"ewald_rtol": 1e-6, "pme_order": 8, "fourierspacing_nm": 0.05},
        "cases": cases,
        "passed": all(case["passed"] for case in cases.values()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
