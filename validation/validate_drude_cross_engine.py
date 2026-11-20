#!/usr/bin/env python3
"""Compare finite CHARMM Drude cellulose energies in NAMD and OpenMM.

The validator evaluates identical native PSF parameters, full-precision binary
coordinates and box vectors.  Its validation-only settings disable switching
and long-range Lennard-Jones corrections in both engines; production inputs
retain their documented settings.
"""

from __future__ import annotations

import argparse
import json
import struct
import subprocess
from pathlib import Path
from typing import Any, cast

from cellulose_forge.cli import main as cellulose_forge_main
from cellulose_forge.energy_drivers.namd import parse_namd_energy_log

KJ_PER_KCAL = 4.184
VECTOR_COMPONENTS = 3
NAMD_BINARY_HEADER_BYTES = 4
PSF_KIND_FIELD_COUNT = 9
PSF_KIND_INDEX = 8
PSF_DRUDE_KIND = -2
CASES = {
    "i_beta_110": ("I-beta", (1, 1, 0)),
    "i_beta_100": ("I-beta", (1, 0, 0)),
    "cellulose_ii_110": ("II", (1, 1, 0)),
}


def _build(case: Path, toppar: Path, allomorph: str, face: tuple[int, int, int]) -> dict[str, Any]:
    arguments = [
        "build",
        "--allomorph",
        allomorph,
        "--shape",
        "slab",
        "--face",
        *(str(value) for value in face),
        "--layers",
        "3",
        "--chains-per-layer",
        "6",
        "--vacuum",
        "30",
        "--force-field",
        "charmm-drude",
        "--md-engine",
        "openmm",
        "--chain-boundary",
        "finite",
        "--degree-polymerization",
        "6",
        "--padding-longitudinal",
        "20",
        "--drude-toppar-directory",
        str(toppar),
        "--output",
        str(case),
    ]
    if cellulose_forge_main(arguments) != 0:
        raise RuntimeError(f"CelluloseForge failed to build {case.name}")
    return cast(
        dict[str, Any],
        json.loads((case / "build_report.json").read_text(encoding="utf-8")),
    )


def _box_from_runner(path: Path) -> list[list[float]]:
    prefix = "BOX_NM = "
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(prefix):
            value = json.loads(line[len(prefix) :])
            if len(value) == VECTOR_COMPONENTS and all(
                len(vector) == VECTOR_COMPONENTS for vector in value
            ):
                return cast(list[list[float]], value)
    raise ValueError("OpenMM Drude runner contains no valid BOX_NM")


def _namd_configuration(box_nm: list[list[float]]) -> str:
    box_a = [[10.0 * value for value in vector] for vector in box_nm]
    center = [0.5 * sum(vector[axis] for vector in box_a) for axis in range(3)]
    vector_lines = [
        f"cellBasisVector{index} {vector[0]:.12f} {vector[1]:.12f} {vector[2]:.12f}"
        for index, vector in enumerate(box_a, 1)
    ]
    return "\n".join(
        [
            "structure system.psf",
            "coordinates structure_namd.pdb",
            "bincoordinates structure_namd.coor",
            "parameters toppar_drude_main_protein_2023a.str",
            "parameters toppar_drude_carbohydrate_2023a.str",
            "paraTypeCharmm on",
            "exclude scaled1-4",
            "oneFourScaling 1.0",
            "switching off",
            "cutoff 12.0",
            "pairlistdist 14.0",
            "LJcorrection off",
            "PME yes",
            "PMEGridSpacing 1.0",
            "PMETolerance 1.0e-6",
            *vector_lines,
            f"cellOrigin {center[0]:.12f} {center[1]:.12f} {center[2]:.12f}",
            "wrapAll off",
            "timestep 1.0",
            "rigidBonds water",
            "useSettle on",
            "mergeCrossterms no",
            "drude on",
            "drudeTemp 1.0",
            "drudeDamping 20.0",
            # The hard wall can reposition run-0 Drudes before the energy is
            # printed; it is a dynamics safeguard, not part of a single point.
            "drudeHardWall off",
            "drudeBondLen 0.25",
            "drudeNbTholeCut 5.0",
            "langevin on",
            "langevinTemp 300",
            "langevinDamping 1.0",
            "langevinHydrogen on",
            # A run-0 energy evaluation must not randomize Drude coordinates.
            "temperature 0",
            "outputEnergies 1",
            "outputName single_point",
            "run 0",
            "",
        ]
    )


_OPENMM_EVALUATOR = r"""
import json, math, struct, sys
import openmm
from openmm import unit
from openmm.app import CharmmParameterSet, CharmmPsfFile, PME

box_nm = json.loads(sys.argv[1])
psf = CharmmPsfFile("system.psf")
parameters = CharmmParameterSet(
    "toppar_drude_main_protein_2023a.str", "toppar_drude_carbohydrate_2023a.str"
)
with open("structure_namd.coor", "rb") as handle:
    count = struct.unpack("=i", handle.read(4))[0]
    values = struct.unpack(f"={count * 3}d", handle.read())
positions = [openmm.Vec3(*values[index:index+3]) for index in range(0, len(values), 3)]
positions = positions * unit.angstrom
vectors = tuple(openmm.Vec3(*row) * unit.nanometer for row in box_nm)
lengths = [math.sqrt(sum(value*value for value in row)) for row in box_nm]
def angle(first, second):
    cosine = sum(a*b for a, b in zip(first, second)) / (
        math.sqrt(sum(value*value for value in first))
        * math.sqrt(sum(value*value for value in second))
    )
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))
angles = [angle(box_nm[1], box_nm[2]), angle(box_nm[0], box_nm[2]),
          angle(box_nm[0], box_nm[1])]
psf.setBox(*(value*unit.nanometer for value in lengths),
           *(value*unit.degree for value in angles))
system = psf.createSystem(
    parameters, nonbondedMethod=PME, nonbondedCutoff=1.2*unit.nanometer,
    constraints=None, ewaldErrorTolerance=1.0e-6,
)
system.setDefaultPeriodicBoxVectors(*vectors)
for force in system.getForces():
    if isinstance(force, openmm.NonbondedForce):
        force.setUseSwitchingFunction(False)
        force.setUseDispersionCorrection(False)
for index, force in enumerate(system.getForces()):
    force.setForceGroup(index)
context = openmm.Context(
    system, openmm.VerletIntegrator(0.001),
    openmm.Platform.getPlatformByName("Reference")
)
context.setPeriodicBoxVectors(*vectors)
context.setPositions(positions)
context.computeVirtualSites()

# NAMD folds the Drude spring and anisotropy into BOND, while OpenMM's
# DrudeForce also includes screened (Thole) pairs.  Evaluate a copy without
# screened pairs so the physically equivalent spring terms can be compared.
spring_system = openmm.System()
for particle in range(system.getNumParticles()):
    spring_system.addParticle(system.getParticleMass(particle))
spring_system.setDefaultPeriodicBoxVectors(*vectors)
spring_force = openmm.DrudeForce()
for force in system.getForces():
    if isinstance(force, openmm.DrudeForce):
        for particle in range(force.getNumParticles()):
            spring_force.addParticle(*force.getParticleParameters(particle))
spring_system.addForce(spring_force)
spring_context = openmm.Context(
    spring_system, openmm.VerletIntegrator(0.001),
    openmm.Platform.getPlatformByName("Reference")
)
spring_context.setPeriodicBoxVectors(*vectors)
spring_context.setPositions(context.getState(getPositions=True).getPositions())
spring_anisotropy = spring_context.getState(getEnergy=True).getPotentialEnergy()
spring_anisotropy = spring_anisotropy.value_in_unit(unit.kilojoule_per_mole)

covalent = {
    tuple(sorted((bond.atom1.idx, bond.atom2.idx))) for bond in psf.bond_list
}
terms = {name: 0.0 for name in (
    "bond", "angle", "dihedral", "improper", "cross", "nonbonded", "drude"
)}
forces = []
for index, force in enumerate(system.getForces()):
    energy = context.getState(getEnergy=True, groups={index}).getPotentialEnergy()
    value = energy.value_in_unit(unit.kilojoule_per_mole)
    kind = force.__class__.__name__
    forces.append({"index": index, "class": kind, "kj_mol": value})
    if kind == "HarmonicBondForce":
        pairs = {
            tuple(sorted(map(int, force.getBondParameters(item)[:2])))
            for item in range(force.getNumBonds())
        }
        target = "bond" if pairs <= covalent else "angle"
        terms[target] += value
    elif kind == "HarmonicAngleForce": terms["angle"] += value
    elif kind == "PeriodicTorsionForce": terms["dihedral"] += value
    elif kind == "CustomTorsionForce": terms["improper"] += value
    elif kind == "CMAPTorsionForce": terms["cross"] += value
    elif kind == "NonbondedForce": terms["nonbonded"] += value
    elif kind == "DrudeForce": terms["drude"] += value
    elif kind != "CMMotionRemover" and value != 0:
        raise SystemExit(f"unsupported nonzero OpenMM force: {kind} {value}")
potential = context.getState(getEnergy=True).getPotentialEnergy()
terms["potential"] = potential.value_in_unit(unit.kilojoule_per_mole)
print(json.dumps({
    "terms_kj_mol": terms,
    "forces": forces,
    "drude_spring_anisotropy_kj_mol": spring_anisotropy,
}, sort_keys=True))
"""


def _drude_indices(psf: Path) -> list[int]:
    """Return zero-based Drude-particle indices from a native extended PSF."""
    lines = psf.read_text(encoding="utf-8").splitlines()
    atom_header = next((index for index, line in enumerate(lines) if "!NATOM" in line), None)
    if atom_header is None:
        raise ValueError("Drude PSF contains no !NATOM section")
    atom_count = int(lines[atom_header].split()[0])
    atom_lines = lines[atom_header + 1 : atom_header + 1 + atom_count]
    if len(atom_lines) != atom_count:
        raise ValueError("Drude PSF atom section is truncated")
    indices = [
        int(fields[0]) - 1
        for fields in (line.split() for line in atom_lines)
        if len(fields) >= PSF_KIND_FIELD_COUNT and int(fields[PSF_KIND_INDEX]) == PSF_DRUDE_KIND
    ]
    if not indices:
        raise ValueError("Drude PSF contains no Drude particles")
    return indices


def _displace_drudes(coordinates: Path, indices: list[int], distance_a: float = 0.01) -> bytes:
    """Displace Drudes deterministically and return the original binary frame."""
    original = coordinates.read_bytes()
    if len(original) < NAMD_BINARY_HEADER_BYTES:
        raise ValueError("NAMD binary coordinate file is truncated")
    particle_count = struct.unpack("=i", original[:NAMD_BINARY_HEADER_BYTES])[0]
    expected = NAMD_BINARY_HEADER_BYTES + particle_count * 3 * 8
    if len(original) != expected:
        raise ValueError("NAMD binary coordinate file has an unexpected size")
    values = list(struct.unpack(f"={particle_count * 3}d", original[NAMD_BINARY_HEADER_BYTES:]))
    for serial, particle in enumerate(indices):
        if not 0 <= particle < particle_count:
            raise ValueError("Drude particle index lies outside the coordinate frame")
        values[3 * particle + serial % 3] += distance_a
    coordinates.write_bytes(
        struct.pack("=i", particle_count) + struct.pack(f"={particle_count * 3}d", *values)
    )
    return original


def _run_engines(
    directory: Path, *, state: str, namd: Path, openmm_python: Path, box_nm: list[list[float]]
) -> tuple[dict[str, float], dict[str, Any]]:
    namd_run = subprocess.run(
        [str(namd), "+p2", "single_point_namd.conf"],
        cwd=directory,
        capture_output=True,
        text=True,
        check=False,
    )
    namd_log = directory / f"single_point_namd_{state}.log"
    namd_log.write_text(namd_run.stdout + namd_run.stderr, encoding="utf-8")
    if namd_run.returncode:
        raise RuntimeError(f"NAMD failed for {directory.name}/{state}: {namd_run.returncode}")
    openmm_run = subprocess.run(
        [str(openmm_python), "single_point_openmm.py", json.dumps(box_nm)],
        cwd=directory,
        capture_output=True,
        text=True,
        check=False,
    )
    if openmm_run.returncode:
        raise RuntimeError(f"OpenMM failed for {directory.name}/{state}: {openmm_run.stderr}")
    openmm = cast(dict[str, Any], json.loads(openmm_run.stdout))
    raw_namd = parse_namd_energy_log(namd_log)
    namd_terms = {
        "bond": raw_namd["BOND"] * KJ_PER_KCAL,
        "angle": raw_namd["ANGLE"] * KJ_PER_KCAL,
        "dihedral": raw_namd["DIHED"] * KJ_PER_KCAL,
        "improper": raw_namd["IMPRP"] * KJ_PER_KCAL,
        "cross": raw_namd.get("CROSS", 0.0) * KJ_PER_KCAL,
        "nonbonded": (raw_namd["ELECT"] + raw_namd["VDW"]) * KJ_PER_KCAL,
        "potential": raw_namd["POTENTIAL"] * KJ_PER_KCAL,
    }
    return namd_terms, openmm


def _comparison(
    namd_value: float, openmm_value: float, *, relative_tolerance: float = 1.0e-5
) -> dict[str, float | bool]:
    difference = openmm_value - namd_value
    scale = max(abs(openmm_value), abs(namd_value), 1.0)
    tolerance = max(0.05, relative_tolerance * scale)
    return {
        "namd_kj_mol": namd_value,
        "openmm_kj_mol": openmm_value,
        "difference_kj_mol": difference,
        "relative_difference": abs(difference) / scale,
        "tolerance_kj_mol": tolerance,
        "passed": abs(difference) <= tolerance,
    }


def _compare_displaced_state(
    namd_terms: dict[str, float], openmm: dict[str, Any]
) -> dict[str, Any]:
    openmm_terms = cast(dict[str, float], openmm["terms_kj_mol"])
    comparisons = {
        name: _comparison(namd_terms[name], openmm_terms[name])
        for name in ("angle", "dihedral", "improper", "cross")
    }
    spring_anisotropy = float(openmm["drude_spring_anisotropy_kj_mol"])
    thole = openmm_terms["drude"] - spring_anisotropy
    comparisons["bond_plus_spring_anisotropy"] = _comparison(
        namd_terms["bond"],
        openmm_terms["bond"] + spring_anisotropy,
        relative_tolerance=5.0e-5,
    )
    comparisons["nonbonded_plus_thole"] = _comparison(
        namd_terms["nonbonded"],
        openmm_terms["nonbonded"] + thole,
        relative_tolerance=1.0e-4,
    )
    comparisons["potential"] = _comparison(namd_terms["potential"], openmm_terms["potential"])
    return {
        "comparisons": comparisons,
        "openmm_drude_force_kj_mol": openmm_terms["drude"],
        "openmm_drude_spring_anisotropy_kj_mol": spring_anisotropy,
        "openmm_thole_kj_mol": thole,
        "openmm_forces": openmm["forces"],
        "passed": all(item["passed"] for item in comparisons.values())
        and abs(openmm_terms["drude"]) > 1.0,
    }


def _compare_polarization_response(
    first_namd: dict[str, float],
    first_openmm: dict[str, Any],
    second_namd: dict[str, float],
    second_openmm: dict[str, Any],
) -> dict[str, Any]:
    first_openmm_terms = cast(dict[str, float], first_openmm["terms_kj_mol"])
    second_openmm_terms = cast(dict[str, float], second_openmm["terms_kj_mol"])
    comparisons = {
        "spring_and_anisotropy": _comparison(
            second_namd["bond"] - first_namd["bond"],
            float(
                second_openmm["drude_spring_anisotropy_kj_mol"]
                - first_openmm["drude_spring_anisotropy_kj_mol"]
            ),
            relative_tolerance=5.0e-5,
        )
    }
    namd_response = (
        second_namd["bond"]
        + second_namd["nonbonded"]
        - first_namd["bond"]
        - first_namd["nonbonded"]
    )
    openmm_response = (
        second_openmm_terms["bond"]
        + second_openmm_terms["nonbonded"]
        + second_openmm_terms["drude"]
        - first_openmm_terms["bond"]
        - first_openmm_terms["nonbonded"]
        - first_openmm_terms["drude"]
    )
    comparisons["total_polarization"] = _comparison(
        namd_response,
        openmm_response,
        relative_tolerance=1.0e-3,
    )
    comparisons["potential"] = _comparison(
        second_namd["potential"] - first_namd["potential"],
        second_openmm_terms["potential"] - first_openmm_terms["potential"],
        relative_tolerance=1.0e-3,
    )
    return {
        "comparisons": comparisons,
        "gap_status": "uncharacterized",
        "passed": all(item["passed"] for item in comparisons.values()),
    }


def _run_case(  # noqa: PLR0913 - explicit executables and scientific case remain visible
    directory: Path,
    *,
    namd: Path,
    openmm_python: Path,
    toppar: Path,
    allomorph: str,
    face: tuple[int, int, int],
) -> dict[str, Any]:
    directory.mkdir(parents=True)
    report = _build(directory, toppar, allomorph, face)
    box_nm = _box_from_runner(directory / "run_drude_openmm.py")
    (directory / "single_point_namd.conf").write_text(_namd_configuration(box_nm), encoding="utf-8")
    evaluator = directory / "single_point_openmm.py"
    evaluator.write_text(_OPENMM_EVALUATOR, encoding="utf-8")
    coordinates = directory / "structure_namd.coor"
    drude_indices = _drude_indices(directory / "system.psf")
    original = _displace_drudes(coordinates, drude_indices, distance_a=0.01)
    try:
        first_namd, first_openmm = _run_engines(
            directory,
            state="drudes_0p01_angstrom",
            namd=namd,
            openmm_python=openmm_python,
            box_nm=box_nm,
        )
        coordinates.write_bytes(original)
        _displace_drudes(coordinates, drude_indices, distance_a=0.02)
        second_namd, second_openmm = _run_engines(
            directory,
            state="drudes_0p02_angstrom",
            namd=namd,
            openmm_python=openmm_python,
            box_nm=box_nm,
        )
    finally:
        coordinates.write_bytes(original)
    first = _compare_displaced_state(first_namd, first_openmm)
    second = _compare_displaced_state(second_namd, second_openmm)
    response = _compare_polarization_response(
        first_namd,
        first_openmm,
        second_namd,
        second_openmm,
    )
    return {
        "allomorph": allomorph,
        "face": list(face),
        "particle_count": report["topology_counts"]["particles"],
        "validation_settings": {
            "switching": "off in both engines",
            "lj_tail_correction": "off in both engines",
            "cutoff_nm": 1.2,
            "pme_tolerance": 1.0e-6,
            "coordinates": "shared full-precision NAMD binary coordinates",
            "drude_displacements_angstrom": [0.01, 0.02],
            "nonbonded_plus_thole_relative_tolerance": 1.0e-4,
            "term_relative_tolerance": 1.0e-5,
            "spring_response_relative_tolerance": 5.0e-5,
            "polarization_response_relative_tolerance": 1.0e-3,
            "absolute_tolerance_kj_mol": 0.05,
            "polarization_response_gap": "uncharacterized",
        },
        "states": {
            "drudes_0p01_angstrom": first,
            "drudes_0p02_angstrom": second,
        },
        "response_0p01_to_0p02_angstrom": response,
        "passed": first["passed"] and second["passed"] and response["passed"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--drude-toppar-directory", type=Path, required=True)
    parser.add_argument("--namd", type=Path, required=True)
    parser.add_argument("--openmm-python", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.work_directory.exists() and any(args.work_directory.iterdir()):
        parser.error("--work-directory must be absent or empty")
    args.work_directory.mkdir(parents=True, exist_ok=True)
    results = {
        name: _run_case(
            args.work_directory / name,
            namd=args.namd,
            openmm_python=args.openmm_python,
            toppar=args.drude_toppar_directory,
            allomorph=allomorph,
            face=face,
        )
        for name, (allomorph, face) in CASES.items()
    }
    summary = {
        "schema": 1,
        "description": "finite CHARMM Drude NAMD/OpenMM matched single points",
        "cases": results,
        "passed": all(case["passed"] for case in results.values()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
