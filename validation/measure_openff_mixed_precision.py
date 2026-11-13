"""Measure the mixed-precision origin of OpenMM/GROMACS bonded-energy gaps.

The OpenFF validators accept each bonded energy term within
``max(1e-5, 4 eps32 sqrt(N))`` relative to the two engines' values, where N
is the number of summed interactions.  This script records the evidence
behind that empirical acceptance criterion for the 36-chain I-alpha finite
fibril, the largest validated case:

1. OpenMM CPU and double-precision Reference platforms agree for bonded
   terms, so the gap is on the GROMACS (mixed-precision) side.
2. Rounding the coordinates to single precision does not reproduce the gap.
3. Summing the double-precision OpenMM per-interaction angle and torsion
   energies sequentially in single precision does reproduce it.

The criterion is empirical, not a general floating-point bound: rounding
error also depends on term magnitudes, ordering and cancellation.  If a case
ever exceeds it, compare against a double-precision GROMACS build.
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import openmm
from openff.units import unit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "validation"))

from _provenance import provenance  # noqa: E402

OPENFF_PACKAGES = (
    "openff-toolkit",
    "openff-interchange",
    "openff-nagl",
    "openff-nagl-models",
    "openmm",
)

import validate_openff_periodic_bulk as bulk  # noqa: E402
from validate_openff_fibrils import CASES, _crystal  # noqa: E402

from cellulose_forge.openff_sage import build_sage_system  # noqa: E402

CASE = "i_alpha_36_finite"
EPSILON32 = 2.0**-24


def _openmm(interchange: object, platform: str, periodic: bool) -> dict[str, float]:
    system = interchange.to_openmm_system(  # type: ignore[attr-defined]
        combine_nonbonded_forces=False, add_constrained_forces=True, ewald_tolerance=1.0e-6
    )
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.setUseDispersionCorrection(False)
            force.setExceptionsUsePeriodicBoundaryConditions(periodic)
        elif isinstance(force, openmm.CustomNonbondedForce):
            force.setUseLongRangeCorrection(False)
        elif periodic and hasattr(force, "setUsesPeriodicBoundaryConditions"):
            force.setUsesPeriodicBoundaryConditions(True)
    raw = bulk._get_openmm_energies(
        system=system,
        box_vectors=interchange.box.to_openmm(),  # type: ignore[attr-defined]
        positions=interchange.positions.to_openmm(),  # type: ignore[attr-defined]
        round_positions=None,
        platform=platform,
    )
    return bulk.energies(
        bulk._process(raw, system=system, combine_nonbonded_forces=False, detailed=True)
    )


def _sequential_float32(values: np.ndarray) -> float:
    total = np.float32(0.0)
    for value in values.astype(np.float32):
        total = np.float32(total + value)
    return float(total)


def _angle_terms(system: openmm.System, xyz: np.ndarray) -> np.ndarray:
    force = next(f for f in system.getForces() if isinstance(f, openmm.HarmonicAngleForce))
    rows = [force.getAngleParameters(index) for index in range(force.getNumAngles())]
    i, j, k = (np.array([row[n] for row in rows]) for n in range(3))
    theta0 = np.array([row[3].value_in_unit(openmm.unit.radian) for row in rows])
    per_radian2 = openmm.unit.kilojoule_per_mole / openmm.unit.radian**2
    constant = np.array([row[4].value_in_unit(per_radian2) for row in rows])
    left, right = xyz[i] - xyz[j], xyz[k] - xyz[j]
    theta = np.arctan2(
        np.linalg.norm(np.cross(left, right), axis=1), np.sum(left * right, axis=1)
    )
    return np.asarray(0.5 * constant * (theta - theta0) ** 2)


def _torsion_terms(system: openmm.System, xyz: np.ndarray) -> np.ndarray:
    values = []
    for force in system.getForces():
        if not isinstance(force, openmm.PeriodicTorsionForce):
            continue
        rows = [force.getTorsionParameters(index) for index in range(force.getNumTorsions())]
        p0, p1, p2, p3 = (xyz[np.array([row[n] for row in rows])] for n in range(4))
        periodicity = np.array([row[4] for row in rows])
        phase = np.array([row[5].value_in_unit(openmm.unit.radian) for row in rows])
        constant = np.array([row[6].value_in_unit(openmm.unit.kilojoule_per_mole) for row in rows])
        b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
        n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
        m1 = np.cross(n1, b1 / np.linalg.norm(b1, axis=1)[:, None])
        phi = np.arctan2(np.sum(m1 * n2, axis=1), np.sum(n1 * n2, axis=1))
        values.append(constant * (1.0 + np.cos(periodicity * phi - phase)))
    return np.concatenate(values)


def measure() -> dict[str, object]:
    allomorph, model, profile, boundary = CASES[CASE]
    periodic = boundary == "periodic"
    build = build_sage_system(_crystal(allomorph, model, profile, boundary), boundary)
    interchange = bulk._interchange(build)
    reference = _openmm(interchange, "Reference", periodic)
    cpu = _openmm(interchange, "CPU", periodic)
    with tempfile.TemporaryDirectory(prefix="openff_precision_") as temporary:
        gromacs = bulk._gromacs_energies(interchange, Path(temporary), periodic=periodic)

    xyz = np.asarray(interchange.positions.m_as(unit.nanometer), dtype=np.float64)  # type: ignore[attr-defined]
    rounded = xyz.astype(np.float32).astype(np.float64)
    interchange.positions = rounded * unit.nanometer  # type: ignore[attr-defined]
    reference_rounded = _openmm(interchange, "Reference", periodic)
    interchange.positions = xyz * unit.nanometer  # type: ignore[attr-defined]

    system = interchange.to_openmm_system(  # type: ignore[attr-defined]
        combine_nonbonded_forces=False, add_constrained_forces=True
    )
    accumulation = {}
    for term, values in (
        ("Angle", _angle_terms(system, xyz)),
        ("Torsion", _torsion_terms(system, xyz)),
    ):
        exact = float(np.sum(values))
        single = _sequential_float32(values)
        accumulation[term] = {
            "interactions": len(values),
            "double_precision_sum_kj_mol": exact,
            "sequential_single_precision_relative_error": abs(single - exact) / abs(exact),
            "error_in_units_of_eps32_sqrt_n": abs(single - exact)
            / abs(exact)
            / (EPSILON32 * np.sqrt(len(values))),
        }

    def relative(left: dict[str, float], right: dict[str, float], term: str) -> float:
        return abs(left[term] - right[term]) / max(abs(left[term]), abs(right[term]), 1.0)

    terms = ("Bond", "Angle", "Torsion", "vdW 1-4", "Electrostatics 1-4")
    return {
        "case": CASE,
        "provenance": provenance(OPENFF_PACKAGES),
        "atom_count": build.mapped.atom_count,
        "relative_differences": {
            term: {
                "openmm_cpu_vs_reference": relative(cpu, reference, term),
                "openmm_reference_vs_gromacs": relative(reference, gromacs, term),
                "openmm_reference_float32_coordinates_vs_gromacs": relative(
                    reference_rounded, gromacs, term
                ),
            }
            for term in terms
        },
        "sequential_single_precision_accumulation": accumulation,
        "acceptance_criterion": (
            "empirical: relative difference <= max(1e-5, 4 eps32 sqrt(N)); "
            "escalate to a double-precision GROMACS comparison if exceeded"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation/results/openff_mixed_precision_accumulation.json",
    )
    args = parser.parse_args()
    result = measure()
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
