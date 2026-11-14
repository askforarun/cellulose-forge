"""Measure OpenMM Reference vs CPU nonbonded energies in tilted periodic boxes.

The OPLS OpenMM builder evaluates its build-time energies on the CPU platform
when the third box vector is tilted (cellulose I-alpha bulk crystals).  This
script records the evidence for that choice.  For each periodic OPLS-AA bulk
crystal it evaluates, with identical coordinates and box:

* the geometric-mixing Lennard-Jones ``CustomNonbondedForce`` (1.0 nm plain
  cutoff, intramolecular 1-2/1-3/1-4 exclusions) on the Reference and CPU
  platforms, with coordinates as built and wrapped into the reduced cell;
* the same energy as an exact NumPy sum over the nearest periodic image of
  every pair (valid because the cutoff is below half of every box height);
* the PME ``NonbondedForce`` Coulomb energy on the Reference and CPU platforms.

Run with an interpreter that provides OpenMM::

    PYTHONPATH=src python validation/measure_openmm_tilted_box.py
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "validation"))

import openmm  # noqa: E402
from _provenance import provenance  # noqa: E402
from openmm import unit  # noqa: E402

from cellulose_forge.crystallography import load_allomorph  # noqa: E402
from cellulose_forge.geometry import build_bulk_crystal  # noqa: E402
from cellulose_forge.opls import build_opls_system  # noqa: E402
from cellulose_forge.opls_engines import chain_template, reduced_box_angstrom  # noqa: E402

CUTOFF_NM = 1.0
LJ = (
    "4*epsilon*((sigma/r)^12-(sigma/r)^6); "
    "sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)"
)
CASES = {
    "i_alpha_3x4x4": ("i_alpha", (3, 4, 4)),
    "i_alpha_3x4x5": ("i_alpha", (3, 4, 5)),
    "i_beta_3x3x3": ("i_beta", (3, 3, 3)),
    "ii_3x3x3": ("ii", (3, 3, 3)),
    "iii_i_5x3x3": ("iii_i", (5, 3, 3)),
}
# A "discrepancy" is a Reference/exact difference far above both the CPU
# single-precision level and the agreement of the untilted controls.
DISCREPANCY_RELATIVE = 1.0e-4
TILT_NM = 1.0e-9


def _system(case: tuple[str, tuple[int, int, int]]) -> dict[str, object]:
    allomorph, cells = case
    build = build_opls_system(build_bulk_crystal(load_allomorph(allomorph), cells), "opls-aa",
                              "periodic")
    template = chain_template(build)
    box = reduced_box_angstrom(build) / 10.0
    positions = np.asarray(build.mapped.cartesian_coordinates, dtype=np.float64) / 10.0
    return {"template": template, "chains": build.mapped.chain_count, "box": box,
            "positions": positions}


def _forces(data: dict[str, object]) -> tuple[openmm.System, openmm.System]:
    template = data["template"]
    per_chain = len(template["masses"])  # type: ignore[index]
    vectors = [openmm.Vec3(*row) for row in data["box"]]  # type: ignore[attr-defined]
    systems = []
    for kind in ("lj", "coulomb"):
        system = openmm.System()
        system.setDefaultPeriodicBoxVectors(*vectors)
        if kind == "lj":
            force = openmm.CustomNonbondedForce(LJ)
            force.addPerParticleParameter("sigma")
            force.addPerParticleParameter("epsilon")
            force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
        else:
            force = openmm.NonbondedForce()
            force.setNonbondedMethod(openmm.NonbondedForce.PME)
            force.setEwaldErrorTolerance(1.0e-6)
        force.setCutoffDistance(CUTOFF_NM)
        for chain in range(data["chains"]):  # type: ignore[call-overload]
            offset = chain * per_chain
            for mass, charge, (sigma, epsilon) in zip(
                template["masses"], template["charges"], template["lj"],  # type: ignore[index]
                strict=True,
            ):
                system.addParticle(mass)
                if kind == "lj":
                    force.addParticle([sigma, epsilon])
                else:
                    force.addParticle(charge, 1.0, 0.0)
            for a, b, *_ in template["exceptions"]:  # type: ignore[index]
                if kind == "lj":
                    force.addExclusion(offset + a, offset + b)
                else:
                    force.addException(offset + a, offset + b, 0.0, 1.0, 0.0)
        system.addForce(force)
        systems.append(system)
    return systems[0], systems[1]


def _energy(system: openmm.System, platform: str, positions: np.ndarray) -> float:
    context = openmm.Context(system, openmm.VerletIntegrator(0.001),
                             openmm.Platform.getPlatformByName(platform))
    context.setPositions([openmm.Vec3(*row) for row in positions])
    energy = context.getState(getEnergy=True).getPotentialEnergy()
    del context
    return float(energy.value_in_unit(unit.kilojoule_per_mole))


def _exact_lj(data: dict[str, object]) -> float:
    template = data["template"]
    per_chain = len(template["masses"])  # type: ignore[index]
    excluded = {(a, b) for a, b, *_ in template["exceptions"]}  # type: ignore[index]
    excluded |= {(b, a) for a, b in excluded}
    lj = np.asarray(template["lj"] * data["chains"])  # type: ignore[index, operator]
    positions = data["positions"]
    box = data["box"]
    shifts = np.array([i * box[0] + j * box[1] + k * box[2]  # type: ignore[index]
                       for i, j, k in itertools.product((-2, -1, 0, 1, 2), repeat=3)])
    count = len(positions)  # type: ignore[arg-type]
    total = 0.0
    for i in range(count - 1):
        delta = positions[i + 1:] - positions[i]  # type: ignore[index]
        distance = np.linalg.norm(delta[:, None, :] + shifts[None, :, :], axis=2).min(axis=1)
        partners = np.arange(i + 1, count)
        keep = distance < CUTOFF_NM
        same = partners // per_chain == i // per_chain
        for index in np.where(same & keep)[0]:
            if (i % per_chain, partners[index] % per_chain) in excluded:
                keep[index] = False
        sigma = np.sqrt(lj[i, 0] * lj[partners[keep], 0])
        epsilon = np.sqrt(lj[i, 1] * lj[partners[keep], 1])
        ratio6 = (sigma / distance[keep]) ** 6
        total += float(np.sum(4.0 * epsilon * (ratio6 * ratio6 - ratio6)))
    return total


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "validation/results/openmm_tilted_box.json")
    args = parser.parse_args()
    cases = {}
    for name, case in CASES.items():
        data = _system(case)
        box = data["box"]
        lj_system, coulomb_system = _forces(data)
        positions = data["positions"]
        fractional = positions @ np.linalg.inv(box)  # type: ignore[operator]
        wrapped = (fractional - np.floor(fractional)) @ box
        exact = _exact_lj(data)
        lj = {
            "reference": _energy(lj_system, "Reference", positions),  # type: ignore[arg-type]
            "reference_wrapped": _energy(lj_system, "Reference", wrapped),
            "cpu": _energy(lj_system, "CPU", positions),  # type: ignore[arg-type]
        }
        coulomb = {
            "reference": _energy(coulomb_system, "Reference", positions),  # type: ignore[arg-type]
            "cpu": _energy(coulomb_system, "CPU", positions),  # type: ignore[arg-type]
        }
        tilted = bool(abs(box[2][0]) > TILT_NM or abs(box[2][1]) > TILT_NM)  # type: ignore[index]
        cases[name] = {
            "box_nm": np.asarray(box).tolist(),
            "third_vector_tilted": tilted,
            "atoms": len(positions),  # type: ignore[arg-type]
            "lennard_jones_kj_mol": {**lj, "exact_nearest_image": exact},
            "lennard_jones_relative_to_exact": {
                key: abs(value - exact) / abs(exact) for key, value in lj.items()
            },
            "pme_coulomb_kj_mol": coulomb,
            "pme_coulomb_reference_vs_cpu_relative": abs(coulomb["reference"] - coulomb["cpu"])
            / abs(coulomb["cpu"]),
        }
        cases[name]["reference_discrepancy"] = (
            cases[name]["lennard_jones_relative_to_exact"]["reference"] > DISCREPANCY_RELATIVE
        )
        print(name, json.dumps(cases[name]["lennard_jones_relative_to_exact"]), flush=True)
    result = {
        "purpose": "evidence for using the OpenMM CPU platform for tilted-third-vector boxes",
        "provenance": provenance(("numpy", "openmm")),
        "openmm_version": openmm.__version__,
        "cutoff_nm": CUTOFF_NM,
        "discrepancy_relative_threshold": DISCREPANCY_RELATIVE,
        "cases": cases,
        "reference_discrepancy_only_in_tilted_boxes": all(
            case["reference_discrepancy"] == case["third_vector_tilted"]
            for case in cases.values()
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
