"""Validate the exact experimental OpenFF files delivered to users."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import openmm
from openmm import app, unit

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

from cellulose_forge.crystallography import MillerIndex, load_allomorph  # noqa: E402
from cellulose_forge.geometry import (  # noqa: E402
    build_36_chain_fibril,
    build_bulk_crystal,
    build_chain_parallel_slab,
)
from cellulose_forge.openff_sage import (  # noqa: E402
    build_sage_system,
    dispersion_correction_for_shape,
    write_sage_outputs,
)

# GROMACS EnerPres and OpenMM's exact long-range correction for a switched
# Lennard-Jones potential use different averaging, giving ~0.4% differences.
LRC_RELATIVE_TOLERANCE = 5.0e-3
# Shipped files use production PME accuracy (OpenMM error tolerance 5e-4,
# GROMACS default PME order 4), so totals agree less tightly than the
# term-by-term validators, which use 1e-6 Ewald tolerance and PME order 8.
TOTAL_RELATIVE_TOLERANCE = 5.0e-4


def _run(command: list[str], directory: Path, *, stdin: str | None = None) -> str:
    result = subprocess.run(
        command,
        cwd=directory,
        input=stdin,
        text=True,
        capture_output=True,
        check=True,
    )
    return result.stdout + result.stderr


def _openmm_energy(system: openmm.System, pdbx: app.PDBxFile) -> float:
    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(
        system,
        integrator,
        openmm.Platform.getPlatformByName("CPU"),
        {"Threads": "1"},
    )
    context.setPositions(pdbx.positions)
    context.setPeriodicBoxVectors(*pdbx.topology.getPeriodicBoxVectors())
    energy = (
        context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
    )
    del context, integrator
    return float(energy)


def _dispersion_settings(system: openmm.System) -> dict[str, bool]:
    """Read the long-range dispersion flags exactly as serialized."""
    settings: dict[str, bool] = {}
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            settings["NonbondedForce"] = bool(force.getUseDispersionCorrection())
        elif isinstance(force, openmm.CustomNonbondedForce):
            settings["CustomNonbondedForce"] = bool(force.getUseLongRangeCorrection())
    if set(settings) != {"NonbondedForce", "CustomNonbondedForce"}:
        raise RuntimeError(f"unexpected shipped OpenMM nonbonded forces: {sorted(settings)}")
    return settings


def _disable_dispersion(system: openmm.System) -> None:
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.setUseDispersionCorrection(False)
        elif isinstance(force, openmm.CustomNonbondedForce):
            force.setUseLongRangeCorrection(False)


def _xvg_values(path: Path) -> list[float]:
    rows = [
        line.split()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith(("#", "@"))
    ]
    if len(rows) != 1:
        raise RuntimeError(f"expected one GROMACS energy row in {path}")
    return [float(value) for value in rows[0][1:]]


def _gromacs_shipped_checks(directory: Path, *, dispersion: bool) -> dict[str, object]:
    logs = [
        _run(
            [
                "gmx",
                "grompp",
                "-f",
                "minim.mdp",
                "-c",
                "structure.gro",
                "-p",
                "system.top",
                "-o",
                "minim.tpr",
                "-maxwarn",
                "0",
            ],
            directory,
        ),
        _run(["gmx", "mdrun", "-deffnm", "minim", "-ntmpi", "1"], directory),
        _run(
            [
                "gmx",
                "grompp",
                "-f",
                "nvt.mdp",
                "-c",
                "minim.gro",
                "-p",
                "system.top",
                "-o",
                "nvt.tpr",
                "-maxwarn",
                "0",
            ],
            directory,
        ),
        _run(["gmx", "mdrun", "-deffnm", "nvt", "-ntmpi", "1"], directory),
    ]
    # Single point from the shipped minimization input: keep every nonbonded
    # setting exactly as shipped and change only the run control.
    point = (directory / "minim.mdp").read_text(encoding="utf-8")
    overrides = {"integrator": "md", "nsteps": "0", "constraints": "none"}
    for key, value in overrides.items():
        point, count = re.subn(
            rf"^{key}\s*=.*$", f"{key} = {value}", point, count=1, flags=re.MULTILINE
        )
        if count != 1:
            raise RuntimeError(f"shipped minim.mdp has no {key} setting")
    point = re.sub(r"^(emtol|emstep)\s*=.*\n", "", point, flags=re.MULTILINE)
    point += "continuation = yes\n"
    (directory / "point.mdp").write_text(point, encoding="utf-8")
    _run(
        [
            "gmx",
            "grompp",
            "-f",
            "point.mdp",
            "-c",
            "structure.gro",
            "-p",
            "system.top",
            "-o",
            "point.tpr",
            "-maxwarn",
            "0",
        ],
        directory,
    )
    _run(
        [
            "gmx",
            "mdrun",
            "-s",
            "point.tpr",
            "-rerun",
            "structure.gro",
            "-deffnm",
            "point",
            "-ntmpi",
            "1",
        ],
        directory,
    )
    _run(
        ["gmx", "energy", "-f", "point.edr", "-o", "potential.xvg"],
        directory,
        stdin="Potential\n0\n",
    )
    potential = _xvg_values(directory / "potential.xvg")[0]
    correction = 0.0
    if dispersion:
        _run(
            ["gmx", "energy", "-f", "point.edr", "-o", "dispersion.xvg"],
            directory,
            stdin="Disper.-corr.\n0\n",
        )
        correction = _xvg_values(directory / "dispersion.xvg")[0]
    combined = "\n".join(logs) + (directory / "nvt.log").read_text(encoding="utf-8")
    markers = [
        marker
        for marker, pattern in {
            "LINCS": r"LINCS WARNING",
            "NaN": r"\bnan\b",
            "non-finite": r"non-finite",
        }.items()
        if re.search(pattern, combined, flags=re.IGNORECASE)
    ]
    return {
        "grompp_maxwarn": 0,
        "minimization_completed": (directory / "minim.gro").is_file(),
        "nvt_reached_final_step": "Statistics over 1001 steps" in combined,
        "failure_markers": markers,
        "point_potential_kj_mol": potential,
        "dispersion_correction_kj_mol": correction,
        "passed": (directory / "minim.gro").is_file()
        and "Statistics over 1001 steps" in combined
        and not markers,
    }


def _case(label: str, crystal: object, shape: str) -> dict[str, object]:
    dispersion = dispersion_correction_for_shape(shape)
    build = build_sage_system(crystal, "periodic")
    with tempfile.TemporaryDirectory(prefix=f"openff_shipped_{label}_") as temporary:
        directory = Path(temporary)
        write_sage_outputs(
            build,
            directory,
            ("gromacs", "openmm"),
            dispersion_correction=dispersion,
        )
        system = openmm.XmlSerializer.deserialize((directory / "system.xml").read_text())
        pdbx = app.PDBxFile(str(directory / "topology.cif"))
        # Read the serialized flags and evaluate the shipped system before any
        # modification; only then remove the correction to measure its size.
        shipped_settings = _dispersion_settings(system)
        openmm_on = _openmm_energy(system, pdbx)
        _disable_dispersion(system)
        openmm_off = _openmm_energy(system, pdbx)
        openmm_correction = openmm_on - openmm_off
        gromacs = _gromacs_shipped_checks(directory, dispersion=dispersion)
    gromacs_correction = float(gromacs["dispersion_correction_kj_mol"])
    correction_relative_difference = (
        abs(openmm_correction - gromacs_correction)
        / max(abs(openmm_correction), abs(gromacs_correction), 1.0)
        if dispersion
        else 0.0
    )
    total_relative_difference = abs(openmm_on - float(gromacs["point_potential_kj_mol"])) / max(
        abs(openmm_on), abs(float(gromacs["point_potential_kj_mol"])), 1.0
    )
    correction_passed = (
        correction_relative_difference <= LRC_RELATIVE_TOLERANCE if dispersion else True
    )
    return {
        "shape": shape,
        "dispersion_correction": dispersion,
        "openmm_force_settings": shipped_settings,
        "gromacs": gromacs,
        "openmm_point_potential_kj_mol": openmm_on,
        "openmm_dispersion_correction_kj_mol": openmm_correction,
        "correction_relative_difference": correction_relative_difference,
        "correction_relative_tolerance": LRC_RELATIVE_TOLERANCE,
        "total_potential_relative_difference": total_relative_difference,
        "total_potential_relative_tolerance": TOTAL_RELATIVE_TOLERANCE,
        "passed": bool(
            all(value is dispersion for value in shipped_settings.values())
            and gromacs["passed"]
            and correction_passed
            and total_relative_difference <= TOTAL_RELATIVE_TOLERANCE
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation/results/openff_shipped_outputs.json",
    )
    args = parser.parse_args()
    cases = {
        "i_beta_110_periodic_slab": _case(
            "slab",
            build_chain_parallel_slab(
                load_allomorph("i_beta"),
                MillerIndex(1, 1, 0),
                layers=2,
                chains_per_layer=4,
                degree_polymerization=6,
                vacuum=30.0,
                chain_boundary="periodic",
            ),
            "slab",
        ),
        "i_beta_periodic_bulk": _case(
            "bulk",
            build_bulk_crystal(load_allomorph("i_beta"), (3, 3, 2)),
            "bulk",
        ),
        "i_beta_36_chain_periodic_fibril": _case(
            "fibril",
            build_36_chain_fibril(
                load_allomorph("i_beta"),
                degree_polymerization=6,
                vacuum=20.0,
                chain_boundary="periodic",
            ),
            "fibril",
        ),
    }
    result = {
        "validation": "exact shipped experimental OpenFF outputs",
        "provenance": provenance(OPENFF_PACKAGES),
        "cases": cases,
        "passed": all(bool(case["passed"]) for case in cases.values()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
