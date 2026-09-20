"""Build and smoke-test periodic GLYCAM bulk crystals for every allomorph."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

from compare_force_fields import ENERGY_TERMS, _parse_xvg, _run, _summarize_rows, _write_protocol

ROOT = Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "validation"
SYSTEMS = {
    "I-alpha": (4, 4, 4),
    "I-beta": (4, 4, 4),
    "II": (4, 4, 4),
    # III-I has a short a axis; six repeats satisfy the PME cutoff condition.
    "III-I": (6, 4, 4),
}


def _safe_name(allomorph: str) -> str:
    return allomorph.lower().replace("-", "_")


def _validate_one(
    args: argparse.Namespace,
    allomorph: str,
    cells: tuple[int, int, int],
) -> dict[str, Any]:
    directory = args.work_directory / _safe_name(allomorph)
    directory.mkdir(parents=True, exist_ok=True)
    build_command = [
        sys.executable,
        "-m",
        "cellulose_forge.cli",
        "build",
        "--allomorph",
        allomorph,
        "--shape",
        "bulk",
        "--cells",
        *(str(value) for value in cells),
        "--chain-boundary",
        "periodic",
        "--force-field",
        "glycam",
        "--tleap-executable",
        args.tleap_executable,
        "--output",
        str(directory),
    ]
    if args.parmed_python:
        build_command.extend(["--parmed-python", str(args.parmed_python)])
    _run(build_command, cwd=ROOT, log_path=directory / "build.log")
    _write_protocol(
        VALIDATION / "periodic_pme_nvt.mdp",
        directory / "nvt.mdp",
        nsteps=args.nvt_steps,
    )
    em_grompp = _run(
        [
            args.gmx,
            "grompp",
            "-f",
            VALIDATION / "periodic_pme_minim.mdp",
            "-c",
            "structure.gro",
            "-p",
            "system.top",
            "-o",
            "em.tpr",
            "-po",
            "em_out.mdp",
        ],
        cwd=directory,
        log_path=directory / "grompp_em.log",
    )
    em_output = _run(
        [args.gmx, "mdrun", "-deffnm", "em", "-nt", str(args.threads)],
        cwd=directory,
        log_path=directory / "mdrun_em.log",
    )
    _run(
        [
            args.gmx,
            "grompp",
            "-f",
            "nvt.mdp",
            "-c",
            "em.gro",
            "-p",
            "system.top",
            "-o",
            "nvt.tpr",
            "-po",
            "nvt_out.mdp",
        ],
        cwd=directory,
        log_path=directory / "grompp_nvt.log",
    )
    _run(
        [args.gmx, "mdrun", "-deffnm", "nvt", "-nt", str(args.threads)],
        cwd=directory,
        log_path=directory / "mdrun_nvt.log",
    )
    observables = directory / "nvt_observables.xvg"
    _run(
        [args.gmx, "energy", "-f", "nvt.edr", "-o", observables.name, "-xvg", "none"],
        cwd=directory,
        log_path=directory / "energy.log",
        stdin="\n".join((*ENERGY_TERMS, "0", "")),
    )
    report = json.loads((directory / "build_report.json").read_text(encoding="utf-8"))
    steps = re.search(r"converged to Fmax < 1000 in (\d+) steps", em_output)
    force = re.search(r"Maximum force\s*=\s*([-+0-9.eE]+)", em_output)
    fudge = re.search(r"Generating 1-4 interactions: fudge = ([0-9.]+)", em_grompp)
    if steps is None or force is None or fudge is None:
        msg = f"could not parse GROMACS validation output for {allomorph}"
        raise ValueError(msg)
    return {
        "allomorph": allomorph,
        "cells": list(cells),
        "atom_count": report["force_field_atom_count"],
        "chain_count": report["chain_count"],
        "degree_polymerization": report["degree_polymerization"],
        "fudge_14": float(fudge.group(1)),
        "minimization_steps": int(steps.group(1)),
        "maximum_force_kJ_mol_nm": float(force.group(1)),
        "nvt": _summarize_rows(_parse_xvg(observables)),
        "generated_force_field_coordinates": report.get(
            "generated_force_field_coordinates"
        ),
        "status": "completed",
    }


def main() -> int:
    """Run the allomorph validation matrix."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tleap-executable", default="tleap")
    parser.add_argument("--parmed-python", type=Path)
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--nvt-steps", type=int, default=1_000)
    parser.add_argument(
        "--work-directory",
        type=Path,
        default=Path("validation/work/glycam_allomorph_matrix"),
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("validation/results/glycam_allomorph_matrix/summary.json"),
    )
    args = parser.parse_args()
    if args.nvt_steps <= 0 or args.threads <= 0:
        parser.error("step and thread counts must be positive")
    args.work_directory.mkdir(parents=True, exist_ok=True)
    results = [
        _validate_one(args, allomorph, cells)
        for allomorph, cells in SYSTEMS.items()
    ]
    summary = {
        "validation": "periodic GLYCAM bulk allomorph matrix",
        "nvt_steps": args.nvt_steps,
        "time_step_ps": 0.001,
        "results": results,
        "interpretation_limit": (
            "Minimization and short NVT runs validate topology and numerical "
            "stability, not equilibrium material properties."
        ),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(args.summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
