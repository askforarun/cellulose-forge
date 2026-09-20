"""Compare periodic GLYCAM and experimental OPLS-CM5 bulk allomorphs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from compare_force_fields import ROOT, _gromacs_version, _run, _run_force_field

SYSTEMS = {
    "I-alpha": (4, 4, 4),
    "I-beta": (4, 4, 4),
    "II": (4, 4, 4),
    "III-I": (6, 4, 4),
}


def _safe_name(allomorph: str) -> str:
    return allomorph.lower().replace("-", "_")


def main() -> int:
    """Run matched NVT/NPT comparisons for all four allomorphs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--opls-template", type=Path, required=True)
    parser.add_argument("--opls-reference-gro", type=Path, required=True)
    parser.add_argument("--tleap-executable", default="tleap")
    parser.add_argument("--parmed-python", type=Path)
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--nvt-steps", type=int, default=10_000)
    parser.add_argument("--npt-steps", type=int, default=10_000)
    parser.add_argument(
        "--work-directory",
        type=Path,
        default=Path("validation/work/allomorph_force_field_comparison"),
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path(
            "validation/results/allomorph_force_field_comparison/summary.json"
        ),
    )
    args = parser.parse_args()
    if args.nvt_steps <= 0 or args.npt_steps <= 0 or args.threads <= 0:
        parser.error("step counts and thread count must be positive")
    args.work_directory.mkdir(parents=True, exist_ok=True)
    version = _gromacs_version(
        _run(
            [args.gmx, "--version"],
            cwd=ROOT,
            log_path=args.work_directory / "gromacs_version.log",
        )
    )
    systems = []
    for allomorph, cells in SYSTEMS.items():
        args.allomorph = allomorph
        args.cells = cells
        directory = args.work_directory / _safe_name(allomorph)
        systems.append(
            {
                "allomorph": allomorph,
                "cells": list(cells),
                "force_fields": [
                    _run_force_field(args, force_field, directory)
                    for force_field in ("glycam", "opls-cm5")
                ],
            }
        )
    summary = {
        "comparison": "matched periodic cellulose bulk allomorph matrix",
        "gromacs": version,
        "protocol": {
            "minimization": "periodic_pme_minim.mdp",
            "nvt_steps": args.nvt_steps,
            "npt_steps": args.npt_steps,
            "time_step_ps": 0.001,
            "pressure_coupling": "isotropic C-rescale",
        },
        "systems": systems,
        "interpretation_limit": (
            "Engineering comparison only; short isotropic NPT trajectories do not "
            "establish force-field accuracy or equilibrium lattice properties."
        ),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(args.summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
