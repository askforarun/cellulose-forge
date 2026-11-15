"""Build and run the three-case experimental NAMD acceptance matrix.

This script requires real ``namd3`` and ``gmx`` executables. It stops on the
first failed build, topology comparison, or single-point energy comparison;
it never synthesizes engine output.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

CASES = (
    ("cellulose_ibeta", "cellulose", "I-beta", (4, 4, 3)),
    ("cellulose_ialpha", "cellulose", "I-alpha", (4, 5, 6)),
    ("chitin_alpha", "chitin", "alpha", (7, 2, 3)),
)


def _executable(value: str) -> str:
    resolved = shutil.which(value)
    if resolved is None:
        raise SystemExit(f"required executable is unavailable: {value}")
    return resolved


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameter-file", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("namd_acceptance_matrix"))
    parser.add_argument("--cellulose-forge", default="cellulose-forge")
    parser.add_argument("--namd", default="namd3")
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--namd-threads", type=int, default=8)
    arguments = parser.parse_args()

    forge = _executable(arguments.cellulose_forge)
    namd = _executable(arguments.namd)
    gmx = _executable(arguments.gmx)
    parameter_file = arguments.parameter_file.resolve()
    if not parameter_file.is_file():
        raise SystemExit(f"CHARMM parameter file does not exist: {parameter_file}")

    root = Path(__file__).resolve().parents[2]
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=root, check=True, capture_output=True, text=True
    ).stdout.strip()
    tracked_status = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=no"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    case_results: dict[str, object] = {}
    results: dict[str, object] = {
        "schema_version": 1,
        "cellulose_forge_commit": commit,
        "tracked_tree_clean": not tracked_status,
        "executables": {"cellulose_forge": forge, "namd": namd, "gromacs": gmx},
        "cases": case_results,
    }
    for name, polymer, allomorph, cells in CASES:
        directory = output / name
        command = [
            forge,
            "build",
            "--polymer",
            polymer,
            "--allomorph",
            allomorph,
            "--shape",
            "bulk",
            "--cells",
            *(str(value) for value in cells),
            "--chain-boundary",
            "periodic",
            "--force-field",
            "charmm36",
            "--md-engine",
            "namd",
            "--md-engine",
            "gromacs",
            "--charmm-parameter-file",
            str(parameter_file),
            "--output",
            str(directory),
        ]
        subprocess.run(command, cwd=root, check=True)
        subprocess.run(
            [
                sys.executable,
                str(root / "validation" / "namd" / "compare_bundle.py"),
                str(directory),
                "--output",
                str(directory / "bundle_comparison.json"),
            ],
            cwd=root,
            check=True,
        )
        environment = os.environ.copy()
        environment.update(
            {
                "NAMD_EXE": namd,
                "GMX_EXE": gmx,
                "PYTHON_EXE": sys.executable,
                "NAMD_THREADS": str(arguments.namd_threads),
            }
        )
        subprocess.run(
            ["bash", "run_single_point_comparison.sh"],
            cwd=directory,
            env=environment,
            check=True,
        )
        comparison = json.loads(
            (directory / "single_point_energy_comparison.json").read_text(encoding="utf-8")
        )
        if comparison.get("passed") is not True:
            raise SystemExit(f"single-point comparison did not pass: {name}")
        case_results[name] = {
            "polymer": polymer,
            "allomorph": allomorph,
            "cells": list(cells),
            "comparison": comparison,
        }

    (output / "acceptance_matrix.json").write_text(
        json.dumps(results, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
