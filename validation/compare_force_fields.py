"""Run matched GLYCAM and experimental OPLS-CM5 bulk validation workflows.

This script builds the same periodic bulk crystal with both force fields, runs
minimization, NVT, and NPT, and writes compact JSON evidence.
Large trajectories remain in the selected work directory and should not be
committed to git.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import statistics
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "validation"
ENERGY_TERMS = (
    "Potential",
    "Temperature",
    "Pressure",
    "Box-X",
    "Box-Y",
    "Box-Z",
    "Volume",
    "Density",
)


def _run(
    command: Sequence[str | Path],
    *,
    cwd: Path,
    log_path: Path,
    stdin: str | None = None,
) -> str:
    rendered = [str(value) for value in command]
    result = subprocess.run(
        rendered,
        cwd=cwd,
        input=stdin,
        text=True,
        capture_output=True,
        check=False,
    )
    output = result.stdout + result.stderr
    log_path.write_text("$ " + " ".join(rendered) + "\n\n" + output, encoding="utf-8")
    if result.returncode:
        msg = f"command failed ({result.returncode}); see {log_path}"
        raise RuntimeError(msg)
    return output


def _write_protocol(source: Path, destination: Path, *, nsteps: int) -> None:
    text = source.read_text(encoding="utf-8")
    rewritten, replacements = re.subn(
        r"(?m)^nsteps\s*=\s*\d+\s*$", f"nsteps                  = {nsteps}", text
    )
    if replacements != 1:
        msg = f"expected exactly one nsteps setting in {source}"
        raise ValueError(msg)
    destination.write_text(rewritten, encoding="utf-8")


def _gromacs_version(output: str) -> str:
    for line in output.splitlines():
        if line.strip().startswith("GROMACS version:"):
            return line.split(":", maxsplit=1)[1].strip()
    msg = "could not identify the GROMACS version"
    raise ValueError(msg)


def _parse_xvg(
    path: Path,
    terms: Sequence[str] = ENERGY_TERMS,
) -> tuple[tuple[float, ...], ...]:
    rows = tuple(
        tuple(float(value) for value in line.split())
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(("#", "@"))
    )
    expected_columns = len(terms) + 1
    if not rows or any(len(row) != expected_columns for row in rows):
        msg = f"unexpected GROMACS energy table in {path}"
        raise ValueError(msg)
    return rows


def _summarize_rows(
    rows: tuple[tuple[float, ...], ...],
    terms: Sequence[str] = ENERGY_TERMS,
) -> dict[str, Any]:
    summary: dict[str, Any] = {
        "start_ps": rows[0][0],
        "end_ps": rows[-1][0],
        "frames": len(rows),
    }
    for column, name in enumerate(terms, start=1):
        values = [row[column] for row in rows]
        summary[name] = {
            "mean": statistics.fmean(values),
            "standard_deviation": statistics.pstdev(values),
            "final_minus_initial": values[-1] - values[0],
        }
    return summary


def _extract_tpr_metadata(
    gmx: str, tpr: Path, directory: Path
) -> tuple[float, float]:
    dump = _run(
        [gmx, "dump", "-s", tpr.name],
        cwd=directory,
        log_path=directory / "dump.log",
    )
    charges = [
        float(match.group(1))
        for line in dump.splitlines()
        if "atom[" in line
        for match in [re.search(r"\bq=\s*([-+0-9.eE]+)", line)]
        if match is not None
    ]
    if not charges:
        msg = f"no atomic charges found in {tpr}"
        raise ValueError(msg)
    fudge_qq = re.search(r"\bfudgeQQ\s*=\s*([-+0-9.eE]+)", dump)
    if fudge_qq is None:
        msg = f"no electrostatic 1-4 scaling found in {tpr}"
        raise ValueError(msg)
    return math.fsum(charges), float(fudge_qq.group(1))


def _analyze_energy(gmx: str, directory: Path) -> dict[str, Any]:
    energy_path = directory / "npt.edr"
    table_path = directory / "npt_observables.xvg"
    _run(
        [gmx, "energy", "-f", energy_path.name, "-o", table_path.name, "-xvg", "none"],
        cwd=directory,
        log_path=directory / "energy.log",
        stdin="\n".join((*ENERGY_TERMS, "0", "")),
    )
    rows = _parse_xvg(table_path)
    midpoint = (rows[0][0] + rows[-1][0]) / 2
    late_rows = tuple(row for row in rows if row[0] >= midpoint)
    return {
        "complete_window": _summarize_rows(rows),
        "late_half_window": _summarize_rows(late_rows),
    }


def _build_command(args: argparse.Namespace, force_field: str, output: Path) -> list[str]:
    command = [
        sys.executable,
        "-m",
        "cellulose_forge.cli",
        "build",
        "--allomorph",
        args.allomorph,
        "--shape",
        "bulk",
        "--cells",
        *(str(value) for value in args.cells),
        "--chain-boundary",
        "periodic",
        "--force-field",
        force_field,
        "--output",
        str(output),
    ]
    if force_field == "glycam":
        command.extend(["--tleap-executable", args.tleap_executable])
        if args.parmed_python:
            command.extend(["--parmed-python", str(args.parmed_python)])
    else:
        command.extend(
            [
                "--opls-template",
                str(args.opls_template),
                "--opls-reference-gro",
                str(args.opls_reference_gro),
            ]
        )
    return command


def _run_force_field(
    args: argparse.Namespace,
    force_field: str,
    work_root: Path,
) -> dict[str, Any]:
    directory = work_root / force_field
    directory.mkdir(parents=True, exist_ok=True)
    _run(
        _build_command(args, force_field, directory),
        cwd=ROOT,
        log_path=directory / "build.log",
    )
    _write_protocol(
        VALIDATION / "periodic_pme_nvt.mdp",
        directory / "nvt.mdp",
        nsteps=args.nvt_steps,
    )
    _write_protocol(
        VALIDATION / "periodic_pme_npt.mdp",
        directory / "npt.mdp",
        nsteps=args.npt_steps,
    )
    commands = (
        (
            "grompp_em",
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
        ),
        ("mdrun_em", [args.gmx, "mdrun", "-deffnm", "em", "-nt", str(args.threads)]),
        (
            "grompp_nvt",
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
        ),
        ("mdrun_nvt", [args.gmx, "mdrun", "-deffnm", "nvt", "-nt", str(args.threads)]),
        (
            "grompp_npt",
            [
                args.gmx,
                "grompp",
                "-f",
                "npt.mdp",
                "-c",
                "nvt.gro",
                "-t",
                "nvt.cpt",
                "-p",
                "system.top",
                "-o",
                "npt.tpr",
                "-po",
                "npt_out.mdp",
            ],
        ),
        ("mdrun_npt", [args.gmx, "mdrun", "-deffnm", "npt", "-nt", str(args.threads)]),
    )
    command_outputs = {
        name: _run(command, cwd=directory, log_path=directory / f"{name}.log")
        for name, command in commands
    }
    report = json.loads((directory / "build_report.json").read_text(encoding="utf-8"))
    fudge_match = re.search(
        r"Generating 1-4 interactions: fudge = ([0-9.]+)",
        command_outputs["grompp_em"],
    )
    if fudge_match is None:
        msg = f"could not identify 1-4 scaling for {force_field}"
        raise ValueError(msg)
    system_charge, fudge_qq = _extract_tpr_metadata(
        args.gmx, directory / "em.tpr", directory
    )
    return {
        "force_field": force_field,
        "build": report,
        "system_charge_e": system_charge,
        "charge_per_chain_e": system_charge / report["chain_count"],
        "fudge_lj_14": float(fudge_match.group(1)),
        "fudge_qq_14": fudge_qq,
        "nvt_steps": args.nvt_steps,
        "npt_steps": args.npt_steps,
        "observables": _analyze_energy(args.gmx, directory),
        "status": "completed",
    }


def main() -> int:
    """Run the matched force-field comparison."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--opls-template", type=Path, required=True)
    parser.add_argument("--opls-reference-gro", type=Path, required=True)
    parser.add_argument(
        "--allomorph",
        choices=("I-alpha", "I-beta", "II", "III-I"),
        default="I-beta",
    )
    parser.add_argument("--cells", type=int, nargs=3, default=(4, 4, 4))
    parser.add_argument("--tleap-executable", default="tleap")
    parser.add_argument("--parmed-python", type=Path)
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--nvt-steps", type=int, default=10_000)
    parser.add_argument("--npt-steps", type=int, default=10_000)
    parser.add_argument(
        "--work-directory",
        type=Path,
        default=Path("validation/work/force_field_comparison"),
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("validation/results/force_field_comparison.json"),
    )
    args = parser.parse_args()
    if (
        args.nvt_steps <= 0
        or args.npt_steps <= 0
        or args.threads <= 0
        or any(value <= 0 for value in args.cells)
    ):
        parser.error("step counts and thread count must be positive")
    args.work_directory.mkdir(parents=True, exist_ok=True)
    version = _gromacs_version(
        _run(
            [args.gmx, "--version"],
            cwd=ROOT,
            log_path=args.work_directory / "gromacs_version.log",
        )
    )
    results = [
        _run_force_field(args, force_field, args.work_directory)
        for force_field in ("glycam", "opls-cm5")
    ]
    summary = {
        "comparison": f"matched periodic cellulose {args.allomorph} bulk crystal",
        "cells": list(args.cells),
        "gromacs": version,
        "protocol": {
            "minimization": "periodic_pme_minim.mdp",
            "nvt_steps": args.nvt_steps,
            "npt_steps": args.npt_steps,
            "time_step_ps": 0.001,
            "pressure_coupling": "isotropic C-rescale",
        },
        "results": results,
        "interpretation_limit": (
            "Engineering validation only; short trajectories do not establish "
            "force-field accuracy or equilibrium material properties."
        ),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(args.summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
