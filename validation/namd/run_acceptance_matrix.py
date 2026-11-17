"""Build and run the three-case experimental NAMD acceptance matrix.

This script requires real ``namd3`` and ``gmx`` executables. It stops on the
first failed build, topology comparison, single-point energy comparison, or
production-input dynamics smoke test; it never synthesizes engine output.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from cellulose_forge.energy_drivers import (  # noqa: E402
    compare_reports,
    run_gromacs_single_point,
    run_namd_single_point,
    validate_potential_sum,
)
from cellulose_forge.energy_drivers.coordinates import (  # noqa: E402
    read_gro_frame,
    write_coordinate_reference,
)

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


def main() -> None:  # noqa: PLR0915 - linear fail-fast acceptance workflow
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameter-file", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("namd_acceptance_matrix"))
    parser.add_argument(
        "--summary-output",
        type=Path,
        help="optional second location for the compact acceptance_matrix.json",
    )
    parser.add_argument("--cellulose-forge", default="cellulose-forge")
    parser.add_argument("--namd", default="namd3")
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument(
        "--audit-python",
        default=sys.executable,
        help="Python interpreter containing ParmEd for compare_bundle.py",
    )
    parser.add_argument("--namd-threads", type=int, default=8)
    parser.add_argument(
        "--case",
        action="append",
        dest="selected_cases",
        choices=[name for name, *_ in CASES],
        help="run only this acceptance case; repeat to select more than one",
    )
    arguments = parser.parse_args()

    forge = _executable(arguments.cellulose_forge)
    namd = _executable(arguments.namd)
    gmx = _executable(arguments.gmx)
    audit_python = _executable(arguments.audit_python)
    subprocess.run(
        [audit_python, "-c", "import parmed"], check=True, capture_output=True, text=True
    )
    parameter_file = arguments.parameter_file.resolve()
    if not parameter_file.is_file():
        raise SystemExit(f"CHARMM parameter file does not exist: {parameter_file}")

    root = ROOT
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
    forge_version = subprocess.run(
        [forge, "--version"], check=True, capture_output=True, text=True
    ).stdout.strip()
    gromacs_version_output = subprocess.run(
        [gmx, "--version"], check=True, capture_output=True, text=True
    ).stdout
    gromacs_version = next(
        line.split(":", 1)[1].strip()
        for line in gromacs_version_output.splitlines()
        if line.startswith("GROMACS version:")
    )
    case_results: dict[str, object] = {}
    software_versions = {
        "cellulose_forge": forge_version,
        "gromacs": gromacs_version,
        "parmed": subprocess.run(
            [audit_python, "-c", "import parmed; print(parmed.__version__)"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip(),
    }
    results: dict[str, object] = {
        "schema_version": 1,
        "provenance": {
            "cellulose_forge_commit": commit,
            "tracked_tree_clean": not tracked_status,
            "software_versions": software_versions,
        },
        "parameter_sha256": hashlib.sha256(parameter_file.read_bytes()).hexdigest(),
        "single_point_protocol": {
            "constraints": "none",
            "lj_cutoff_angstrom": 12.0,
            "lj_switching": "off in both engines",
            "pme_tolerance": 1.0e-6,
        },
        "cases": case_results,
    }
    selected = set(arguments.selected_cases or (name for name, *_ in CASES))
    for name, polymer, allomorph, cells in CASES:
        if name not in selected:
            continue
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
                audit_python,
                str(root / "validation" / "namd" / "compare_bundle.py"),
                str(directory),
                "--output",
                str(directory / "bundle_comparison.json"),
            ],
            cwd=root,
            check=True,
        )
        positions, box = read_gro_frame(directory / "single_point_coordinates.gro")
        coordinate_reference = directory / "coordinate_reference.json"
        coordinate_sha256 = write_coordinate_reference(coordinate_reference, positions, box)
        gromacs_report = run_gromacs_single_point(
            directory,
            coordinate_file=directory / "single_point_coordinates.gro",
            coordinate_reference_file=coordinate_reference,
            topology_file=directory / "system.top",
            mdp_file=directory / "single_point_gromacs.mdp",
            coordinate_sha256=coordinate_sha256,
            executable=gmx,
        )
        namd_report = run_namd_single_point(
            directory,
            configuration_file=directory / "single_point_namd.conf",
            coordinate_reference_file=coordinate_reference,
            binary_coordinate_file=directory / "single_point.coor",
            coordinate_sha256=coordinate_sha256,
            executable=namd,
            threads=arguments.namd_threads,
        )
        for energy_report in (gromacs_report, namd_report):
            validate_potential_sum(energy_report, absolute_kj_mol=0.1, relative=1.0e-6)
        raw_gromacs = gromacs_report.extra["raw_terms"]
        electrostatic_scale = sum(
            abs(float(raw_gromacs.get(term, 0.0)))
            for term in ("Coulomb-14", "Coulomb (SR)", "Coul. recip.")
        )
        comparison_rows = compare_reports(
            gromacs_report,
            namd_report,
            absolute_kj_mol=0.1,
            relative=5.0e-5,
            term_tolerances={"electrostatics": (max(0.1, 5.0e-5 * electrostatic_scale), 5.0e-5)},
        )
        comparison = {
            "schema_version": 1,
            "coordinate_sha256": coordinate_sha256,
            "reference_engine": "gromacs",
            "energy_reports": {
                "gromacs": gromacs_report.to_dict(),
                "namd": namd_report.to_dict(),
            },
            "comparison": comparison_rows,
            "passed": comparison_rows["passed"],
        }
        (directory / "single_point_energy_comparison.json").write_text(
            json.dumps(comparison, indent=2) + "\n", encoding="utf-8"
        )
        with (directory / "production_dynamics.log").open("w", encoding="utf-8") as log:
            subprocess.run(
                [namd, f"+p{arguments.namd_threads}", "run_namd.conf"],
                cwd=directory,
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
        subprocess.run(
            [
                sys.executable,
                str(root / "validation" / "namd" / "analyze_dynamics.py"),
                str(directory),
                "--output",
                str(directory / "dynamics_validation.json"),
            ],
            cwd=root,
            check=True,
        )
        bundle_comparison = json.loads(
            (directory / "bundle_comparison.json").read_text(encoding="utf-8")
        )
        dynamics_validation = json.loads(
            (directory / "dynamics_validation.json").read_text(encoding="utf-8")
        )
        if comparison.get("passed") is not True:
            raise SystemExit(f"single-point comparison did not pass: {name}")
        if bundle_comparison["atom_row_mismatches"] or not all(
            bundle_comparison["term_counters_equal"].values()
        ):
            raise SystemExit(f"topology comparison did not pass: {name}")
        if dynamics_validation.get("passed") is not True:
            raise SystemExit(f"production dynamics validation did not pass: {name}")
        namd_log = (directory / "single_point_namd.log").read_text(
            encoding="utf-8", errors="replace"
        )
        namd_match = re.search(r"^Info: NAMD\s+([^\s]+)", namd_log, flags=re.MULTILINE)
        if namd_match is None:
            raise SystemExit(f"cannot identify NAMD version: {name}")
        software_versions["namd"] = namd_match.group(1)
        case_results[name] = {
            "polymer": polymer,
            "allomorph": allomorph,
            "cells": list(cells),
            "bundle_comparison": bundle_comparison,
            "comparison": comparison,
            "production_dynamics": dynamics_validation,
        }

    rendered = json.dumps(results, indent=2) + "\n"
    (output / "acceptance_matrix.json").write_text(rendered, encoding="utf-8")
    if arguments.summary_output is not None:
        summary_output = arguments.summary_output.resolve()
        summary_output.parent.mkdir(parents=True, exist_ok=True)
        summary_output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
