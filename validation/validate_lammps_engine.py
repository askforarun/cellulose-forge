"""Run the dependency-gated LAMMPS cross-engine validation suite.

The tests regenerate finite and periodic I-beta and periodic cellulose II
systems, evaluate them with LAMMPS, and compare energies with independent
OpenMM or GROMACS calculations.
The retained reference inputs and outputs live under
``validation/results/lammps_cross_engine_files``.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    """Run the cross-engine tests using explicitly selectable executables."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--openmm-python",
        help="Python interpreter containing OpenMM; defaults to the test interpreter",
    )
    parser.add_argument(
        "--lammps-executable",
        default="lmp",
        help="LAMMPS executable name or path (default: lmp)",
    )
    parser.add_argument(
        "--bundle-directory",
        type=Path,
        default=ROOT / "validation/results/lammps_cross_engine_files",
        help="destination for the trimmed, path-redacted evidence bundle",
    )
    arguments = parser.parse_args()
    environment = os.environ.copy()
    environment["CELLULOSE_FORGE_LAMMPS_EXECUTABLE"] = arguments.lammps_executable
    if arguments.openmm_python:
        environment["CELLULOSE_FORGE_OPENMM_PYTHON"] = arguments.openmm_python
    with tempfile.TemporaryDirectory(prefix="cellulose_forge_lammps_bundle_") as temporary:
        stage = Path(temporary)
        environment["CELLULOSE_FORGE_LAMMPS_BUNDLE_DIR"] = str(stage)
        completed = subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                "tests/integration/test_lammps_workflows.py",
                "-m",
                "lammps and openmm and gromacs and scientific_integration",
                "-q",
            ],
            cwd=ROOT,
            env=environment,
            check=False,
        )
        if completed.returncode != 0:
            return completed.returncode
        _publish_bundle(stage, arguments.bundle_directory)
    print(f"Reference files: {arguments.bundle_directory}")
    return 0


def _publish_bundle(stage: Path, destination: Path) -> None:
    """Replace generated cases and write a path-redacted comparison summary."""
    destination.mkdir(parents=True, exist_ok=True)
    comparisons: dict[str, dict[str, object]] = {}
    generated_on = date.today().isoformat()
    for source_case in sorted(path for path in stage.iterdir() if path.is_dir()):
        comparison_path = source_case / ".comparison.json"
        comparison = json.loads(comparison_path.read_text(encoding="utf-8"))
        generated_on = str(comparison.pop("generated_on", generated_on))
        comparisons[source_case.name] = comparison
        comparison_path.unlink()
        target_case = destination / source_case.name
        if target_case.exists():
            shutil.rmtree(target_case)
        shutil.copytree(source_case, target_case)

    lammps_log = next(destination.glob("*/lammps/log.lammps"))
    lammps_match = re.search(r"LAMMPS \(([^)]+)\)", lammps_log.read_text(encoding="utf-8"))
    gromacs_log = next(destination.glob("*/gromacs/*single_point.log"))
    gromacs_match = re.search(
        r"GROMACS version:\s+([^\s]+)", gromacs_log.read_text(encoding="utf-8")
    )
    openmm_report = next(destination.glob("glycam_*/common/build_report.json"))
    report = json.loads(openmm_report.read_text(encoding="utf-8"))
    summary = {
        "generated_on": generated_on,
        "engines": {
            "gromacs": gromacs_match.group(1) if gromacs_match else "unknown",
            "lammps": lammps_match.group(1) if lammps_match else "unknown",
            "openmm": report["openmm_validation"]["openmm_version"],
        },
        "cases": comparisons,
    }
    (destination / "comparison.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    forbidden = re.compile(r"/Users|/private|/home/|miniconda|pytest-")
    for path in destination.rglob("*"):
        if path.is_file() and forbidden.search(path.read_text(encoding="utf-8", errors="replace")):
            raise RuntimeError(f"local path was not redacted from {path}")


if __name__ == "__main__":
    raise SystemExit(main())
