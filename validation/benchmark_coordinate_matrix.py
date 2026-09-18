"""Benchmark deterministic batch construction of cellulose bulk crystals.

Run from the repository root with::

    uv run python validation/benchmark_coordinate_matrix.py

The benchmark invokes the public command-line interface for every build. It
uses temporary coordinate files and retains only a machine-readable summary.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

from cellulose_forge import __version__

ALLOMORPHS = ("I-alpha", "I-beta", "II", "III-I")
CELL_MATRICES = ((1, 1, 1), (2, 2, 2), (4, 4, 5))
MINIMUM_DETERMINISM_REPEATS = 2


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atom_count(path: Path) -> int:
    with path.open(encoding="utf-8") as stream:
        return sum(line.startswith(("ATOM ", "HETATM ")) for line in stream)


def _run_build(allomorph: str, cells: tuple[int, int, int], output: Path) -> float:
    command = [
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
        "--force-field",
        "none",
        "--output",
        str(output),
    ]
    started = time.perf_counter()
    subprocess.run(command, check=True, capture_output=True, text=True)
    return time.perf_counter() - started


def run_benchmark(repeats: int) -> dict[str, Any]:
    """Run the coordinate matrix and return its summary."""
    if repeats < MINIMUM_DETERMINISM_REPEATS:
        msg = "repeats must be at least 2 so determinism can be checked"
        raise ValueError(msg)

    results: list[dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="cellulose_forge_matrix_") as temporary:
        root = Path(temporary)
        for allomorph in ALLOMORPHS:
            for cells in CELL_MATRICES:
                timings: list[float] = []
                hashes: list[str] = []
                atom_counts: list[int] = []
                for repeat in range(repeats):
                    output = root / f"{allomorph}_{cells}_{repeat}.cif"
                    timings.append(_run_build(allomorph, cells, output))
                    hashes.append(_sha256(output))
                    atom_counts.append(_atom_count(output))
                results.append(
                    {
                        "allomorph": allomorph,
                        "cells": list(cells),
                        "atom_count": atom_counts[0],
                        "repeats": repeats,
                        "deterministic": len(set(hashes)) == 1,
                        "sha256": hashes[0],
                        "timings_seconds": timings,
                        "median_seconds": statistics.median(timings),
                    }
                )

    return {
        "benchmark": "coordinate-only bulk construction matrix",
        "command": "uv run python validation/benchmark_coordinate_matrix.py",
        "cellulose_forge_version": __version__,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "system_count": len(results),
        "build_count": len(results) * repeats,
        "all_deterministic": all(result["deterministic"] for result in results),
        "total_atoms_per_matrix": sum(result["atom_count"] for result in results),
        "total_measured_seconds": sum(
            sum(result["timings_seconds"]) for result in results
        ),
        "results": results,
    }


def main() -> int:
    """Run the benchmark and write its JSON summary."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("validation/results/high_throughput_coordinate_matrix.json"),
    )
    args = parser.parse_args()
    summary = run_benchmark(args.repeats)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(
        f"Built {summary['build_count']} files for {summary['system_count']} systems; "
        f"deterministic={summary['all_deterministic']}; "
        f"elapsed={summary['total_measured_seconds']:.3f} s"
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
