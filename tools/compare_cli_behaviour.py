"""Compare ``cellulose-forge build`` behaviour between two source trees.

Each combination of polymer, shape, force field, MD-engine set and chain
boundary is run once against each tree. A combination matches when the exit
codes are equal and, for refusals, the normalized final error line is equal or,
for successful builds, every output file is byte-identical after the output
directory path, the source-tree path and tool write-time stamps are replaced by
placeholders. Outputs whose files contain the
same lines in a different order are reported separately as
``same_up_to_line_order`` and do not fail the comparison. Run from the repository root::

    git worktree add --detach /tmp/cf-old <old-commit>
    python tools/compare_cli_behaviour.py --old-src /tmp/cf-old/src --new-src src \
        --output comparison.json

Sizes are chosen so that force-field builds pass the cutoff check (bulk boxes,
DP8 chains, 20-angstrom vacuum, 12-angstrom padding for finite chains). The
interpreter must provide the optional dependencies (OpenFF, OpenMM, AmberTools)
that the compared builds need; a missing dependency affects both trees equally
and appears as a matching refusal. The script exits 1 if any combination
differs, so intended changes are reviewed from the JSON report.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

_RUN = (
    "import sys, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0, sys.argv[1]); "
    "from cellulose_forge.cli import main; raise SystemExit(main(sys.argv[2:]))"
)
POLYMERS = {
    "cellulose": (["--allomorph", "I-beta"], ["--cells", "3", "3", "4"]),
    "chitin": (["--polymer", "chitin", "--allomorph", "beta"], ["--cells", "5", "3", "4"]),
    "mannan": (["--polymer", "mannan", "--allomorph", "I"], ["--cells", "3", "3", "3"]),
}
SHAPES = ("unit-cell", "bulk", "crystal", "slab", "fibril")
FORCE_FIELDS = ("none", "glycam", "charmm36", "openff", "opls-aa", "opls-aa-sei")
ENGINE_SETS = ((), ("gromacs",), ("openmm",), ("lammps",), ("amber",), ("gromacs", "amber"))
BOUNDARIES = ("finite", "periodic")
_SLAB = [
    "--face", "1", "0", "0", "--layers", "3", "--chains-per-layer", "4",
    "--degree-polymerization", "8", "--vacuum", "20",
]  # fmt: skip
_FIBRIL = ["--fibril-model", "36-chain", "--degree-polymerization", "8", "--vacuum", "20"]
_PADDING = ["--padding-longitudinal", "12"]
_PLACEHOLDER = "<OUTPUT>"
_SOURCE_PLACEHOLDER = "<SOURCE>"
# Write-time stamps that tleap (Amber prmtop) and ParmEd (GROMACS topology
# header) embed in GLYCAM output; they differ between otherwise identical runs.
_WRITE_TIME_STAMPS = (
    (re.compile(rb"(%VERSION\s+VERSION_STAMP = \S+\s+DATE = )[^\n]*"), rb"\1<DATE>"),
    (re.compile(rb"(;\s+At date:)[^\n]*"), rb"\1<DATE>"),
)


def combination_arguments(
    polymer: str, shape: str, force_field: str, engines: tuple[str, ...], boundary: str
) -> list[str]:
    """Return the build arguments (without --output) for one combination."""
    identity, cells = POLYMERS[polymer]
    size = {"unit-cell": [], "bulk": cells, "crystal": cells, "slab": _SLAB, "fibril": _FIBRIL}
    arguments = ["build", *identity, "--shape", shape, *size[shape]]
    arguments += ["--force-field", force_field, "--chain-boundary", boundary]
    if boundary == "finite" and force_field != "none" and shape in {"slab", "fibril"}:
        arguments += _PADDING
    for engine in engines:
        arguments += ["--md-engine", engine]
    return arguments


def normalize_message(stderr: str, output: str) -> str:
    """Return the final error line with temporary paths and whitespace normalized."""
    lines = [line for line in stderr.strip().splitlines() if line.strip()]
    message = lines[-1] if lines else ""
    message = message.replace(output, _PLACEHOLDER)
    message = re.sub(r"/(?:private/)?(?:var|tmp)/\S+", "<TMP>", message)
    return " ".join(message.split())


def output_fingerprint(output: Path, source: Path | None = None) -> dict[str, str] | None:
    """Hash every file under ``output`` with the output and source paths replaced.

    ParmEd records the converter command line, including the source-tree path,
    in GLYCAM GROMACS topology headers, so ``source`` is replaced as well.

    ``exact`` hashes file bytes. ``line_multiset`` hashes each file's sorted
    lines, so two outputs that differ only in line order (OpenFF Interchange
    does not write GROMACS ``[ pairs ]`` in a stable order) can be told apart
    from outputs whose content differs.
    """
    if not output.exists():
        return None
    files = [output] if output.is_file() else sorted(p for p in output.rglob("*") if p.is_file())
    exact = hashlib.sha256()
    multiset = hashlib.sha256()
    for path in files:
        relative = "." if path == output else str(path.relative_to(output))
        content = path.read_bytes().replace(str(output).encode(), _PLACEHOLDER.encode())
        if source is not None:
            content = content.replace(str(source).encode(), _SOURCE_PLACEHOLDER.encode())
        for stamp, replacement in _WRITE_TIME_STAMPS:
            content = stamp.sub(replacement, content)
        exact.update(f"{relative}\0{len(content)}\0".encode())
        exact.update(content)
        multiset.update(f"{relative}\0".encode())
        multiset.update(b"\n".join(sorted(content.splitlines())))
    return {"exact": exact.hexdigest(), "line_multiset": multiset.hexdigest()}


def run_once(
    python: str, source: Path, arguments: list[str], force_field: str
) -> dict[str, object]:
    """Run one build against ``source`` and summarize the outcome."""
    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory) / ("model.pdb" if force_field == "none" else "build")
        environment = dict(os.environ)
        environment["PATH"] = os.pathsep.join((str(Path(python).parent), environment["PATH"]))
        process = subprocess.run(
            [python, "-c", _RUN, str(source), *arguments, "--output", str(output)],
            capture_output=True,
            text=True,
            timeout=1800,
            check=False,
            env=environment,
        )
        return {
            "exit_code": process.returncode,
            "message": normalize_message(process.stderr, str(output)) if process.returncode else "",
            "fingerprint": output_fingerprint(output, source) if process.returncode == 0 else None,
        }


def compare(
    python: str, old: Path, new: Path, combination: tuple[str, str, str, tuple[str, ...], str]
) -> dict[str, object]:
    polymer, shape, force_field, engines, boundary = combination
    arguments = combination_arguments(polymer, shape, force_field, engines, boundary)
    before = run_once(python, old, arguments, force_field)
    after = run_once(python, new, arguments, force_field)
    if before["exit_code"] != after["exit_code"]:
        verdict = "outcome_changed"
    elif before["exit_code"] == 0:
        old_print, new_print = before["fingerprint"], after["fingerprint"]
        if old_print == new_print:
            verdict = "same"
        elif old_print and new_print and old_print["line_multiset"] == new_print["line_multiset"]:  # type: ignore[index]
            verdict = "same_up_to_line_order"
        else:
            verdict = "output_changed"
    else:
        verdict = "same" if before["message"] == after["message"] else "message_changed"
    return {
        "combination": {
            "polymer": polymer,
            "shape": shape,
            "force_field": force_field,
            "md_engines": list(engines),
            "chain_boundary": boundary,
        },
        "verdict": verdict,
        "old": before,
        "new": after,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--old-src", type=Path, required=True, help="src/ of the reference tree")
    parser.add_argument("--new-src", type=Path, required=True, help="src/ of the tree under test")
    parser.add_argument("--python", default=sys.executable, help="interpreter for both trees")
    parser.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    parser.add_argument("--output", type=Path, help="write the full comparison as JSON")
    arguments = parser.parse_args(argv)
    combinations = list(itertools.product(POLYMERS, SHAPES, FORCE_FIELDS, ENGINE_SETS, BOUNDARIES))
    old, new = arguments.old_src.resolve(), arguments.new_src.resolve()
    with ThreadPoolExecutor(max_workers=arguments.jobs) as pool:
        results = list(
            pool.map(lambda item: compare(arguments.python, old, new, item), combinations)
        )
    counts: dict[str, int] = {}
    built = 0
    for result in results:
        counts[str(result["verdict"])] = counts.get(str(result["verdict"]), 0) + 1
        built += result["old"]["exit_code"] == 0  # type: ignore[index]
    summary = {"combinations": len(results), "built_successfully_on_old": built, **counts}
    print(json.dumps(summary))
    for result in results:
        if result["verdict"] not in {"same", "same_up_to_line_order"}:
            combo = result["combination"]
            print(
                f"{result['verdict']}: {combo}\n  old: {result['old']['message'] or 'built'}"  # type: ignore[index]
                f"\n  new: {result['new']['message'] or 'built'}"  # type: ignore[index]
            )
    if arguments.output is not None:
        arguments.output.write_text(
            json.dumps({"summary": summary, "results": results}, indent=1) + "\n",
            encoding="utf-8",
        )
    matching = counts.get("same", 0) + counts.get("same_up_to_line_order", 0)
    return 0 if matching == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
