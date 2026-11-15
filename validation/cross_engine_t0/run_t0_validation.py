#!/usr/bin/env python3
"""Run and compare normalized t=0 energy evaluations across MD engines."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

CANONICAL_TERMS = {
    "bond",
    "angle",
    "urey_bradley",
    "dihedral",
    "proper_dihedral",
    "improper",
    "improper_dihedral",
    "lj_14",
    "coulomb_14",
    "lj_nonbonded",
    "coulomb_nonbonded",
    "vdw",
    "electrostatics",
    "potential",
}


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("manifest.json"))
    parser.add_argument("--output", type=Path, default=Path("t0_energy_comparison.json"))
    parser.add_argument("--no-run", action="store_true", help="compare existing energy.json files")
    parser.add_argument(
        "--engine",
        action="append",
        dest="selected_engines",
        help="run only this engine subfolder; repeat to select more than one",
    )
    return parser.parse_args()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_energy(
    path: Path,
    engine: str,
    coordinate_sha256: str,
    required_terms: set[str],
) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema_version") != 1:
        raise ValueError(f"{engine}: unsupported energy schema")
    if data.get("engine") != engine:
        raise ValueError(f"{engine}: energy.json names engine {data.get('engine')!r}")
    if data.get("units") != "kJ/mol":
        raise ValueError(f"{engine}: energies must be normalized to kJ/mol")
    if data.get("coordinate_sha256") != coordinate_sha256:
        raise ValueError(f"{engine}: coordinate checksum differs from the shared coordinates")
    terms = data.get("terms")
    if not isinstance(terms, dict) or not terms:
        raise ValueError(f"{engine}: energy.json contains no terms")
    unknown = set(terms) - CANONICAL_TERMS
    if unknown:
        raise ValueError(f"{engine}: unknown canonical terms: {sorted(unknown)}")
    for name, value in terms.items():
        if not isinstance(value, (int, float)):
            raise TypeError(f"{engine}: term {name!r} is not numeric")
    if "potential" not in terms:
        raise ValueError(f"{engine}: potential energy is required")
    missing = required_terms - set(terms)
    if missing:
        raise ValueError(f"{engine}: required comparable terms are missing: {sorted(missing)}")
    if not isinstance(data.get("corrections", []), list):
        raise TypeError(f"{engine}: corrections must be a list")
    return data


def _comparison(
    reference: dict[str, Any],
    candidate: dict[str, Any],
    absolute: float,
    relative: float,
) -> dict[str, Any]:
    common = sorted(set(reference["terms"]) & set(candidate["terms"]))
    if not common:
        raise ValueError("engines report no common energy terms")
    terms: dict[str, Any] = {}
    passed = True
    for name in common:
        expected = float(reference["terms"][name])
        observed = float(candidate["terms"][name])
        difference = observed - expected
        scale = max(abs(expected), abs(observed), 1.0)
        allowed = max(absolute, relative * scale)
        term_passed = abs(difference) <= allowed
        passed &= term_passed
        terms[name] = {
            "reference_kj_mol": expected,
            "candidate_kj_mol": observed,
            "difference_kj_mol": difference,
            "relative_difference": abs(difference) / scale,
            "allowed_difference_kj_mol": allowed,
            "passed": term_passed,
        }
    return {"passed": passed, "common_terms": terms}


def _selected_engines(arguments: argparse.Namespace, manifest: dict[str, Any]) -> set[str]:
    selected = set(arguments.selected_engines or manifest["engines"])
    unknown = selected - set(manifest["engines"])
    if unknown:
        raise ValueError(f"selected engines are absent from the manifest: {sorted(unknown)}")
    reference = manifest["reference_engine"]
    if arguments.selected_engines and reference not in selected:
        raise ValueError(f"selected engines must include reference engine {reference!r}")
    return selected


def _collect_energies(
    arguments: argparse.Namespace,
    root: Path,
    manifest: dict[str, Any],
    checksum: str,
    required_terms: set[str],
) -> tuple[dict[str, dict[str, Any]], dict[str, str]]:
    selected = _selected_engines(arguments, manifest)
    energies: dict[str, dict[str, Any]] = {}
    skipped: dict[str, str] = {}
    for engine, specification in manifest["engines"].items():
        if engine not in selected:
            continue
        directory = root / specification["directory"]
        if not directory.is_dir():
            if specification.get("optional", False):
                skipped[engine] = "engine directory is absent"
                continue
            raise FileNotFoundError(f"{engine}: missing directory {directory}")
        if not arguments.no_run:
            try:
                subprocess.run(specification["command"], cwd=directory, check=True)
            except (FileNotFoundError, subprocess.CalledProcessError) as error:
                if specification.get("optional", False):
                    skipped[engine] = str(error)
                    continue
                raise
        energies[engine] = _load_energy(directory / "energy.json", engine, checksum, required_terms)
    return energies, skipped


def main() -> int:
    args = _arguments()
    manifest_path = args.manifest.resolve()
    root = manifest_path.parent
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    coordinate = root / manifest["coordinate_file"]
    actual_checksum = _sha256(coordinate)
    expected_checksum = manifest["coordinate_sha256"]
    if actual_checksum != expected_checksum:
        raise ValueError("shared coordinate checksum does not match manifest")

    required_terms = set(manifest["required_terms"])
    unknown_required = required_terms - CANONICAL_TERMS
    if unknown_required:
        raise ValueError(f"manifest contains unknown required terms: {sorted(unknown_required)}")
    if "potential" not in required_terms:
        raise ValueError("manifest required_terms must include potential")

    reference_name = manifest["reference_engine"]
    energies, skipped = _collect_energies(args, root, manifest, expected_checksum, required_terms)

    if reference_name not in energies:
        raise ValueError(f"reference engine {reference_name!r} did not produce energies")
    absolute = float(manifest["absolute_tolerance_kj_mol"])
    relative = float(manifest["relative_tolerance"])
    comparisons = {
        engine: _comparison(energies[reference_name], values, absolute, relative)
        for engine, values in energies.items()
        if engine != reference_name
    }
    passed = bool(comparisons) and all(item["passed"] for item in comparisons.values())
    report = {
        "schema_version": 1,
        "system": manifest["system"],
        "coordinate_file": manifest["coordinate_file"],
        "coordinate_sha256": expected_checksum,
        "reference_engine": reference_name,
        "tolerances": {"absolute_kj_mol": absolute, "relative": relative},
        "required_terms": sorted(required_terms),
        "energies": energies,
        "comparisons": comparisons,
        "skipped_optional_engines": skipped,
        "passed": passed,
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"t=0 cross-engine validation: {'PASS' if passed else 'FAIL'}")
    print(f"report: {args.output}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
