#!/usr/bin/env python3
"""Run and compare normalized single-point energies across MD engines."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

from cellulose_forge.energy_drivers.report import (
    CANONICAL_TERMS,
    EnergyReport,
    compare_reports,
    load_energy_report,
    validate_potential_sum,
)


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("manifest.json"))
    parser.add_argument("--output", type=Path, default=Path("single_point_energy_comparison.json"))
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
) -> EnergyReport:
    return load_energy_report(
        path,
        expected_engine=engine,
        coordinate_sha256=coordinate_sha256,
        required_terms=required_terms,
    )


def _selected_engines(arguments: argparse.Namespace, manifest: dict[str, Any]) -> set[str]:
    selected = set(arguments.selected_engines or manifest["engines"])
    unknown = selected - set(manifest["engines"])
    if unknown:
        raise ValueError(f"selected engines are absent from the manifest: {sorted(unknown)}")
    reference = manifest["reference_engine"]
    if arguments.selected_engines and reference not in selected:
        raise ValueError(f"selected engines must include reference engine {reference!r}")
    required = {
        engine
        for engine, specification in manifest["engines"].items()
        if not specification.get("optional", False)
    }
    omitted = required - selected
    if omitted:
        raise ValueError(f"selected engines omit required engines: {sorted(omitted)}")
    return selected


def _collect_energies(
    arguments: argparse.Namespace,
    root: Path,
    manifest: dict[str, Any],
    checksum: str,
    required_terms: set[str],
) -> tuple[dict[str, EnergyReport], dict[str, str]]:
    selected = _selected_engines(arguments, manifest)
    energies: dict[str, EnergyReport] = {}
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
                command = list(specification["command"])
                if command and command[0] in {"python", "python3"}:
                    command[0] = sys.executable
                subprocess.run(command, cwd=directory, check=True)
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

    absolute = float(manifest["absolute_tolerance_kj_mol"])
    relative = float(manifest["relative_tolerance"])
    if not math.isfinite(absolute) or absolute < 0:
        raise ValueError("absolute tolerance must be finite and non-negative")
    if not math.isfinite(relative) or relative < 0:
        raise ValueError("relative tolerance must be finite and non-negative")

    reference_name = manifest["reference_engine"]
    energies, skipped = _collect_energies(args, root, manifest, expected_checksum, required_terms)

    if reference_name not in energies:
        raise ValueError(f"reference engine {reference_name!r} did not produce energies")
    for energy in energies.values():
        validate_potential_sum(energy, absolute_kj_mol=absolute, relative=relative)
    comparisons = {
        engine: compare_reports(
            energies[reference_name],
            values,
            absolute_kj_mol=absolute,
            relative=relative,
        )
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
        "energies": {engine: report.to_dict() for engine, report in energies.items()},
        "comparisons": comparisons,
        "skipped_optional_engines": skipped,
        "passed": passed,
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"single-point cross-engine validation: {'PASS' if passed else 'FAIL'}")
    print(f"report: {args.output}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
