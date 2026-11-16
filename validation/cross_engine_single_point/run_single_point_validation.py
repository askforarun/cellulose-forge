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

from cellulose_forge.energy_drivers import (
    OpenMMForceClassification,
    run_gromacs_single_point,
    run_namd_single_point,
    run_openmm_single_point,
)
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


def _path(directory: Path, value: object, description: str) -> Path:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{description} must be a non-empty path string")
    return (directory / value).resolve()


def _interaction_map(value: object, description: str) -> dict[int, tuple[int, ...]]:
    if not isinstance(value, dict):
        raise TypeError(f"{description} must be an object")
    result: dict[int, tuple[int, ...]] = {}
    for raw_force, raw_interactions in value.items():
        try:
            force = int(raw_force)
        except (TypeError, ValueError) as error:
            raise ValueError(f"{description} force indices must be integers") from error
        if str(force) != str(raw_force) or force < 0:
            raise ValueError(f"{description} force indices must be non-negative integers")
        if not isinstance(raw_interactions, list) or any(
            type(index) is not int or index < 0 for index in raw_interactions
        ):
            raise TypeError(f"{description} interaction indices must be non-negative integers")
        result[force] = tuple(raw_interactions)
    return result


def _force_roles(value: object) -> dict[int, str]:
    if not isinstance(value, dict):
        raise TypeError("OpenMM force_roles must be an object")
    result: dict[int, str] = {}
    for raw_force, role in value.items():
        try:
            force = int(raw_force)
        except (TypeError, ValueError) as error:
            raise ValueError("OpenMM force-role indices must be integers") from error
        if str(force) != str(raw_force) or force < 0 or not isinstance(role, str):
            raise ValueError(
                "OpenMM force roles require non-negative integer keys and string values"
            )
        result[force] = role
    return result


def _classified_terms(value: object) -> dict[str, str]:
    if value is None:
        return {}
    if not isinstance(value, dict) or any(
        not isinstance(name, str) or not isinstance(role, str) for name, role in value.items()
    ):
        raise TypeError("NAMD classified_terms must be a string-to-string object")
    return dict(value)


def _run_adapter(
    engine: str,
    specification: dict[str, Any],
    directory: Path,
    coordinate_sha256: str,
) -> EnergyReport:
    adapter = specification.get("adapter")
    inputs = specification.get("inputs")
    if not isinstance(inputs, dict):
        raise TypeError(f"{engine}: adapter inputs must be an object")
    if adapter == "gromacs":
        if engine != "gromacs":
            raise ValueError(f"{engine}: gromacs adapter cannot write another engine's report")
        report = run_gromacs_single_point(
            directory,
            coordinate_file=_path(directory, inputs.get("coordinate_file"), "coordinate_file"),
            topology_file=_path(directory, inputs.get("topology_file"), "topology_file"),
            mdp_file=_path(directory, inputs.get("mdp_file"), "mdp_file"),
            coordinate_sha256=coordinate_sha256,
            executable=str(inputs.get("executable", "gmx")),
        )
    elif adapter == "openmm":
        if engine != "openmm":
            raise ValueError(f"{engine}: openmm adapter cannot write another engine's report")
        raw_classification = inputs.get("classification")
        if not isinstance(raw_classification, dict):
            raise TypeError("openmm: classification must be an object")
        classification = OpenMMForceClassification(
            urey_bradley_bonds=_interaction_map(
                raw_classification.get("urey_bradley_bonds"), "OpenMM Urey--Bradley map"
            ),
            periodic_impropers=_interaction_map(
                raw_classification.get("periodic_impropers"), "OpenMM periodic-improper map"
            ),
            rb_impropers=_interaction_map(
                raw_classification.get("rb_impropers"), "OpenMM RB-improper map"
            ),
            force_roles=_force_roles(raw_classification.get("force_roles")),
        )
        report = run_openmm_single_point(
            system_file=_path(directory, inputs.get("system_file"), "system_file"),
            coordinate_file=_path(directory, inputs.get("coordinate_file"), "coordinate_file"),
            coordinate_sha256=coordinate_sha256,
            classification=classification,
            platform_name=str(inputs.get("platform_name", "Reference")),
        )
    elif adapter == "namd":
        if engine != "namd":
            raise ValueError(f"{engine}: namd adapter cannot write another engine's report")
        threads = inputs.get("threads", 1)
        if type(threads) is not int:
            raise TypeError("namd: threads must be an integer")
        report = run_namd_single_point(
            directory,
            configuration_file=_path(
                directory, inputs.get("configuration_file"), "configuration_file"
            ),
            shared_coordinate_file=_path(
                directory, inputs.get("coordinate_file"), "coordinate_file"
            ),
            binary_coordinate_file=_path(
                directory, inputs.get("binary_coordinate_file"), "binary_coordinate_file"
            ),
            coordinate_sha256=coordinate_sha256,
            executable=str(inputs.get("executable", "namd3")),
            threads=threads,
            classified_terms=_classified_terms(inputs.get("classified_terms")),
        )
    else:
        raise ValueError(f"{engine}: unsupported built-in energy adapter {adapter!r}")
    (directory / "energy.json").write_text(
        json.dumps(report.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return report


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
                if "adapter" in specification:
                    _run_adapter(engine, specification, directory, checksum)
                else:
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
