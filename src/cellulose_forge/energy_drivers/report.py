"""Typed, engine-neutral single-point energy reports and comparisons."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any

CANONICAL_TERMS = frozenset(
    {"bond", "angle", "dihedral", "improper", "vdw", "electrostatics", "potential"}
)


@dataclass(frozen=True)
class EnergyReport:
    """One engine's normalized single-point energies in kJ/mol."""

    engine: str
    coordinate_sha256: str
    terms: Mapping[str, float]
    corrections: tuple[Mapping[str, Any], ...]
    extra: Mapping[str, Any]
    schema_version: int = 1
    units: str = "kJ/mol"

    def to_dict(self) -> dict[str, Any]:
        """Return the migration-compatible JSON representation."""
        return {
            "schema_version": self.schema_version,
            "engine": self.engine,
            "units": self.units,
            "coordinate_sha256": self.coordinate_sha256,
            "terms": dict(self.terms),
            "corrections": [dict(item) for item in self.corrections],
            **dict(self.extra),
        }


def _validate_required_terms(required_terms: set[str]) -> None:
    unknown = required_terms - CANONICAL_TERMS
    if unknown:
        raise ValueError(f"unknown required canonical terms: {sorted(unknown)}")
    if "potential" not in required_terms:
        raise ValueError("required terms must include potential")


def energy_report_from_dict(
    data: Mapping[str, Any],
    *,
    expected_engine: str,
    coordinate_sha256: str,
    required_terms: set[str],
) -> EnergyReport:
    """Validate one schema-1 adapter result during the migration period."""
    _validate_required_terms(required_terms)
    if data.get("schema_version") != 1:
        raise ValueError(f"{expected_engine}: unsupported energy schema")
    if data.get("engine") != expected_engine:
        raise ValueError(f"{expected_engine}: report names engine {data.get('engine')!r}")
    if data.get("units") != "kJ/mol":
        raise ValueError(f"{expected_engine}: energies must be normalized to kJ/mol")
    if data.get("coordinate_sha256") != coordinate_sha256:
        raise ValueError(f"{expected_engine}: coordinate checksum differs from shared coordinates")
    raw_terms = data.get("terms")
    if not isinstance(raw_terms, dict) or not raw_terms:
        raise ValueError(f"{expected_engine}: report contains no terms")
    unknown = set(raw_terms) - CANONICAL_TERMS
    if unknown:
        raise ValueError(f"{expected_engine}: unknown canonical terms: {sorted(unknown)}")
    terms: dict[str, float] = {}
    for name, value in raw_terms.items():
        if isinstance(value, bool) or not isinstance(value, int | float):
            raise TypeError(f"{expected_engine}: term {name!r} is not numeric")
        number = float(value)
        if not math.isfinite(number):
            raise ValueError(f"{expected_engine}: term {name!r} is not finite")
        terms[name] = number
    missing = required_terms - set(terms)
    if missing:
        raise ValueError(
            f"{expected_engine}: required comparable terms are missing: {sorted(missing)}"
        )
    raw_corrections = data.get("corrections", [])
    if not isinstance(raw_corrections, list) or not all(
        isinstance(item, dict) for item in raw_corrections
    ):
        raise TypeError(f"{expected_engine}: corrections must be a list of objects")
    known = {"schema_version", "engine", "units", "coordinate_sha256", "terms", "corrections"}
    extra = {key: value for key, value in data.items() if key not in known}
    return EnergyReport(
        engine=expected_engine,
        coordinate_sha256=coordinate_sha256,
        terms=MappingProxyType(terms),
        corrections=tuple(MappingProxyType(dict(item)) for item in raw_corrections),
        extra=MappingProxyType(extra),
    )


def load_energy_report(
    path: Path,
    *,
    expected_engine: str,
    coordinate_sha256: str,
    required_terms: set[str],
) -> EnergyReport:
    """Load and validate one normalized engine report."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise TypeError(f"{expected_engine}: energy report must be an object")
    return energy_report_from_dict(
        data,
        expected_engine=expected_engine,
        coordinate_sha256=coordinate_sha256,
        required_terms=required_terms,
    )


def compare_reports(
    reference: EnergyReport,
    candidate: EnergyReport,
    *,
    absolute_kj_mol: float,
    relative: float,
) -> dict[str, Any]:
    """Compare common canonical terms with CelluloseForge's mixed tolerance."""
    if reference.coordinate_sha256 != candidate.coordinate_sha256:
        raise ValueError("engine reports use different coordinates")
    if absolute_kj_mol < 0 or relative < 0:
        raise ValueError("energy tolerances must be non-negative")
    common = sorted(set(reference.terms) & set(candidate.terms))
    if not common:
        raise ValueError("engines report no common energy terms")
    terms: dict[str, Any] = {}
    passed = True
    for name in common:
        expected = reference.terms[name]
        observed = candidate.terms[name]
        difference = observed - expected
        scale = max(abs(expected), abs(observed), 1.0)
        allowed = max(absolute_kj_mol, relative * scale)
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
