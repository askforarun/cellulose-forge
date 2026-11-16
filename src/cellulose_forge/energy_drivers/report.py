"""Typed, engine-neutral single-point energy reports and comparisons."""

from __future__ import annotations

import copy
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any

CANONICAL_TERMS = frozenset(
    {
        "bond",
        "angle",
        "dihedral",
        "improper",
        "vdw",
        "electrostatics",
        "restraint",
        "dispersion_correction",
        "potential",
    }
)
_RESERVED_FIELDS = frozenset(
    {"schema_version", "engine", "units", "coordinate_sha256", "terms", "corrections"}
)
EnergyValue = int | float


def _finite_number(value: object, description: str) -> EnergyValue:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise TypeError(f"{description} is not numeric")
    try:
        finite = math.isfinite(value)
    except OverflowError as error:
        raise ValueError(f"{description} is outside the supported numeric range") from error
    if not finite:
        raise ValueError(f"{description} is not finite")
    return value


def validate_tolerance(value: float, description: str) -> None:
    """Require a finite, non-negative energy tolerance."""
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"{description} must be finite and non-negative")


def _freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType({key: _freeze(copy.deepcopy(item)) for key, item in value.items()})
    if isinstance(value, list | tuple):
        return tuple(_freeze(copy.deepcopy(item)) for item in value)
    return copy.deepcopy(value)


def _thaw(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _thaw(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_thaw(item) for item in value]
    return copy.deepcopy(value)


@dataclass(frozen=True)
class EnergyReport:
    """One engine's normalized single-point energies in kJ/mol."""

    engine: str
    coordinate_sha256: str
    terms: Mapping[str, EnergyValue]
    corrections: tuple[Mapping[str, Any], ...] | None = None
    # A factory is required for Python 3.11, whose dataclasses reject a
    # MappingProxyType instance as a field default even though it is immutable.
    extra: Mapping[str, Any] = field(default_factory=lambda: MappingProxyType({}))
    schema_version: int = 1
    units: str = "kJ/mol"

    def __post_init__(self) -> None:
        if type(self.schema_version) is not int or self.schema_version != 1:
            raise ValueError("unsupported energy schema")
        if not isinstance(self.engine, str) or not self.engine:
            raise ValueError("engine is required")
        if self.units != "kJ/mol":
            raise ValueError("energies must be normalized to kJ/mol")
        if not isinstance(self.coordinate_sha256, str) or not self.coordinate_sha256:
            raise ValueError("coordinate checksum is required")
        if not isinstance(self.terms, Mapping) or not self.terms:
            raise ValueError(f"{self.engine}: report contains no terms")
        unknown = set(self.terms) - CANONICAL_TERMS
        if unknown:
            raise ValueError(f"{self.engine}: unknown canonical terms: {sorted(unknown)}")
        if "potential" not in self.terms:
            raise ValueError(f"{self.engine}: potential energy is required")
        terms = {
            name: _finite_number(value, f"{self.engine}: term {name!r}")
            for name, value in self.terms.items()
        }
        if self.corrections is not None and not all(
            isinstance(item, Mapping) for item in self.corrections
        ):
            raise TypeError(f"{self.engine}: corrections must be a sequence of objects")
        overlap = set(self.extra) & _RESERVED_FIELDS
        if overlap:
            raise ValueError(f"extra fields collide with canonical fields: {sorted(overlap)}")
        object.__setattr__(self, "terms", MappingProxyType(terms))
        object.__setattr__(
            self,
            "corrections",
            None if self.corrections is None else tuple(_freeze(item) for item in self.corrections),
        )
        object.__setattr__(self, "extra", _freeze(self.extra))

    def to_dict(self) -> dict[str, Any]:
        """Return the byte-compatible schema-1 JSON representation."""
        result = {
            "schema_version": self.schema_version,
            "engine": self.engine,
            "units": self.units,
            "coordinate_sha256": self.coordinate_sha256,
            "terms": dict(self.terms),
        }
        if self.corrections is not None:
            result["corrections"] = _thaw(self.corrections)
        result.update(_thaw(self.extra))
        return result


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
    if type(data.get("schema_version")) is not int or data.get("schema_version") != 1:
        raise ValueError(f"{expected_engine}: unsupported energy schema")
    if data.get("engine") != expected_engine:
        raise ValueError(f"{expected_engine}: report names engine {data.get('engine')!r}")
    if data.get("coordinate_sha256") != coordinate_sha256:
        raise ValueError(f"{expected_engine}: coordinate checksum differs from shared coordinates")
    raw_terms = data.get("terms")
    if not isinstance(raw_terms, dict):
        raise TypeError(f"{expected_engine}: terms must be an object")
    if "potential" not in raw_terms:
        raise ValueError(f"{expected_engine}: potential energy is required")
    raw_corrections = data.get("corrections") if "corrections" in data else None
    if raw_corrections is not None and (
        not isinstance(raw_corrections, list)
        or not all(isinstance(item, dict) for item in raw_corrections)
    ):
        raise TypeError(f"{expected_engine}: corrections must be a list of objects")
    extra = {key: value for key, value in data.items() if key not in _RESERVED_FIELDS}
    report = EnergyReport(
        engine=expected_engine,
        coordinate_sha256=coordinate_sha256,
        terms=raw_terms,
        corrections=None if raw_corrections is None else tuple(raw_corrections),
        extra=extra,
        schema_version=data["schema_version"],
        units=data.get("units", ""),
    )
    missing = required_terms - set(report.terms)
    if missing:
        raise ValueError(
            f"{expected_engine}: required comparable terms are missing: {sorted(missing)}"
        )
    return report


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


def validate_potential_sum(
    report: EnergyReport, *, absolute_kj_mol: float, relative: float
) -> None:
    """Require potential to equal all reported aggregate components."""
    validate_tolerance(absolute_kj_mol, "absolute tolerance")
    validate_tolerance(relative, "relative tolerance")
    components = sum(value for name, value in report.terms.items() if name != "potential")
    potential = report.terms["potential"]
    scale = max(abs(components), abs(potential), 1.0)
    allowed = max(absolute_kj_mol, relative * scale)
    if abs(potential - components) > allowed:
        raise ValueError(
            f"{report.engine}: reported potential differs from its component sum by "
            f"{potential - components:.6g} kJ/mol"
        )


def compare_reports(
    reference: EnergyReport,
    candidate: EnergyReport,
    *,
    absolute_kj_mol: float,
    relative: float,
    term_tolerances: Mapping[str, tuple[float, float]] | None = None,
) -> dict[str, Any]:
    """Compare canonical terms with mixed global or per-term tolerances."""
    if reference.coordinate_sha256 != candidate.coordinate_sha256:
        raise ValueError("engine reports use different coordinates")
    validate_tolerance(absolute_kj_mol, "absolute tolerance")
    validate_tolerance(relative, "relative tolerance")
    term_tolerances = term_tolerances or {}
    unknown_tolerances = set(term_tolerances) - CANONICAL_TERMS
    if unknown_tolerances:
        raise ValueError(f"tolerances name unknown terms: {sorted(unknown_tolerances)}")
    for name, (absolute, relative_term) in term_tolerances.items():
        validate_tolerance(absolute, f"{name} absolute tolerance")
        validate_tolerance(relative_term, f"{name} relative tolerance")
    reference_only = sorted(set(reference.terms) - set(candidate.terms))
    candidate_only = sorted(set(candidate.terms) - set(reference.terms))
    common = sorted(set(reference.terms) & set(candidate.terms))
    if not common:
        raise ValueError("engines report no common energy terms")
    terms: dict[str, Any] = {}
    passed = not reference_only and not candidate_only
    for name in common:
        # Comparison output historically used JSON floating-point values even
        # when an adapter supplied an exact integer (most commonly a zero
        # term).  Keep the source report lossless, but preserve that comparison
        # schema for byte-compatible migration results.
        expected = float(reference.terms[name])
        observed = float(candidate.terms[name])
        difference = observed - expected
        scale = max(abs(expected), abs(observed), 1.0)
        absolute, relative_term = term_tolerances.get(name, (absolute_kj_mol, relative))
        allowed = max(absolute, relative_term * scale)
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
    result: dict[str, Any] = {"passed": passed, "common_terms": terms}
    if reference_only or candidate_only:
        result["unmatched_terms"] = {
            "reference_only": reference_only,
            "candidate_only": candidate_only,
        }
    return result
