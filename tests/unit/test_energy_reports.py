from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from cellulose_forge.energy_drivers.report import (
    EnergyReport,
    compare_reports,
    load_energy_report,
    validate_potential_sum,
)


def _data(engine: str, terms: dict[str, float]) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "engine": engine,
        "units": "kJ/mol",
        "coordinate_sha256": "abc",
        "terms": terms,
        "corrections": [],
    }


def _load(tmp_path: Path, engine: str, terms: dict[str, float]):
    path = tmp_path / f"{engine}.json"
    path.write_text(json.dumps(_data(engine, terms)), encoding="utf-8")
    return load_energy_report(
        path,
        expected_engine=engine,
        coordinate_sha256="abc",
        required_terms={"bond", "potential"},
    )


def test_report_round_trip_preserves_schema_one_data(tmp_path: Path) -> None:
    source = _data("gromacs", {"bond": 1.0, "potential": -2.0})
    source["engine_version"] = "2025.4"
    path = tmp_path / "energy.json"
    path.write_text(json.dumps(source), encoding="utf-8")

    report = load_energy_report(
        path,
        expected_engine="gromacs",
        coordinate_sha256="abc",
        required_terms={"bond", "potential"},
    )

    assert report.to_dict() == source


def test_round_trip_preserves_integer_terms_and_absent_corrections(tmp_path: Path) -> None:
    source = _data("gromacs", {"bond": 1, "potential": 1})
    del source["corrections"]
    path = tmp_path / "energy.json"
    path.write_text(json.dumps(source), encoding="utf-8")

    report = load_energy_report(
        path,
        expected_engine="gromacs",
        coordinate_sha256="abc",
        required_terms={"bond", "potential"},
    )

    assert report.to_dict() == source
    assert type(report.to_dict()["terms"]["bond"]) is int


def test_comparison_matches_mixed_absolute_relative_rule(tmp_path: Path) -> None:
    reference = _load(tmp_path, "gromacs", {"bond": 1000.0, "potential": -2.0})
    candidate = _load(tmp_path, "openmm", {"bond": 1000.05, "potential": -1.995})

    result = compare_reports(reference, candidate, absolute_kj_mol=0.01, relative=1.0e-4)

    assert result["passed"] is True
    assert result["common_terms"]["bond"]["allowed_difference_kj_mol"] == 0.100005
    assert result["common_terms"]["potential"]["allowed_difference_kj_mol"] == 0.01


@pytest.mark.parametrize(
    ("field", "value", "error"),
    [
        ("units", "kcal/mol", "normalized"),
        ("coordinate_sha256", "wrong", "coordinate checksum"),
        ("terms", {"bond": True, "potential": 1.0}, "not numeric"),
        ("terms", {"bond": float("nan"), "potential": 1.0}, "not finite"),
        ("terms", {"unknown": 1.0, "potential": 1.0}, "unknown canonical"),
    ],
)
def test_invalid_reports_are_refused(tmp_path: Path, field: str, value: object, error: str) -> None:
    data = _data("gromacs", {"bond": 1.0, "potential": 1.0})
    data[field] = value
    path = tmp_path / "energy.json"
    path.write_text(json.dumps(data), encoding="utf-8")

    with pytest.raises((TypeError, ValueError), match=error):
        load_energy_report(
            path,
            expected_engine="gromacs",
            coordinate_sha256="abc",
            required_terms={"bond", "potential"},
        )


def test_comparison_refuses_different_coordinates(tmp_path: Path) -> None:
    reference = _load(tmp_path, "gromacs", {"bond": 1.0, "potential": 1.0})
    path = tmp_path / "openmm.json"
    data = _data("openmm", {"bond": 1.0, "potential": 1.0})
    data["coordinate_sha256"] = "different"
    path.write_text(json.dumps(data), encoding="utf-8")
    candidate = load_energy_report(
        path,
        expected_engine="openmm",
        coordinate_sha256="different",
        required_terms={"bond", "potential"},
    )

    with pytest.raises(ValueError, match="different coordinates"):
        compare_reports(reference, candidate, absolute_kj_mol=0.0, relative=0.0)


def test_direct_construction_is_validated_and_deeply_immutable() -> None:
    terms = {"bond": 1.0, "potential": 1.0}
    nested = {"notes": {"values": [1]}}
    report = EnergyReport("gromacs", "abc", terms, extra=nested)
    terms["bond"] = 8.0
    nested["notes"]["values"].append(2)

    assert report.to_dict()["terms"]["bond"] == 1.0
    assert report.to_dict()["notes"] == {"values": [1]}
    with pytest.raises(ValueError, match="collide"):
        EnergyReport(
            "gromacs",
            "abc",
            {"potential": 1.0},
            extra={"engine": "wrong"},
        )
    with pytest.raises(ValueError, match="schema"):
        EnergyReport("gromacs", "abc", {"potential": 1.0}, schema_version=True)


def test_potential_must_equal_component_sum() -> None:
    report = EnergyReport("gromacs", "abc", {"bond": 1.0, "potential": 2.0})

    with pytest.raises(ValueError, match="component sum"):
        validate_potential_sum(report, absolute_kj_mol=0.0, relative=0.0)


def test_unmatched_terms_fail_comparison() -> None:
    reference = EnergyReport("gromacs", "abc", {"bond": 1.0, "potential": 1.0})
    candidate = EnergyReport("openmm", "abc", {"bond": 1.0, "improper": 2.0, "potential": 3.0})

    result = compare_reports(reference, candidate, absolute_kj_mol=0.0, relative=0.0)

    assert result["passed"] is False
    assert result["unmatched_terms"]["candidate_only"] == ["improper"]


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -1.0])
def test_invalid_tolerances_are_refused(value: float) -> None:
    reference = EnergyReport("gromacs", "abc", {"potential": 1.0})
    candidate = EnergyReport("openmm", "abc", {"potential": 1.0})

    with pytest.raises(ValueError, match="finite and non-negative"):
        compare_reports(reference, candidate, absolute_kj_mol=value, relative=0.0)
