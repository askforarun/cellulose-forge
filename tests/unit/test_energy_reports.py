from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from cellulose_forge.energy_drivers.report import compare_reports, load_energy_report


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
