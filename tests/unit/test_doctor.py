"""Tests for environment diagnostics."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from cellulose_forge.doctor import render_json, render_text, run_diagnostics


def test_python_and_numpy_checks_pass_in_the_test_environment() -> None:
    checks = run_diagnostics()
    by_name = {check.name: check for check in checks}

    assert by_name["python"].ok
    assert by_name["numpy"].ok
    assert by_name["cellulose-forge"].ok


def test_missing_tleap_executable_is_reported_with_a_fix() -> None:
    checks = run_diagnostics(tleap_executable="definitely-not-a-real-tleap-binary")
    by_name = {check.name: check for check in checks}

    assert by_name["tleap"].ok is False
    assert by_name["tleap"].fix is not None


def test_parmed_check_uses_the_current_interpreter_when_available() -> None:
    checks = run_diagnostics(parmed_python=sys.executable)
    by_name = {check.name: check for check in checks}

    assert by_name["parmed"].name == "parmed"
    assert by_name["parmed"].detail.endswith(sys.executable)


def test_missing_parmed_interpreter_is_reported_instead_of_raising() -> None:
    missing = "/definitely/not/a/real/python-interpreter"

    checks = run_diagnostics(parmed_python=missing)
    by_name = {check.name: check for check in checks}

    assert by_name["parmed"].ok is False
    assert missing in by_name["parmed"].detail
    assert by_name["parmed"].fix is not None


def test_missing_openmm_interpreter_is_reported_instead_of_raising() -> None:
    missing = "/definitely/not/a/real/openmm-python-interpreter"

    checks = run_diagnostics(openmm_python=missing)
    by_name = {check.name: check for check in checks}

    assert by_name["openmm"].ok is False
    assert missing in by_name["openmm"].detail
    assert by_name["openmm"].fix is not None


def test_render_text_reports_readiness_tiers() -> None:
    checks = run_diagnostics(tleap_executable="definitely-not-a-real-tleap-binary")

    text = render_text(checks)

    assert "coordinate-only builds ready: True" in text
    assert "GLYCAM/GROMACS-file builds ready: False" in text


def test_render_json_is_valid_and_matches_check_count() -> None:
    checks = run_diagnostics()

    payload = json.loads(render_json(checks))

    assert len(payload["checks"]) == len(checks)
    assert set(payload) >= {
        "checks",
        "coordinate_only_ready",
        "glycam_gromacs_files_ready",
        "charmm36_gromacs_files_ready",
        "gromacs_simulation_ready",
        "glycam_openmm_ready",
        "charmm36_openmm_ready",
    }


def test_charmm36_openmm_readiness_does_not_require_glycam_tools(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(
        "cellulose_forge.doctor.discover_charmm36_directory", lambda: tmp_path
    )
    checks = run_diagnostics(
        tleap_executable="/definitely/not/tleap", openmm_python=sys.executable
    )
    payload = json.loads(render_json(checks))
    openmm_ok = {check.name: check.ok for check in checks}["openmm"]

    assert payload["glycam_gromacs_files_ready"] is False
    assert payload["glycam_openmm_ready"] is False
    assert payload["charmm36_gromacs_files_ready"] is True
    assert payload["charmm36_openmm_ready"] is openmm_ok


def test_missing_charmm36_release_is_reported_with_install_command(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def missing() -> Path:
        raise FileNotFoundError("CHARMM36 is not installed")

    monkeypatch.setattr("cellulose_forge.doctor.discover_charmm36_directory", missing)
    by_name = {check.name: check for check in run_diagnostics()}

    assert by_name["charmm36"].ok is False
    assert by_name["charmm36"].fix == "cellulose-forge forcefield install charmm36"
