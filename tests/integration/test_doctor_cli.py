"""Integration tests for the ``cellulose-forge doctor`` subcommand."""

from __future__ import annotations

import json

import pytest

from cellulose_forge.cli import main


def test_doctor_text_output_reports_python_and_numpy(capsys: pytest.CaptureFixture[str]) -> None:
    result = main(["doctor"])

    captured = capsys.readouterr()
    assert result in {0, 1}
    assert "[OK] python" in captured.out
    assert "[OK] numpy" in captured.out


def test_doctor_json_output_is_parseable(capsys: pytest.CaptureFixture[str]) -> None:
    main(["doctor", "--json"])

    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert isinstance(payload["coordinate_only_ready"], bool)


def test_doctor_reports_missing_tleap_with_nonzero_exit(
    capsys: pytest.CaptureFixture[str],
) -> None:
    result = main(["doctor", "--tleap-executable", "definitely-not-a-real-tleap-binary"])

    captured = capsys.readouterr()
    assert result == 1
    assert "[MISSING] tleap" in captured.out


def test_doctor_json_survives_a_missing_parmed_interpreter(
    capsys: pytest.CaptureFixture[str],
) -> None:
    result = main(
        [
            "doctor",
            "--json",
            "--parmed-python",
            "/definitely/not/a/real/python-interpreter",
        ]
    )

    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    parmed = next(check for check in payload["checks"] if check["name"] == "parmed")
    assert result == 1
    assert parmed["ok"] is False
    assert payload["coordinate_only_ready"] is True
    assert payload["glycam_gromacs_files_ready"] is False
