"""Tests for the matched force-field validation driver."""

from pathlib import Path
from runpy import run_path
from typing import Any, cast

import pytest

_SCRIPT = Path(__file__).resolve().parents[2] / "validation" / "compare_force_fields.py"
_NAMESPACE = run_path(str(_SCRIPT), run_name="validation_compare_force_fields")
_parse_xvg = cast(Any, _NAMESPACE["_parse_xvg"])
_summarize_rows = cast(Any, _NAMESPACE["_summarize_rows"])
_write_protocol = cast(Any, _NAMESPACE["_write_protocol"])


def test_write_protocol_replaces_only_nsteps(tmp_path: Path) -> None:
    source = tmp_path / "source.mdp"
    destination = tmp_path / "destination.mdp"
    source.write_text("integrator = md\nnsteps = 10\ndt = 0.001\n", encoding="utf-8")

    _write_protocol(source, destination, nsteps=25)

    assert destination.read_text(encoding="utf-8") == (
        "integrator = md\nnsteps                  = 25\ndt = 0.001\n"
    )


def test_parse_and_summarize_gromacs_observables(tmp_path: Path) -> None:
    table = tmp_path / "energy.xvg"
    table.write_text(
        "# comment\n@ metadata\n"
        "0 10 300 1 3 4 5 60 1500\n"
        "1 14 302 3 5 6 7 70 1400\n",
        encoding="utf-8",
    )

    summary = _summarize_rows(_parse_xvg(table))

    assert summary["frames"] == 2
    assert summary["Potential"]["mean"] == pytest.approx(12.0)
    assert summary["Temperature"]["mean"] == pytest.approx(301.0)
    assert summary["Density"]["final_minus_initial"] == pytest.approx(-100.0)


def test_parse_rejects_incomplete_observable_table(tmp_path: Path) -> None:
    table = tmp_path / "energy.xvg"
    table.write_text("0 1 2\n", encoding="utf-8")

    with pytest.raises(ValueError, match="unexpected GROMACS energy table"):
        _parse_xvg(table)
