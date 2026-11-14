"""Tests for declarative polymer capabilities and early CLI refusal."""

from __future__ import annotations

from pathlib import Path

import pytest

from cellulose_forge import cli
from cellulose_forge.polymer_profiles import POLYMER_PROFILES, capability_error


def test_capability_table_records_supported_force_fields_and_engines() -> None:
    assert set(POLYMER_PROFILES) == {"cellulose", "chitin", "mannan"}
    assert set(POLYMER_PROFILES["cellulose"].force_fields) == {
        "none",
        "glycam",
        "charmm36",
        "openff",
        "opls-aa",
        "opls-aa-sei",
    }
    assert set(POLYMER_PROFILES["chitin"].force_fields) == {
        "none",
        "glycam",
        "charmm36",
        "openff",
    }
    assert set(POLYMER_PROFILES["mannan"].force_fields) == {"none"}
    assert POLYMER_PROFILES["chitin"].force_fields["openff"].engines == {
        "gromacs",
        "openmm",
    }


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        (
            [
                "--polymer",
                "chitin",
                "--allomorph",
                "beta",
                "--shape",
                "bulk",
                "--force-field",
                "opls-aa",
            ],
            "no published OPLS GlcNAc parameters",
        ),
        (
            [
                "--polymer",
                "chitin",
                "--allomorph",
                "beta",
                "--shape",
                "fibril",
                "--force-field",
                "openff",
            ],
            "supports --shape bulk, crystal, slab",
        ),
        (
            [
                "--polymer",
                "mannan",
                "--allomorph",
                "I",
                "--shape",
                "bulk",
                "--force-field",
                "glycam",
            ],
            "mannan force-field output is not implemented",
        ),
        (
            [
                "--allomorph",
                "I-beta",
                "--shape",
                "bulk",
                "--force-field",
                "charmm36",
                "--md-engine",
                "amber",
            ],
            "supports --md-engine gromacs, openmm, and lammps",
        ),
    ],
)
def test_cli_refuses_unsupported_capability_before_construction(
    arguments: list[str],
    message: str,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
) -> None:
    def construction_must_not_start(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("crystal construction started before capability validation")

    monkeypatch.setattr(cli, "_load_definition", construction_must_not_start)
    with pytest.raises(SystemExit):
        cli.main(["build", *arguments, "--output", str(tmp_path / "result")])
    assert message in capsys.readouterr().err


def test_native_amber_periodic_chains_are_refused_by_profile() -> None:
    assert (
        capability_error("chitin", "glycam", "slab", ("amber",), "periodic")
        == "native Amber output is available only for finite chains"
    )
