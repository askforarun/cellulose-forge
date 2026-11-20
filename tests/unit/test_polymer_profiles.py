"""Tests for declarative polymer capabilities and early CLI refusal."""

from __future__ import annotations

from pathlib import Path

import pytest

from cellulose_forge import cli
from cellulose_forge.polymer_profiles import (
    POLYMER_PROFILES,
    capability_error,
    engine_capability_error,
)


def test_capability_table_records_supported_force_fields_and_engines() -> None:
    assert set(POLYMER_PROFILES) == {"cellulose", "chitin", "mannan"}
    assert set(POLYMER_PROFILES["cellulose"].force_fields) == {
        "none",
        "glycam",
        "charmm36",
        "charmm-drude",
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
            "chitin fibrils are not implemented",
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
            "CHARMM36 does not support --md-engine amber",
        ),
    ],
)
def test_cli_refuses_unsupported_capability_before_construction(
    arguments: list[str],
    message: str,
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
) -> None:
    with pytest.raises(SystemExit):
        cli.main(["build", *arguments, "--output", str(tmp_path / "result")])
    assert message in capsys.readouterr().err


def test_native_amber_periodic_chains_are_refused_by_profile() -> None:
    assert capability_error("chitin", "glycam", "slab", ("amber",), "periodic") is None


@pytest.mark.parametrize(
    ("polymer", "force_field", "engines", "message"),
    [
        ("cellulose", "charmm36", ("amber",), "does not support --md-engine amber"),
        ("chitin", "glycam", ("namd",), "amber for finite chains"),
        ("cellulose", "openff", ("lammps",), "use gromacs or openmm only"),
        ("cellulose", "opls-aa", (), "support --md-engine gromacs, openmm, and lammps"),
    ],
)
def test_engine_refusals_name_the_request_and_preserve_qualifiers(
    polymer: str, force_field: str, engines: tuple[str, ...], message: str
) -> None:
    assert message in str(engine_capability_error(polymer, force_field, engines))


def test_unit_cell_reaches_physical_box_validation_instead_of_profile_refusal() -> None:
    assert capability_error("cellulose", "glycam", "unit-cell", ("gromacs",), "periodic") is None
    assert capability_error("chitin", "charmm36", "unit-cell", ("gromacs",), "periodic") is None
