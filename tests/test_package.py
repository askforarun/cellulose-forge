"""Package-foundation tests."""

from cellulose_forge import __version__
from cellulose_forge.cli import build_parser, main


def test_version_is_defined() -> None:
    assert __version__


def test_cli_has_expected_program_name() -> None:
    assert build_parser().prog == "cellulose-forge"


def test_empty_cli_succeeds() -> None:
    assert main([]) == 0
