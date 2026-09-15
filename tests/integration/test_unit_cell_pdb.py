"""Integration tests for unit-cell PDB generation through the CLI."""

from pathlib import Path

from cellulose_forge.cli import main


def test_cli_writes_annotated_i_beta_unit_cell(tmp_path: Path) -> None:
    output = tmp_path / "i_beta_unit_cell.pdb"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "unit-cell",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("HETATM")]
    assert len(atom_lines) == 84
    assert any(line.startswith("CRYST1") for line in lines)
    assert {line[21] for line in atom_lines} == {"A", "B"}
    assert {(line[21], int(line[22:26])) for line in atom_lines} == {
        ("A", 1),
        ("A", 2),
        ("B", 1),
        ("B", 2),
    }
