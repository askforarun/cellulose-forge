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


def test_cli_writes_i_beta_bulk_supercell(tmp_path: Path) -> None:
    output = tmp_path / "i_beta_2x2x2.pdb"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "bulk",
            "--cells",
            "2",
            "2",
            "2",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("HETATM")]
    cryst1 = next(line for line in lines if line.startswith("CRYST1"))
    assert len(atom_lines) == 84 * 2 * 2 * 2
    assert (float(cryst1[6:15]), float(cryst1[15:24]), float(cryst1[24:33])) == (
        15.568,
        16.402,
        20.760,
    )
    assert {line[21] for line in atom_lines} == set("ABCDEFGH")
    assert {int(line[22:26]) for line in atom_lines} == {1, 2, 3, 4}


def test_cli_writes_i_alpha_unit_cell(tmp_path: Path) -> None:
    output = tmp_path / "i_alpha_unit_cell.pdb"

    result = main(
        [
            "build",
            "--allomorph",
            "I-alpha",
            "--shape",
            "unit-cell",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("HETATM")]
    cryst1 = next(line for line in lines if line.startswith("CRYST1"))
    assert len(atom_lines) == 42
    assert {line[21] for line in atom_lines} == {"A"}
    assert {int(line[22:26]) for line in atom_lines} == {1, 2}
    assert (float(cryst1[6:15]), float(cryst1[15:24]), float(cryst1[24:33])) == (
        10.400,
        6.717,
        5.962,
    )


def test_cli_writes_iii_i_unit_cell(tmp_path: Path) -> None:
    output = tmp_path / "iii_i_unit_cell.pdb"

    assert (
        main(
            [
                "build",
                "--allomorph",
                "III-I",
                "--shape",
                "unit-cell",
                "--output",
                str(output),
            ]
        )
        == 0
    )
    lines = output.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("HETATM")]
    assert len(atom_lines) == 42
    assert {line[21] for line in atom_lines} == {"A"}
    assert {int(line[22:26]) for line in atom_lines} == {1, 2}


def test_cli_writes_ii_unit_cell_without_inventing_hydrogens(tmp_path: Path) -> None:
    output = tmp_path / "ii_unit_cell.pdb"

    assert (
        main(
            [
                "build",
                "--allomorph",
                "II",
                "--shape",
                "unit-cell",
                "--output",
                str(output),
            ]
        )
        == 0
    )
    lines = output.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("HETATM")]
    assert len(atom_lines) == 72
    assert {line[21] for line in atom_lines} == {"A", "B"}
    assert not any(line[12:16].strip() in {"HO2", "HO3", "HO6"} for line in atom_lines)


def test_cli_writes_i_beta_110_slab(tmp_path: Path) -> None:
    output = tmp_path / "i_beta_110_slab.pdb"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "slab",
            "--face",
            "1",
            "1",
            "0",
            "--layers",
            "4",
            "--chains-per-layer",
            "20",
            "--degree-polymerization",
            "12",
            "--chain-boundary",
            "finite",
            "--vacuum",
            "30",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("HETATM")]
    assert len(atom_lines) == 20_160
    assert any("shape=slab face=(1,1,0) layers=4" in line for line in lines)
    assert any("chains_per_layer=20 total_chains=80" in line for line in lines)
    assert any("degree_polymerization=12" in line for line in lines)
    assert any("vacuum=30.000000 angstrom" in line for line in lines)
    assert atom_lines[0][72:76] == "C001"
    assert {line[72:76] for line in atom_lines} == {f"C{chain:03d}" for chain in range(1, 81)}


def test_cli_writes_slab_from_physical_dimensions(tmp_path: Path) -> None:
    output = tmp_path / "i_beta_110_sized_slab.pdb"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "slab",
            "--face",
            "1",
            "1",
            "0",
            "--slab-width",
            "40",
            "--slab-thickness",
            "18",
            "--slab-length",
            "25",
            "--chain-boundary",
            "periodic",
            "--vacuum",
            "10",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    assert any("size_mode=physical_dimensions" in line for line in lines)
    assert any("requested_width=40.000000" in line for line in lines)
    assert any("requested_thickness=18.000000" in line for line in lines)
    assert any("requested_length=25.000000" in line for line in lines)
    assert any("chain_boundary=periodic" in line for line in lines)


def test_cli_writes_literature_i_beta_36_chain_fibril(tmp_path: Path) -> None:
    output = tmp_path / "i_beta_36_chain_fibril.pdb"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "fibril",
            "--fibril-model",
            "36-chain",
            "--degree-polymerization",
            "6",
            "--chain-boundary",
            "periodic",
            "--vacuum",
            "20",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("HETATM")]
    assert len(atom_lines) == 36 * 6 * 21
    assert {line[72:76] for line in atom_lines} == {
        f"C{chain:03d}" for chain in range(1, 37)
    }
    assert any("shape=fibril fibril_model=36-chain total_chains=36" in line for line in lines)
