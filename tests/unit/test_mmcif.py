"""Tests for unrestricted-size PDBx/mmCIF coordinate output."""

from pathlib import Path

import numpy as np

from cellulose_forge.coordinates import write_coordinates
from cellulose_forge.geometry import BulkCrystal
from cellulose_forge.mmcif import write_mmcif
from cellulose_forge.pdb import MAXIMUM_PDB_ATOMS


def _crystal(atom_count: int) -> BulkCrystal:
    return BulkCrystal(
        allomorph="I-beta",
        space_group="P2_1",
        atom_labels=("C1",) * atom_count,
        cartesian_coordinates=np.zeros((atom_count, 3)),
        basis_vectors=np.diag([7.784, 8.201, 10.38]),
        replication=(1, 1, 1),
        chain_numbers=(1,) * atom_count,
        residue_numbers=(1,) * atom_count,
    )


def test_write_mmcif_preserves_cell_and_unrestricted_identifiers(tmp_path: Path) -> None:
    crystal = BulkCrystal(
        allomorph="I-beta",
        space_group="P2_1",
        atom_labels=("C1", "O4"),
        cartesian_coordinates=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        basis_vectors=np.diag([7.784, 8.201, 10.38]),
        replication=(1, 1, 1),
        chain_numbers=(1, 12_345),
        residue_numbers=(1, 123_456),
    )
    output = tmp_path / "structure.cif"

    write_mmcif(crystal, output)

    text = output.read_text(encoding="utf-8")
    assert "_cell.length_a 7.784000" in text
    assert "_symmetry.space_group_name_H-M 'P2 1'" in text
    assert "HETATM 2 O O4 . BGC C12345 1 123456 ? 4.000000 5.000000 6.000000" in text


def test_coordinate_writer_keeps_small_requested_pdb(tmp_path: Path) -> None:
    requested = tmp_path / "small.pdb"

    actual = write_coordinates(_crystal(2), requested)

    assert actual == requested
    assert requested.is_file()


def test_coordinate_writer_uses_mmcif_above_pdb_atom_limit(tmp_path: Path) -> None:
    requested = tmp_path / "large.pdb"

    actual = write_coordinates(_crystal(MAXIMUM_PDB_ATOMS + 1), requested)

    assert actual == tmp_path / "large.cif"
    assert actual.is_file()
    assert not requested.exists()
    text = actual.read_text(encoding="utf-8")
    assert f"HETATM {MAXIMUM_PDB_ATOMS + 1} " in text


def test_coordinate_writer_accepts_explicit_mmcif_extension(tmp_path: Path) -> None:
    requested = tmp_path / "structure.mmcif"

    actual = write_coordinates(_crystal(2), requested)

    assert actual == requested
    assert requested.read_text(encoding="utf-8").startswith("data_cellulose_forge_I-beta")
