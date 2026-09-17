"""Coordinate-format selection for generated structures."""

from __future__ import annotations

from pathlib import Path

from cellulose_forge.geometry import BulkCrystal
from cellulose_forge.mmcif import write_mmcif
from cellulose_forge.pdb import MAXIMUM_PDB_ATOMS, write_pdb


def write_coordinates(crystal: BulkCrystal, requested_path: str | Path) -> Path:
    """Write PDB or mmCIF, automatically avoiding the PDB atom limit."""
    path = Path(requested_path)
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        write_mmcif(crystal, path)
        return path
    if suffix != ".pdb":
        msg = "coordinate output must use a .pdb, .cif, or .mmcif extension"
        raise ValueError(msg)
    if crystal.atom_count <= MAXIMUM_PDB_ATOMS:
        write_pdb(crystal, path)
        return path

    fallback = path.with_suffix(".cif")
    write_mmcif(crystal, fallback)
    return fallback
