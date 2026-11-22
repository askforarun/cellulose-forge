"""Coordinate-format selection for generated structures."""

from __future__ import annotations

from pathlib import Path

from cellulose_forge.geometry import BulkCrystal
from cellulose_forge.glycam_mapping import GlycamMappedSystem
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


def write_mapped_coordinates(
    mapped: GlycamMappedSystem,
    crystal: BulkCrystal,
    requested_path: str | Path,
) -> Path:
    """Write chemically complete mapped coordinates with crystal metadata."""
    chain_id_by_number = dict(zip(crystal.chain_numbers, crystal.chain_ids, strict=True))
    chain_ids = tuple(chain_id_by_number.get(number, "A") for number in mapped.chain_numbers)
    cap_annotation: tuple[str, ...]
    atom_names = set(mapped.atom_names)
    if {"O1", "HO1", "HO4"} <= atom_names:
        cap_annotation = ("finite-chain terminal caps: O1/HO1 and HO4",)
    elif {"O1", "HO1", "H4O"} <= atom_names:
        cap_annotation = ("finite-chain terminal caps: O1/HO1 and O4/H4O (GLYCAM naming)",)
    else:
        cap_annotation = ()
    structure = BulkCrystal(
        allomorph=crystal.allomorph,
        space_group=crystal.space_group,
        atom_labels=mapped.atom_names,
        cartesian_coordinates=mapped.cartesian_coordinates,
        basis_vectors=mapped.basis_vectors,
        replication=crystal.replication,
        chain_ids=chain_ids,
        residue_numbers=mapped.residue_numbers,
        chain_numbers=mapped.chain_numbers,
        annotations=(
            *crystal.annotations,
            "chemically complete mapped coordinates",
            *cap_annotation,
            "PDB residue name uses the three-character crystallographic alias",
        ),
        residue_name=crystal.residue_name,
    )
    return write_coordinates(structure, requested_path)
