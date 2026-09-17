"""PDBx/mmCIF coordinate output for large generated structures."""

from __future__ import annotations

import re
from pathlib import Path

from cellulose_forge.crystallography import basis_lengths_and_angles
from cellulose_forge.geometry import BulkCrystal


def _element_from_label(label: str) -> str:
    return label[0].upper()


def _data_block_name(allomorph: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]", "_", allomorph)
    return f"cellulose_forge_{cleaned}"


def write_mmcif(crystal: BulkCrystal, path: str | Path) -> None:
    """Write an unrestricted-size diagnostic PDBx/mmCIF coordinate file."""
    a, b, c, alpha, beta, gamma = basis_lengths_and_angles(crystal.basis_vectors)
    lines = [
        f"data_{_data_block_name(crystal.allomorph)}",
        "#",
        "_entry.id cellulose_forge",
        f"_cell.length_a {a:.6f}",
        f"_cell.length_b {b:.6f}",
        f"_cell.length_c {c:.6f}",
        f"_cell.angle_alpha {alpha:.6f}",
        f"_cell.angle_beta {beta:.6f}",
        f"_cell.angle_gamma {gamma:.6f}",
        f"_symmetry.space_group_name_H-M '{crystal.space_group.replace('_', ' ')}'",
    ]
    lines.extend(f"# cellulose_forge_annotation {annotation}" for annotation in crystal.annotations)
    lines.extend(
        (
            "#",
            "loop_",
            "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_atom_id",
            "_atom_site.pdbx_PDB_model_num",
        )
    )
    chain_numbers = crystal.chain_numbers or (1,) * crystal.atom_count
    residue_numbers = crystal.residue_numbers or (1,) * crystal.atom_count
    records = zip(
        crystal.atom_labels,
        crystal.cartesian_coordinates,
        chain_numbers,
        residue_numbers,
        strict=True,
    )
    for atom_id, (label, coords, chain_number, residue_number) in enumerate(records, start=1):
        x, y, z = coords
        asym_id = f"C{chain_number}"
        lines.append(
            f"HETATM {atom_id} {_element_from_label(label)} {label} . BGC "
            f"{asym_id} 1 {residue_number} ? {x:.6f} {y:.6f} {z:.6f} "
            f"1.00 0.00 {residue_number} BGC {asym_id} {label} 1"
        )
    lines.extend(("#", ""))
    Path(path).write_text("\n".join(lines), encoding="utf-8")
