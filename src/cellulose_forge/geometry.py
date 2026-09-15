"""Crystallographic unit-cell and bulk-crystal geometry generation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from cellulose_forge.crystallography import (
    AllomorphDefinition,
    FloatArray,
    fractional_to_cartesian,
)

_ARRAY_DIMENSIONS = 2
_COORDINATE_DIMENSIONS = 3
_PDB_CHAIN_IDENTIFIERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


@dataclass(frozen=True)
class BulkCrystal:
    """A generated crystallographic structure in Cartesian coordinates."""

    allomorph: str
    space_group: str
    atom_labels: tuple[str, ...]
    cartesian_coordinates: FloatArray
    basis_vectors: FloatArray
    replication: tuple[int, int, int]
    chain_ids: tuple[str, ...] = ()
    residue_numbers: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        coordinates = np.asarray(self.cartesian_coordinates, dtype=np.float64)
        basis = np.asarray(self.basis_vectors, dtype=np.float64)
        if coordinates.ndim != _ARRAY_DIMENSIONS or coordinates.shape[1] != _COORDINATE_DIMENSIONS:
            msg = f"cartesian_coordinates must have shape (n, 3), got {coordinates.shape}"
            raise ValueError(msg)
        if basis.shape != (3, 3):
            msg = f"basis_vectors must have shape (3, 3), got {basis.shape}"
            raise ValueError(msg)
        if len(self.atom_labels) != len(coordinates):
            msg = "atom labels and Cartesian coordinates must have equal lengths"
            raise ValueError(msg)
        if self.chain_ids and len(self.chain_ids) != len(coordinates):
            msg = "chain IDs and Cartesian coordinates must have equal lengths"
            raise ValueError(msg)
        if self.residue_numbers and len(self.residue_numbers) != len(coordinates):
            msg = "residue numbers and Cartesian coordinates must have equal lengths"
            raise ValueError(msg)
        if any(value < 1 for value in self.replication):
            msg = "replication values must be positive integers"
            raise ValueError(msg)
        coordinates.setflags(write=False)
        basis.setflags(write=False)
        object.__setattr__(self, "cartesian_coordinates", coordinates)
        object.__setattr__(self, "basis_vectors", basis)

    @property
    def atom_count(self) -> int:
        """Return the number of atoms in the generated geometry."""
        return len(self.atom_labels)


def _complete_cell_fractional(allomorph: AllomorphDefinition) -> FloatArray:
    asymmetric = allomorph.asymmetric_unit.fractional_coordinates
    return np.vstack([operation.apply(asymmetric) for operation in allomorph.symmetry_operations])


def build_unit_cell(allomorph: AllomorphDefinition) -> BulkCrystal:
    """Build one complete crystallographic unit cell from its asymmetric unit."""
    asymmetric = allomorph.asymmetric_unit
    fractional = _complete_cell_fractional(allomorph)
    operation_count = len(allomorph.symmetry_operations)
    atom_labels = asymmetric.atom_labels * operation_count

    chain_names = list(dict.fromkeys(asymmetric.chain_labels))
    if len(chain_names) > len(_PDB_CHAIN_IDENTIFIERS):
        msg = "diagnostic PDB output supports at most 62 crystallographic chains"
        raise ValueError(msg)
    chain_letters = {name: _PDB_CHAIN_IDENTIFIERS[index] for index, name in enumerate(chain_names)}

    chain_ids: list[str] = []
    residue_numbers: list[int] = []
    residues_per_asymmetric_unit = len(asymmetric.chain_labels)
    for operation_index in range(operation_count):
        for residue_index in range(residues_per_asymmetric_unit):
            chain_name = asymmetric.chain_labels[residue_index]
            chain_ids.extend([chain_letters[chain_name]] * asymmetric.residue_atom_count)
            residue_numbers.extend([operation_index + 1] * asymmetric.residue_atom_count)

    return BulkCrystal(
        allomorph=allomorph.name,
        space_group=allomorph.space_group,
        atom_labels=atom_labels,
        cartesian_coordinates=fractional_to_cartesian(fractional, allomorph.unit_cell),
        basis_vectors=allomorph.unit_cell.vectors(),
        replication=(1, 1, 1),
        chain_ids=tuple(chain_ids),
        residue_numbers=tuple(residue_numbers),
    )


def build_bulk_crystal(
    allomorph: AllomorphDefinition, replication: tuple[int, int, int]
) -> BulkCrystal:
    """Replicate a complete unit cell into a periodic bulk supercell.

    Replication follows the crystallographic ``a``, ``b``, and ``c`` basis
    vectors. Chains translated along ``c`` retain their chain identifier and
    receive consecutive residue numbers. Transverse translations create new
    chains.
    """
    if len(replication) != _COORDINATE_DIMENSIONS or any(
        not isinstance(value, int) or isinstance(value, bool) or value < 1 for value in replication
    ):
        msg = "replication must contain three positive integers"
        raise ValueError(msg)

    nx, ny, nz = replication
    chains_per_cell = len(dict.fromkeys(allomorph.asymmetric_unit.chain_labels))
    transverse_chain_count = nx * ny * chains_per_cell
    if transverse_chain_count > len(_PDB_CHAIN_IDENTIFIERS):
        msg = (
            "diagnostic PDB chain identifiers support at most 62 transverse "
            f"chains; requested {transverse_chain_count}"
        )
        raise ValueError(msg)

    unit_cell = build_unit_cell(allomorph)
    unit_vectors = allomorph.unit_cell.vectors()
    unit_chain_order = tuple(dict.fromkeys(unit_cell.chain_ids))
    unit_chain_indices = {chain_id: index for index, chain_id in enumerate(unit_chain_order)}
    residues_per_chain_per_cell = max(unit_cell.residue_numbers)

    coordinate_blocks: list[FloatArray] = []
    atom_labels: list[str] = []
    chain_ids: list[str] = []
    residue_numbers: list[int] = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                translation = i * unit_vectors[0] + j * unit_vectors[1] + k * unit_vectors[2]
                coordinate_blocks.append(unit_cell.cartesian_coordinates + translation)
                atom_labels.extend(unit_cell.atom_labels)

                transverse_cell_index = j * nx + i
                for chain_id, residue_number in zip(
                    unit_cell.chain_ids, unit_cell.residue_numbers, strict=True
                ):
                    chain_index = (
                        transverse_cell_index * chains_per_cell + unit_chain_indices[chain_id]
                    )
                    chain_ids.append(_PDB_CHAIN_IDENTIFIERS[chain_index])
                    residue_numbers.append(residue_number + k * residues_per_chain_per_cell)

    scale = np.asarray(replication, dtype=np.float64)[:, np.newaxis]
    return BulkCrystal(
        allomorph=allomorph.name,
        space_group=allomorph.space_group,
        atom_labels=tuple(atom_labels),
        cartesian_coordinates=np.vstack(coordinate_blocks),
        basis_vectors=unit_vectors * scale,
        replication=replication,
        chain_ids=tuple(chain_ids),
        residue_numbers=tuple(residue_numbers),
    )
