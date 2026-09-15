"""Crystallographic unit-cell and bulk-crystal geometry generation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from cellulose_forge.crystallography import (
    AllomorphDefinition,
    FloatArray,
    fractional_to_cartesian,
)

_I_BETA_ATOMS_PER_RESIDUE = 21
_I_BETA_ATOMS_PER_CELL = 84
_MINIMUM_TRANSVERSE_REPLICATION = 2
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


def _validate_i_beta_layout(allomorph: AllomorphDefinition) -> None:
    if allomorph.name != "I-beta":
        msg = "build_i_beta_parallelepiped requires the I-beta allomorph"
        raise ValueError(msg)
    if allomorph.asymmetric_unit.residue_atom_count != _I_BETA_ATOMS_PER_RESIDUE:
        msg = "the I-beta legacy replication requires 21 atoms per glucose residue"
        raise ValueError(msg)
    if len(_complete_cell_fractional(allomorph)) != _I_BETA_ATOMS_PER_CELL:
        msg = "the I-beta legacy replication requires 84 atoms per complete unit cell"
        raise ValueError(msg)


def _build_layer_block(
    allomorph: AllomorphDefinition, i: int, j: int
) -> tuple[FloatArray, tuple[str, ...]]:
    """Build one translated complete I-beta cell in legacy atom order."""
    block = _complete_cell_fractional(allomorph).copy()
    block[:, 0] += i
    block[:, 1] += j
    labels = allomorph.asymmetric_unit.atom_labels * len(allomorph.symmetry_operations)
    return block, labels


def _trim_edge_rings(
    block: FloatArray,
    labels: tuple[str, ...],
    cell_index: tuple[int, int],
    grid_shape: tuple[int, int],
) -> tuple[FloatArray, tuple[str, ...]]:
    """Delete boundary glucose rings using the legacy Cellulose-Builder rules."""
    i, j = cell_index
    xsize, ysize = grid_shape
    drop_center_image = (i == 0 and j != ysize - 1) or (j == 0 and i not in (0, xsize - 1))
    drop_center = (i == xsize - 1 and j != 0) or (j == ysize - 1 and i not in (0, xsize - 1))
    if (i == 0 and j == ysize - 1) or (j == 0 and i == xsize - 1):
        drop_center = True
        drop_center_image = True

    keep = np.ones(len(block), dtype=np.bool_)
    if drop_center:
        keep[21:42] = False
    if drop_center_image:
        keep[63:84] = False
    kept_labels = tuple(label for label, include in zip(labels, keep, strict=True) if include)
    return block[keep], kept_labels


def build_i_beta_parallelepiped(
    allomorph: AllomorphDefinition, xsize: int, ysize: int, zsize: int
) -> BulkCrystal:
    """Build a legacy-compatible I-beta parallelepiped with transverse PBC disabled."""
    _validate_i_beta_layout(allomorph)
    if xsize < _MINIMUM_TRANSVERSE_REPLICATION or ysize < _MINIMUM_TRANSVERSE_REPLICATION:
        msg = "xsize and ysize must each be at least 2 (legacy CLI constraint)"
        raise ValueError(msg)
    if zsize < 1:
        msg = "zsize must be at least 1"
        raise ValueError(msg)

    layer_blocks: list[FloatArray] = []
    layer_labels: list[str] = []
    for j in range(ysize):
        for i in range(xsize):
            block, labels = _build_layer_block(allomorph, i, j)
            trimmed, trimmed_labels = _trim_edge_rings(block, labels, (i, j), (xsize, ysize))
            layer_blocks.append(trimmed)
            layer_labels.extend(trimmed_labels)

    layer_cartesian = fractional_to_cartesian(np.vstack(layer_blocks), allomorph.unit_cell)
    cell_c = allomorph.unit_cell.vectors()[2]
    all_layers = [layer_cartesian + cell_c * k for k in range(zsize)]
    cartesian = np.vstack(all_layers)
    atom_labels = tuple(layer_labels) * zsize

    scale = np.asarray([xsize, ysize, zsize], dtype=np.float64)[:, np.newaxis]
    basis_vectors = allomorph.unit_cell.vectors() * scale
    return BulkCrystal(
        allomorph=allomorph.name,
        space_group=allomorph.space_group,
        atom_labels=atom_labels,
        cartesian_coordinates=cartesian,
        basis_vectors=basis_vectors,
        replication=(xsize, ysize, zsize),
    )
