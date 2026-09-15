"""Crystallographic unit-cell and bulk-crystal geometry generation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from cellulose_forge.crystallography import (
    AllomorphDefinition,
    FloatArray,
    MillerIndex,
    fractional_to_cartesian,
)

_ARRAY_DIMENSIONS = 2
_COORDINATE_DIMENSIONS = 3
_ORTHOGONAL_TOLERANCE = 1e-10
_NONZERO_VECTOR_TOLERANCE = 1e-12
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
    annotations: tuple[str, ...] = ()

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


def build_chain_parallel_slab(  # noqa: PLR0913, PLR0915
    allomorph: AllomorphDefinition,
    face: MillerIndex,
    *,
    layers: int,
    surface_repeats: int,
    chain_repeats: int,
    vacuum: float,
) -> BulkCrystal:
    """Build a complete-residue slab for a chain-parallel ``(h k 0)`` face.

    This first surface implementation constructs the surface directly from
    translated complete unit cells. The first surface vector is ``k*a-h*b``;
    the second follows the cellulose-chain vector ``c``. The requested face
    normal is oriented along Cartesian z in the returned structure.
    """
    if face.l != 0:
        msg = "the initial slab builder supports chain-parallel faces with l = 0"
        raise ValueError(msg)
    if allomorph.name != "I-beta":
        msg = "the initial slab builder supports cellulose I-beta only"
        raise ValueError(msg)
    if any(
        not isinstance(value, int) or isinstance(value, bool) or value < 1
        for value in (layers, surface_repeats, chain_repeats)
    ):
        msg = "layers and repeat counts must be positive integers"
        raise ValueError(msg)
    if not np.isfinite(vacuum) or vacuum < 0:
        msg = "vacuum must be a finite, nonnegative angstrom value"
        raise ValueError(msg)

    vectors = allomorph.unit_cell.vectors()
    in_plane_vector = face.k * vectors[0] - face.h * vectors[1]
    in_plane_length = float(np.linalg.norm(in_plane_vector))
    if in_plane_length < _NONZERO_VECTOR_TOLERANCE:
        msg = "face does not define a nonzero transverse surface vector"
        raise ValueError(msg)

    normal = face.normal(allomorph.unit_cell)
    surface_x = in_plane_vector / in_plane_length
    chain_axis = vectors[2] / np.linalg.norm(vectors[2])
    if abs(float(np.dot(surface_x, chain_axis))) > _ORTHOGONAL_TOLERANCE:
        msg = (
            "initial slab basis requires the chain and transverse surface vectors to be orthogonal"
        )
        raise ValueError(msg)

    unit_cell = build_unit_cell(allomorph)
    axes = np.vstack([surface_x, chain_axis, normal])
    coordinate_blocks: list[FloatArray] = []
    atom_labels: list[str] = []
    chain_ids: list[str] = []
    residue_numbers: list[int] = []
    chains_per_cell = len(set(unit_cell.chain_ids))
    chain_count = layers * surface_repeats * chains_per_cell
    if chain_count > len(_PDB_CHAIN_IDENTIFIERS):
        msg = f"diagnostic PDB supports at most 62 slab chains; requested {chain_count}"
        raise ValueError(msg)

    unit_chain_order = tuple(dict.fromkeys(unit_cell.chain_ids))
    unit_chain_indices = {chain_id: index for index, chain_id in enumerate(unit_chain_order)}
    residues_per_chain = max(unit_cell.residue_numbers)
    for chain_cell in range(chain_repeats):
        for layer in range(layers):
            for surface_cell in range(surface_repeats):
                # k*a-h*b lies in the surface; b advances adjacent (hkl) planes
                translation = (
                    surface_cell * in_plane_vector + layer * vectors[1] + chain_cell * vectors[2]
                )
                coordinate_blocks.append(unit_cell.cartesian_coordinates + translation)
                atom_labels.extend(unit_cell.atom_labels)
                transverse_index = layer * surface_repeats + surface_cell
                for chain_id, residue_number in zip(
                    unit_cell.chain_ids, unit_cell.residue_numbers, strict=True
                ):
                    chain_index = transverse_index * chains_per_cell + unit_chain_indices[chain_id]
                    chain_ids.append(_PDB_CHAIN_IDENTIFIERS[chain_index])
                    residue_numbers.append(residue_number + chain_cell * residues_per_chain)

    projected = np.vstack(coordinate_blocks) @ axes.T
    projected[:, 0] -= float(np.min(projected[:, 0]))
    projected[:, 1] -= float(np.min(projected[:, 1]))
    projected[:, 2] -= float(np.min(projected[:, 2])) - vacuum / 2.0
    slab_extent = float(np.ptp(projected[:, 2]))
    chain_length = float(np.linalg.norm(vectors[2])) * chain_repeats
    basis = np.diag([in_plane_length * surface_repeats, chain_length, slab_extent + vacuum])
    return BulkCrystal(
        allomorph=allomorph.name,
        space_group="P1",
        atom_labels=tuple(atom_labels),
        cartesian_coordinates=projected,
        basis_vectors=basis,
        replication=(surface_repeats, layers, chain_repeats),
        chain_ids=tuple(chain_ids),
        residue_numbers=tuple(residue_numbers),
        annotations=(
            f"shape=slab face=({face.h},{face.k},{face.l}) layers={layers}",
            f"d_hkl={face.spacing(allomorph.unit_cell):.6f} angstrom",
            f"vacuum={vacuum:.6f} angstrom",
            "surface normal aligned with Cartesian z",
        ),
    )
