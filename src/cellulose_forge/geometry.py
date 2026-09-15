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
_LAYER_TOLERANCE = 1e-6


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
    chain_numbers: tuple[int, ...] = ()
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
        if self.chain_numbers and len(self.chain_numbers) != len(coordinates):
            msg = "chain numbers and Cartesian coordinates must have equal lengths"
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
        chain_numbers=tuple(0 if chain_id == "A" else 1 for chain_id in chain_ids),
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
    chain_numbers: list[int] = []
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
                    chain_numbers.append(chain_index)
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
        chain_numbers=tuple(chain_numbers),
    )


def build_chain_parallel_slab(  # noqa: PLR0913, PLR0915
    allomorph: AllomorphDefinition,
    face: MillerIndex,
    *,
    layers: int,
    chains_per_layer: int,
    degree_polymerization: int,
    vacuum: float,
) -> BulkCrystal:
    """Build a complete-residue slab for a chain-parallel ``(h k 0)`` face.

    This implementation follows the chain-layer construction in the cellulose
    wetting reference systems. Chains in a layer are spaced by ``(a-b)/2``;
    successive layers alternate origin and center chains and advance by
    ``(a+b)/2``. The chain direction follows ``c`` and the requested face
    normal is oriented along Cartesian z.
    """
    if allomorph.name != "I-beta":
        msg = "the initial slab builder supports cellulose I-beta only"
        raise ValueError(msg)
    if face.l != 0:
        msg = "I-beta slab faces must be parallel to the cellulose chain axis (l = 0)"
        raise ValueError(msg)
    if any(
        not isinstance(value, int) or isinstance(value, bool) or value < 1
        for value in (layers, chains_per_layer, degree_polymerization)
    ):
        msg = "layers, chains per layer, and degree of polymerization must be positive integers"
        raise ValueError(msg)
    if degree_polymerization % 2:
        msg = "the initial crystallographic slab requires an even degree of polymerization"
        raise ValueError(msg)
    if not np.isfinite(vacuum) or vacuum < 0:
        msg = "vacuum must be a finite, nonnegative angstrom value"
        raise ValueError(msg)

    if face.values not in {(1, 1, 0), (-1, -1, 0)}:
        return _build_exploratory_chain_parallel_slab(
            allomorph,
            face,
            layers=layers,
            chains_per_layer=chains_per_layer,
            degree_polymerization=degree_polymerization,
            vacuum=vacuum,
        )

    vectors = allomorph.unit_cell.vectors()
    in_plane_vector = (vectors[0] - vectors[1]) / 2.0
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
    chain_repeats = degree_polymerization // 2
    axes = np.vstack([surface_x, chain_axis, normal])
    coordinate_blocks: list[FloatArray] = []
    atom_labels: list[str] = []
    chain_ids: list[str] = []
    chain_numbers: list[int] = []
    residue_numbers: list[int] = []
    chain_count = layers * chains_per_layer
    residue_atoms = allomorph.asymmetric_unit.residue_atom_count
    origin_indices = np.r_[0:residue_atoms, 2 * residue_atoms : 3 * residue_atoms]
    center_indices = np.r_[residue_atoms : 2 * residue_atoms, 3 * residue_atoms : 4 * residue_atoms]
    origin_template = unit_cell.cartesian_coordinates[origin_indices]
    center_template = unit_cell.cartesian_coordinates[center_indices].copy()
    center_template[residue_atoms:] += vectors[0] + vectors[1]
    origin_labels = tuple(unit_cell.atom_labels[index] for index in origin_indices)
    center_labels = tuple(unit_cell.atom_labels[index] for index in center_indices)

    for layer in range(layers):
        is_center_layer = bool(layer % 2)
        template = center_template if is_center_layer else origin_template
        template_labels = center_labels if is_center_layer else origin_labels
        layer_translation = (layer // 2) * (vectors[0] + vectors[1])
        for surface_chain in range(chains_per_layer):
            chain_index = layer * chains_per_layer + surface_chain
            chain_id = _PDB_CHAIN_IDENTIFIERS[chain_index % len(_PDB_CHAIN_IDENTIFIERS)]
            surface_translation = surface_chain * in_plane_vector
            for chain_cell in range(chain_repeats):
                translation = layer_translation + surface_translation + chain_cell * vectors[2]
                coordinate_blocks.append(template + translation)
                atom_labels.extend(template_labels)
                chain_ids.extend([chain_id] * len(template))
                chain_numbers.extend([chain_index] * len(template))
                first_residue = 2 * chain_cell + 1
                residue_numbers.extend([first_residue] * residue_atoms)
                residue_numbers.extend([first_residue + 1] * residue_atoms)

    projected = np.vstack(coordinate_blocks) @ axes.T
    projected[:, 0] -= float(np.min(projected[:, 0]))
    projected[:, 1] -= float(np.min(projected[:, 1]))
    projected[:, 2] -= float(np.min(projected[:, 2])) - vacuum / 2.0
    slab_extent = float(np.ptp(projected[:, 2]))
    chain_length = float(np.linalg.norm(vectors[2])) * chain_repeats
    basis = np.diag([in_plane_length * chains_per_layer, chain_length, slab_extent + vacuum])
    return BulkCrystal(
        allomorph=allomorph.name,
        space_group="P1",
        atom_labels=tuple(atom_labels),
        cartesian_coordinates=projected,
        basis_vectors=basis,
        replication=(chains_per_layer, layers, chain_repeats),
        chain_ids=tuple(chain_ids),
        residue_numbers=tuple(residue_numbers),
        chain_numbers=tuple(chain_numbers),
        annotations=(
            f"shape=slab face=({face.h},{face.k},{face.l}) layers={layers}",
            f"chains_per_layer={chains_per_layer} total_chains={chain_count}",
            f"degree_polymerization={degree_polymerization}",
            f"d_hkl={face.spacing(allomorph.unit_cell):.6f} angstrom",
            f"vacuum={vacuum:.6f} angstrom",
            "surface normal aligned with Cartesian z",
        ),
    )


def _build_exploratory_chain_parallel_slab(  # noqa: PLR0913, PLR0915
    allomorph: AllomorphDefinition,
    face: MillerIndex,
    *,
    layers: int,
    chains_per_layer: int,
    degree_polymerization: int,
    vacuum: float,
) -> BulkCrystal:
    """Build an exploratory I-beta slab by selecting complete lattice chains.

    Candidate origin and center chains are generated on the two-dimensional
    ``a``/``b`` lattice, grouped into planes normal to ``(h k 0)``, and cut as
    complete chains. Unlike the dedicated (110) path, these faces have not yet
    been compared with an independent published atomic reference.
    """
    vectors = allomorph.unit_cell.vectors()
    normal = face.normal(allomorph.unit_cell)
    chain_axis = vectors[2] / np.linalg.norm(vectors[2])
    surface_x = np.cross(chain_axis, normal)
    surface_x /= np.linalg.norm(surface_x)
    axes = np.vstack([surface_x, chain_axis, normal])

    unit_cell = build_unit_cell(allomorph)
    residue_atoms = allomorph.asymmetric_unit.residue_atom_count
    origin_indices = np.r_[0:residue_atoms, 2 * residue_atoms : 3 * residue_atoms]
    center_indices = np.r_[residue_atoms : 2 * residue_atoms, 3 * residue_atoms : 4 * residue_atoms]
    origin_template = unit_cell.cartesian_coordinates[origin_indices]
    center_template = unit_cell.cartesian_coordinates[center_indices].copy()
    center_template[residue_atoms:] += vectors[0] + vectors[1]
    templates = (origin_template, center_template)
    template_labels = tuple(
        tuple(unit_cell.atom_labels[index] for index in indices)
        for indices in (origin_indices, center_indices)
    )

    search_radius = max(8, layers + chains_per_layer + abs(face.h) + abs(face.k))
    candidates: list[tuple[float, float, int, int, int]] = []
    for i in range(-search_radius, search_radius + 1):
        for j in range(-search_radius, search_radius + 1):
            translation = i * vectors[0] + j * vectors[1]
            for template_index, template in enumerate(templates):
                center = np.mean(template, axis=0) + translation
                candidates.append(
                    (
                        float(np.dot(center, normal)),
                        float(np.dot(center, surface_x)),
                        template_index,
                        i,
                        j,
                    )
                )

    candidates.sort(key=lambda item: (item[0], item[1], item[2], item[3], item[4]))
    rows: list[list[tuple[float, float, int, int, int]]] = []
    for candidate in candidates:
        if not rows or abs(candidate[0] - rows[-1][0][0]) > _LAYER_TOLERANCE:
            rows.append([candidate])
        else:
            rows[-1].append(candidate)
    eligible_rows = [row for row in rows if len(row) >= chains_per_layer]
    if len(eligible_rows) < layers:
        msg = f"could not construct {layers} complete molecular layers for face {face.values}"
        raise ValueError(msg)

    middle = len(eligible_rows) // 2
    start = max(0, min(middle - layers // 2, len(eligible_rows) - layers))
    selected_rows = eligible_rows[start : start + layers]
    selected: list[tuple[float, float, int, int, int]] = []
    for row in selected_rows:
        row.sort(key=lambda item: (item[1], item[2], item[3], item[4]))
        row_middle = len(row) // 2
        row_start = max(0, min(row_middle - chains_per_layer // 2, len(row) - chains_per_layer))
        selected.extend(row[row_start : row_start + chains_per_layer])

    chain_repeats = degree_polymerization // 2
    coordinate_blocks: list[FloatArray] = []
    atom_labels: list[str] = []
    chain_ids: list[str] = []
    chain_numbers: list[int] = []
    residue_numbers: list[int] = []
    for chain_index, (_, _, template_index, i, j) in enumerate(selected):
        chain_id = _PDB_CHAIN_IDENTIFIERS[chain_index % len(_PDB_CHAIN_IDENTIFIERS)]
        transverse_translation = i * vectors[0] + j * vectors[1]
        template = templates[template_index]
        labels = template_labels[template_index]
        for chain_cell in range(chain_repeats):
            translation = transverse_translation + chain_cell * vectors[2]
            coordinate_blocks.append(template + translation)
            atom_labels.extend(labels)
            chain_ids.extend([chain_id] * len(template))
            chain_numbers.extend([chain_index] * len(template))
            first_residue = 2 * chain_cell + 1
            residue_numbers.extend([first_residue] * residue_atoms)
            residue_numbers.extend([first_residue + 1] * residue_atoms)

    projected = np.vstack(coordinate_blocks) @ axes.T
    projected[:, 0] -= float(np.min(projected[:, 0]))
    projected[:, 1] -= float(np.min(projected[:, 1]))
    projected[:, 2] -= float(np.min(projected[:, 2])) - vacuum / 2.0
    extents = np.ptp(projected, axis=0)
    chain_length = float(np.linalg.norm(vectors[2])) * chain_repeats
    basis = np.diag([float(extents[0]), chain_length, float(extents[2]) + vacuum])
    return BulkCrystal(
        allomorph=allomorph.name,
        space_group="P1",
        atom_labels=tuple(atom_labels),
        cartesian_coordinates=projected,
        basis_vectors=basis,
        replication=(chains_per_layer, layers, chain_repeats),
        chain_ids=tuple(chain_ids),
        residue_numbers=tuple(residue_numbers),
        chain_numbers=tuple(chain_numbers),
        annotations=(
            f"shape=slab face=({face.h},{face.k},{face.l}) layers={layers}",
            f"chains_per_layer={chains_per_layer} total_chains={len(selected)}",
            f"degree_polymerization={degree_polymerization}",
            f"d_hkl={face.spacing(allomorph.unit_cell):.6f} angstrom",
            f"vacuum={vacuum:.6f} angstrom",
            "surface normal aligned with Cartesian z",
            "validation=exploratory; independent surface reference pending",
        ),
    )
