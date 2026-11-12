"""Crystallographic unit-cell and bulk-crystal geometry generation."""

from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import pairwise

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
_MAX_SLAB_DISCRETE_COUNT = 10_000
_PROPER_ROTATION_TOLERANCE = 1e-10
_CHIRALITY_VOLUME_TOLERANCE = 1e-6
_I_BETA_TERMINATION_BY_FACE = {
    (1, 0, 0): 0,
    (-1, 0, 0): 0,
    (2, 0, 0): 1,
    (-2, 0, 0): 1,
    (0, 1, 0): 0,
    (0, -1, 0): 0,
    (0, 2, 0): 1,
    (0, -2, 0): 1,
}
_I_BETA_CHAIN_TYPE_NAMES = ("origin", "center")


def _proper_slab_axes(
    surface_x: FloatArray,
    chain_axis: FloatArray,
    normal: FloatArray,
) -> FloatArray:
    """Return a right-handed slab basis with chains along -y and normal along -z."""
    axes = np.vstack([surface_x, -chain_axis, -normal])
    if np.linalg.det(axes) < 0:
        axes[0] *= -1
    determinant = float(np.linalg.det(axes))
    if not np.isclose(determinant, 1.0, atol=_PROPER_ROTATION_TOLERANCE):
        msg = f"slab orientation is not a proper rotation (determinant={determinant:.12g})"
        raise ValueError(msg)
    return axes


def validate_cellulose_chirality(crystal: BulkCrystal) -> dict[str, object]:
    """Validate the C1 configuration of every complete glucose residue.

    The signed volume uses the ordered C1 substituents O5, C2, and H1.  The
    crystallographic D-cellulose definitions supplied with CelluloseForge have
    a negative sign under this convention.  Checking the molecular geometry,
    rather than only a transformation determinant, also detects an already
    mirrored input structure.
    """
    labels = np.asarray(crystal.atom_labels)
    chains = np.asarray(crystal.chain_numbers)
    residues = np.asarray(crystal.residue_numbers)
    signed_volumes: list[float] = []
    for chain, residue in sorted(set(zip(chains, residues, strict=True))):
        mask = (chains == chain) & (residues == residue)
        lookup = dict(zip(labels[mask], crystal.cartesian_coordinates[mask], strict=True))
        required = {"C1", "O5", "C2", "H1"}
        if not required.issubset(lookup):
            msg = (
                "cannot validate cellulose chirality for "
                f"chain {chain}, residue {residue}: missing "
                f"{', '.join(sorted(required - lookup.keys()))}"
            )
            raise ValueError(msg)
        center = lookup["C1"]
        signed_volumes.append(
            float(
                np.linalg.det(
                    np.vstack(
                        [
                            lookup["O5"] - center,
                            lookup["C2"] - center,
                            lookup["H1"] - center,
                        ]
                    )
                )
            )
        )

    if not signed_volumes:
        raise ValueError("cannot validate cellulose chirality: no residues found")
    volumes = np.asarray(signed_volumes)
    if np.any(np.abs(volumes) <= _CHIRALITY_VOLUME_TOLERANCE):
        raise ValueError("cellulose chirality validation found a degenerate C1 centre")
    if np.any(volumes > 0.0):
        raise ValueError("cellulose chirality validation found mirrored C1 geometry")
    return {
        "status": "passed",
        "configuration": "D-cellulose",
        "stereocentre": "C1",
        "substituent_order": ["O5", "C2", "H1"],
        "expected_signed_volume_sign": "negative",
        "residue_count": len(signed_volumes),
        "minimum_absolute_signed_volume_angstrom3": float(np.min(np.abs(volumes))),
    }


# Signed-volume signs of the beta-D-pyranose ring stereocentres, using the same
# convention as the cellulose C1 check. N-acetyl-beta-D-glucosamine shares the
# configuration of beta-D-glucose at every ring carbon.
_PYRANOSE_STEREOCENTRES: tuple[tuple[str, tuple[str, str, str], int], ...] = (
    ("C1", ("O5", "C2", "H1"), -1),
    ("C2", ("C1", "C3", "H2"), 1),
    ("C3", ("C2", "C4", "H3"), -1),
    ("C4", ("C3", "C5", "H4"), 1),
    ("C5", ("C4", "O5", "H5"), -1),
)

# beta-D-mannose is the C2 epimer of beta-D-glucose. The other four ring
# stereocentres therefore retain the signs used above.
_MANNAN_STEREOCENTRES: tuple[tuple[str, tuple[str, str, str], int], ...] = (
    ("C1", ("O5", "C2", "H1"), -1),
    ("C2", ("C1", "C3", "H2"), -1),
    ("C3", ("C2", "C4", "H3"), -1),
    ("C4", ("C3", "C5", "H4"), 1),
    ("C5", ("C4", "O5", "H5"), -1),
)


def validate_pyranose_chirality(
    crystal: BulkCrystal,
    configuration: str,
    *,
    stereocentres: tuple[tuple[str, tuple[str, str, str], int], ...] = _PYRANOSE_STEREOCENTRES,
) -> dict[str, object]:
    """Validate all five ring stereocentres of every beta-D-pyranose residue."""
    labels = np.asarray(crystal.atom_labels)
    chains = np.asarray(crystal.chain_numbers)
    residues = np.asarray(crystal.residue_numbers)
    minimum: dict[str, float] = {}
    residue_count = 0
    for chain, residue in sorted(set(zip(chains, residues, strict=True))):
        mask = (chains == chain) & (residues == residue)
        lookup = dict(zip(labels[mask], crystal.cartesian_coordinates[mask], strict=True))
        residue_count += 1
        for centre, substituents, expected in stereocentres:
            required = {centre, *substituents}
            if not required.issubset(lookup):
                msg = (
                    f"cannot validate {configuration} chirality for chain {chain}, "
                    f"residue {residue}: missing "
                    f"{', '.join(sorted(required - lookup.keys()))}"
                )
                raise ValueError(msg)
            volume = float(
                np.linalg.det(np.vstack([lookup[name] - lookup[centre] for name in substituents]))
            )
            if abs(volume) <= _CHIRALITY_VOLUME_TOLERANCE:
                raise ValueError(f"{configuration} chirality: degenerate {centre} centre")
            if np.sign(volume) != expected:
                raise ValueError(
                    f"{configuration} chirality: inverted {centre} in chain {chain}, "
                    f"residue {residue}"
                )
            minimum[centre] = min(minimum.get(centre, abs(volume)), abs(volume))
    if not residue_count:
        raise ValueError(f"cannot validate {configuration} chirality: no residues found")
    return {
        "status": "passed",
        "configuration": configuration,
        "stereocentres": [centre for centre, _, _ in stereocentres],
        "substituent_orders": {
            centre: list(substituents) for centre, substituents, _ in stereocentres
        },
        "expected_signed_volume_signs": {
            centre: "positive" if expected > 0 else "negative"
            for centre, _, expected in stereocentres
        },
        "residue_count": residue_count,
        "minimum_absolute_signed_volume_angstrom3": minimum,
    }


def validate_mannose_chirality(crystal: BulkCrystal) -> dict[str, object]:
    """Validate all ring stereocentres as beta-D-mannose."""
    return validate_pyranose_chirality(
        crystal,
        "beta-D-mannose",
        stereocentres=_MANNAN_STEREOCENTRES,
    )


_FIBRIL_FRAGMENT_IDS = {
    "I-alpha": (
        2,
        3,
        4,
        5,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        36,
        37,
        38,
        39,
    ),
    "I-beta": (
        18,
        27,
        36,
        10,
        19,
        28,
        37,
        11,
        20,
        29,
        38,
        47,
        3,
        12,
        21,
        30,
        39,
        48,
        4,
        13,
        22,
        31,
        40,
        49,
        5,
        14,
        23,
        32,
        41,
        15,
        24,
        33,
        42,
        16,
        25,
        34,
    ),
    "II": (
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        41,
        42,
        43,
        44,
        45,
        46,
        47,
    ),
}

# Ordered bottom-to-top row populations in the published 36-chain
# cross-sections after applying the documented counterclockwise rotation.
_FIBRIL_LAYER_SIZES = {
    "I-alpha": (4, 5, 6, 6, 6, 5, 4),
    "I-beta": (3, 4, 5, 6, 6, 5, 4, 3),
    "II": (7, 11, 11, 7),
}

I_BETA_FIBRIL_PROFILES: dict[str, tuple[int, ...]] = {
    "18-chain": (2, 3, 4, 4, 3, 2),
    "24-chain": (3, 4, 5, 5, 4, 3),
}


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
    layer_numbers: tuple[int, ...] = ()
    annotations: tuple[str, ...] = ()
    residue_name: str = "BGC"

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
        if self.layer_numbers and len(self.layer_numbers) != len(coordinates):
            msg = "layer numbers and Cartesian coordinates must have equal lengths"
            raise ValueError(msg)
        if self.chain_numbers and min(self.chain_numbers) < 1:
            msg = "chain numbers must be one-based positive integers"
            raise ValueError(msg)
        if self.layer_numbers and min(self.layer_numbers) < 1:
            msg = "layer numbers must be one-based positive integers"
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
    residues_per_chain = {
        chain_name: asymmetric.chain_labels.count(chain_name) for chain_name in chain_names
    }
    for operation_index in range(operation_count):
        chain_residue_counts = dict.fromkeys(chain_names, 0)
        for residue_index in range(residues_per_asymmetric_unit):
            chain_name = asymmetric.chain_labels[residue_index]
            chain_residue_counts[chain_name] += 1
            chain_ids.extend([chain_letters[chain_name]] * asymmetric.residue_atom_count)
            residue_number = (
                operation_index * residues_per_chain[chain_name] + chain_residue_counts[chain_name]
            )
            residue_numbers.extend([residue_number] * asymmetric.residue_atom_count)

    return BulkCrystal(
        allomorph=allomorph.name,
        space_group=allomorph.space_group,
        atom_labels=atom_labels,
        cartesian_coordinates=fractional_to_cartesian(fractional, allomorph.unit_cell),
        basis_vectors=allomorph.unit_cell.vectors(),
        replication=(1, 1, 1),
        chain_ids=tuple(chain_ids),
        residue_numbers=tuple(residue_numbers),
        chain_numbers=tuple(1 if chain_id == "A" else 2 for chain_id in chain_ids),
        annotations=_allomorph_chain_annotations(allomorph),
        residue_name=allomorph.residue_name,
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
    cell_counts = replication
    chains_per_cell = len(dict.fromkeys(allomorph.asymmetric_unit.chain_labels))
    transverse_axes = tuple(axis for axis in range(3) if axis != allomorph.chain_axis)
    transverse_chain_count = (
        cell_counts[transverse_axes[0]] * cell_counts[transverse_axes[1]] * chains_per_cell
    )
    if transverse_chain_count > len(_PDB_CHAIN_IDENTIFIERS):
        msg = (
            "diagnostic PDB chain identifiers support at most 62 transverse "
            f"chains; requested {transverse_chain_count}"
        )
        raise ValueError(msg)

    unit_cell = build_unit_cell(allomorph)
    unit_vectors = allomorph.unit_cell.vectors()
    unit_coordinates = _unwrap_unit_cell_chains(
        unit_cell,
        unit_vectors,
        chain_axis=allomorph.chain_axis,
    )
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
                coordinate_blocks.append(unit_coordinates + translation)
                atom_labels.extend(unit_cell.atom_labels)

                cell_indices = (i, j, k)
                transverse_cell_index = (
                    cell_indices[transverse_axes[0]]
                    + cell_counts[transverse_axes[0]] * cell_indices[transverse_axes[1]]
                )
                for chain_id, residue_number in zip(
                    unit_cell.chain_ids, unit_cell.residue_numbers, strict=True
                ):
                    chain_index = (
                        transverse_cell_index * chains_per_cell + unit_chain_indices[chain_id]
                    )
                    chain_ids.append(_PDB_CHAIN_IDENTIFIERS[chain_index])
                    chain_numbers.append(chain_index + 1)
                    residue_numbers.append(
                        residue_number
                        + cell_indices[allomorph.chain_axis] * residues_per_chain_per_cell
                    )

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
        annotations=_allomorph_chain_annotations(allomorph),
        residue_name=allomorph.residue_name,
    )


def _allomorph_chain_annotations(allomorph: AllomorphDefinition) -> tuple[str, ...]:
    """Return machine-readable chain-axis and direction annotations."""
    axis_name = ("a", "b", "c")[allomorph.chain_axis]
    directions = dict(
        zip(
            allomorph.asymmetric_unit.chain_labels,
            allomorph.asymmetric_unit.chain_directions,
            strict=True,
        )
    )
    direction_text = ",".join(
        f"{label}:{'+' if direction > 0 else '-'}{axis_name}"
        for label, direction in directions.items()
    )
    return f"chain_axis={axis_name}", f"chain_directions={direction_text}"


def _unwrap_unit_cell_chains(
    unit_cell: BulkCrystal,
    unit_vectors: FloatArray,
    *,
    chain_axis: int,
) -> FloatArray:
    """Make sequential residues whole across transverse unit-cell boundaries.

    Crystallographic symmetry can place consecutive residues of one cellulose
    chain in neighbouring transverse images.  Molecular topology generation
    requires one consistent image of each chain, so choose integer ``a``/``b``
    translations that minimize the lateral displacement between consecutive
    C1 atoms.  Replication along ``c`` is applied only after this operation.
    """
    coordinates = unit_cell.cartesian_coordinates.copy()
    labels = np.asarray(unit_cell.atom_labels)
    chains = np.asarray(unit_cell.chain_numbers)
    residues = np.asarray(unit_cell.residue_numbers)
    chain_vector = unit_vectors[chain_axis] / np.linalg.norm(unit_vectors[chain_axis])
    transverse_axes = tuple(axis for axis in range(3) if axis != chain_axis)
    transverse_translations = tuple(
        i * unit_vectors[transverse_axes[0]] + j * unit_vectors[transverse_axes[1]]
        for i in range(-2, 3)
        for j in range(-2, 3)
    )

    for chain_number in dict.fromkeys(unit_cell.chain_numbers):
        residue_order = tuple(dict.fromkeys(residues[chains == chain_number]))
        previous_anchor: FloatArray | None = None
        for residue_number in residue_order:
            residue_mask = (chains == chain_number) & (residues == residue_number)
            anchor_indices = np.flatnonzero(residue_mask & (labels == "C1"))
            if len(anchor_indices) != 1:
                msg = "each crystallographic cellulose residue must contain exactly one C1 atom"
                raise ValueError(msg)
            anchor = coordinates[anchor_indices[0]]
            if previous_anchor is not None:
                translation = min(
                    transverse_translations,
                    key=lambda candidate: np.linalg.norm(
                        (anchor + candidate - previous_anchor)
                        - np.dot(anchor + candidate - previous_anchor, chain_vector) * chain_vector
                    ),
                )
                coordinates[residue_mask] += translation
                anchor = anchor + translation
            previous_anchor = anchor

    coordinates.setflags(write=False)
    return coordinates


def build_36_chain_fibril(  # noqa: PLR0912, PLR0913, PLR0915
    allomorph: AllomorphDefinition,
    *,
    degree_polymerization: int,
    vacuum: float,
    padding_x: float | None = None,
    padding_y: float | None = None,
    padding_z: float | None = None,
    chain_boundary: str = "finite",
) -> BulkCrystal:
    """Build the literature-defined 36-chain elementary fibril.

    The transverse chain selections reproduce the explicit fragment lists in
    Cellulose-Builder. They are discrete, allomorph-specific arrangements and
    are deliberately not approximated as regular hexagons.
    """
    if allomorph.name not in _FIBRIL_FRAGMENT_IDS:
        msg = (
            f"the 36-chain fibril is not defined for cellulose {allomorph.name}; "
            "the reference Cellulose-Builder supports I-alpha, I-beta, and II"
        )
        raise ValueError(msg)
    if (
        not isinstance(degree_polymerization, int)
        or isinstance(degree_polymerization, bool)
        or degree_polymerization < 1
    ):
        msg = "degree of polymerization must be a positive integer"
        raise ValueError(msg)
    if degree_polymerization % 2:
        msg = "the crystallographic 36-chain fibril requires an even degree of polymerization"
        raise ValueError(msg)
    if chain_boundary not in {"finite", "periodic"}:
        msg = "chain boundary must be 'finite' or 'periodic'"
        raise ValueError(msg)
    if not np.isfinite(vacuum) or vacuum < 0:
        msg = "vacuum must be a finite, nonnegative angstrom value"
        raise ValueError(msg)
    for name, value in (
        ("padding_x", padding_x),
        ("padding_y", padding_y),
        ("padding_z", padding_z),
    ):
        if value is not None and (not np.isfinite(value) or value < 0):
            msg = f"{name} must be a finite, nonnegative angstrom value"
            raise ValueError(msg)
    if chain_boundary == "periodic" and padding_z not in {None, 0}:
        msg = "padding_z must be zero for a fibril with periodic chain boundaries"
        raise ValueError(msg)
    transverse_padding_x = vacuum if padding_x is None else padding_x
    transverse_padding_z = vacuum if padding_y is None else padding_y
    longitudinal_padding = 0.0 if padding_z is None else padding_z

    cellobiose_units = degree_polymerization // 2
    vectors = allomorph.unit_cell.vectors()
    unit_cell = build_unit_cell(allomorph)
    unit_coordinates = _unwrap_unit_cell_chains(
        unit_cell,
        vectors,
        chain_axis=allomorph.chain_axis,
    )


    unit_chain_order = tuple(dict.fromkeys(unit_cell.chain_numbers))
    templates = {
        chain_number: (
            unit_coordinates[np.asarray(unit_cell.chain_numbers) == chain_number],
            tuple(
                label
                for label, number in zip(
                    unit_cell.atom_labels, unit_cell.chain_numbers, strict=True
                )
                if number == chain_number
            ),
            tuple(
                residue
                for residue, number in zip(
                    unit_cell.residue_numbers, unit_cell.chain_numbers, strict=True
                )
                if number == chain_number
            ),
        )
        for chain_number in unit_chain_order
    }
    direction_by_chain = {
        index + 1: direction
        for index, direction in enumerate(
            dict(
                zip(
                    allomorph.asymmetric_unit.chain_labels,
                    allomorph.asymmetric_unit.chain_directions,
                    strict=True,
                )
            ).values()
        )
    }

    candidates: list[tuple[int, FloatArray]] = []
    if allomorph.name == "I-alpha":
        for k in range(7):
            candidates.extend((1, j * vectors[1] + k * vectors[2]) for j in range(6))
    else:
        nx, ny = (5, 7) if allomorph.name == "I-beta" else (7, 5)
        for j in range(ny):
            for i in range(nx):
                translation = i * vectors[0] + j * vectors[1]
                candidates.append((1, translation))
                if i < nx - 1 and j < ny - 1:
                    candidates.append((2, translation))

    fragment_ids = _FIBRIL_FRAGMENT_IDS[allomorph.name]
    selected = [candidates[fragment_id] for fragment_id in fragment_ids]
    coordinate_blocks: list[FloatArray] = []
    atom_labels: list[str] = []
    chain_ids: list[str] = []
    chain_numbers: list[int] = []
    residue_numbers: list[int] = []
    for chain_index, (template_number, transverse_translation) in enumerate(selected):
        template, labels, template_residues = templates[template_number]
        direction = direction_by_chain[template_number]
        for repeat in range(cellobiose_units):
            translation = (
                transverse_translation + repeat * direction * vectors[allomorph.chain_axis]
            )
            coordinate_blocks.append(template + translation)
            atom_labels.extend(labels)
            chain_number = chain_index + 1
            chain_ids.extend([_PDB_CHAIN_IDENTIFIERS[chain_index]] * len(template))
            chain_numbers.extend([chain_number] * len(template))
            residue_numbers.extend(residue + 2 * repeat for residue in template_residues)

    chain_vector = vectors[allomorph.chain_axis]
    chain_axis = chain_vector / np.linalg.norm(chain_vector)
    transverse_axis = next(axis for axis in range(3) if axis != allomorph.chain_axis)
    surface_x = vectors[transverse_axis] - np.dot(vectors[transverse_axis], chain_axis) * chain_axis
    surface_x /= np.linalg.norm(surface_x)
    surface_y = np.cross(chain_axis, surface_x)
    # Rotate the transverse cross-section 90 degrees counterclockwise relative
    # to the crystallographic projection used by the legacy builder.
    axes = np.vstack([-surface_y, -chain_axis, surface_x])
    projected = np.vstack(coordinate_blocks) @ axes.T
    projected[:, 0] -= float(np.min(projected[:, 0])) - transverse_padding_x / 2.0
    projected[:, 1] -= float(np.min(projected[:, 1])) - longitudinal_padding / 2.0
    projected[:, 2] -= float(np.min(projected[:, 2])) - transverse_padding_z / 2.0
    extents = np.ptp(projected, axis=0)
    basis = np.diag(
        [
            float(extents[0]) + transverse_padding_x,
            float(np.linalg.norm(chain_vector)) * cellobiose_units + longitudinal_padding,
            float(extents[2]) + transverse_padding_z,
        ]
    )
    replication_values = [1, 1, 1]
    replication_values[allomorph.chain_axis] = cellobiose_units
    replication = (
        replication_values[0],
        replication_values[1],
        replication_values[2],
    )
    original_layer_by_chain: dict[int, int] = {}
    first_chain = 1
    for layer_number, layer_size in enumerate(_FIBRIL_LAYER_SIZES[allomorph.name], start=1):
        for chain_number in range(first_chain, first_chain + layer_size):
            original_layer_by_chain[chain_number] = layer_number
        first_chain += layer_size
    original_chain_numbers = np.asarray(chain_numbers)
    new_chain_by_original: dict[int, int] = {}
    next_chain_number = 1
    for layer_number in range(1, len(_FIBRIL_LAYER_SIZES[allomorph.name]) + 1):
        layer_chains = [
            chain_number
            for chain_number, assigned_layer in original_layer_by_chain.items()
            if assigned_layer == layer_number
        ]
        layer_chains.sort(
            key=lambda chain_number: float(
                projected[original_chain_numbers == chain_number, 0].mean()
            )
        )
        for original_chain_number in layer_chains:
            new_chain_by_original[original_chain_number] = next_chain_number
            next_chain_number += 1
    renumbered_chains = tuple(new_chain_by_original[value] for value in chain_numbers)
    renumbered_chain_ids = tuple(
        _PDB_CHAIN_IDENTIFIERS[chain_number - 1] for chain_number in renumbered_chains
    )
    layer_numbers = tuple(original_layer_by_chain[chain_number] for chain_number in chain_numbers)
    return BulkCrystal(
        allomorph=allomorph.name,
        space_group="P1",
        atom_labels=tuple(atom_labels),
        cartesian_coordinates=projected,
        basis_vectors=basis,
        replication=replication,
        chain_ids=renumbered_chain_ids,
        residue_numbers=tuple(residue_numbers),
        chain_numbers=renumbered_chains,
        layer_numbers=layer_numbers,
        annotations=(
            "shape=fibril fibril_model=36-chain total_chains=36",
            f"degree_polymerization={degree_polymerization} cellobiose_units={cellobiose_units}",
            f"chain_boundary={chain_boundary}",
            f"vacuum={vacuum:.6f} angstrom transverse_padding_default",
            f"padding_x={transverse_padding_x:.6f} angstrom total",
            f"padding_y={transverse_padding_z:.6f} angstrom total transverse-z padding",
            f"padding_z={longitudinal_padding:.6f} angstrom total longitudinal padding",
            "chains aligned with Cartesian -y",
            "fibril cross-section rotated 90 degrees counterclockwise",
            "fibril restraint layers numbered from low Cartesian z to high Cartesian z",
            "fibril chains numbered from low Cartesian x to high Cartesian x within each layer",
            "legacy_fragment_ids=" + ",".join(str(value) for value in fragment_ids),
            "source=Gomes_and_Skaf_2012_Cellulose-Builder doi:10.1002/jcc.22959",
        ),
        residue_name=allomorph.residue_name,
    )


def build_i_beta_profile_fibril(  # noqa: PLR0913, PLR0915
    allomorph: AllomorphDefinition,
    *,
    row_counts: tuple[int, ...],
    degree_polymerization: int,
    vacuum: float,
    padding_x: float | None = None,
    padding_y: float | None = None,
    padding_z: float | None = None,
    chain_boundary: str = "finite",
    model_name: str = "custom",
) -> BulkCrystal:
    """Build an I-beta fibril by selecting centred rows from the 36-chain lattice.

    Rows are supplied from low to high Cartesian z, matching CelluloseForge's
    layer-numbering convention.  The current profile builder deliberately uses
    the validated I-beta 36-chain lattice as its envelope; extending the same
    syntax to another allomorph requires an independently defined lattice-site
    mapping rather than reusing I-beta row semantics.
    """
    if allomorph.name != "I-beta":
        msg = "custom and 18-/24-chain fibril profiles currently support cellulose I-beta only"
        raise ValueError(msg)
    if not row_counts or any(
        not isinstance(value, int) or isinstance(value, bool) or value < 1
        for value in row_counts
    ):
        msg = "fibril row counts must be a nonempty sequence of positive integers"
        raise ValueError(msg)
    source_profile = _FIBRIL_LAYER_SIZES["I-beta"]
    if len(row_counts) > len(source_profile):
        msg = f"I-beta fibril profiles may contain at most {len(source_profile)} rows"
        raise ValueError(msg)
    if sum(row_counts) > len(_PDB_CHAIN_IDENTIFIERS):
        msg = f"fibril profiles may contain at most {len(_PDB_CHAIN_IDENTIFIERS)} chains"
        raise ValueError(msg)

    source_start = (len(source_profile) - len(row_counts)) // 2
    source_layers = tuple(
        range(source_start + 1, source_start + len(row_counts) + 1)
    )
    available_counts = tuple(source_profile[layer - 1] for layer in source_layers)
    oversized = [
        (index + 1, requested, available)
        for index, (requested, available) in enumerate(
            zip(row_counts, available_counts, strict=True)
        )
        if requested > available
    ]
    if oversized:
        details = ", ".join(
            f"row {row}: requested {requested}, maximum {available}"
            for row, requested, available in oversized
        )
        msg = f"profile does not fit the validated I-beta lattice envelope ({details})"
        raise ValueError(msg)

    full = build_36_chain_fibril(
        allomorph,
        degree_polymerization=degree_polymerization,
        vacuum=vacuum,
        padding_x=padding_x,
        padding_y=padding_y,
        padding_z=padding_z,
        chain_boundary=chain_boundary,
    )
    coordinates = np.asarray(full.cartesian_coordinates)
    chain_numbers = np.asarray(full.chain_numbers)
    layer_numbers = np.asarray(full.layer_numbers)

    selected_chains: list[int] = []
    for profile_row, (source_layer, requested) in enumerate(
        zip(source_layers, row_counts, strict=True)
    ):
        layer_chains = sorted(
            set(chain_numbers[layer_numbers == source_layer]),
            key=lambda chain: float(coordinates[chain_numbers == chain, 0].mean()),
        )
        excess = len(layer_chains) - requested
        selection_start = excess // 2
        if excess % 2 and profile_row % 2:
            selection_start += 1
        selected_chains.extend(
            layer_chains[selection_start : selection_start + requested]
        )

    selected_mask = np.isin(chain_numbers, selected_chains)
    selected_coordinates = coordinates[selected_mask].copy()
    selected_old_chains = chain_numbers[selected_mask]
    selected_old_layers = layer_numbers[selected_mask]

    transverse_padding_x = vacuum if padding_x is None else padding_x
    transverse_padding_z = vacuum if padding_y is None else padding_y
    selected_coordinates[:, 0] -= (
        float(np.min(selected_coordinates[:, 0])) - transverse_padding_x / 2.0
    )
    selected_coordinates[:, 2] -= (
        float(np.min(selected_coordinates[:, 2])) - transverse_padding_z / 2.0
    )
    basis = np.asarray(full.basis_vectors).copy()
    basis[0, 0] = float(np.ptp(selected_coordinates[:, 0])) + transverse_padding_x
    basis[2, 2] = float(np.ptp(selected_coordinates[:, 2])) + transverse_padding_z

    chain_mapping: dict[int, int] = {}
    layer_mapping: dict[int, int] = {
        source_layer: layer for layer, source_layer in enumerate(source_layers, start=1)
    }
    next_chain = 1
    for source_layer in source_layers:
        layer_selected = [
            chain
            for chain in selected_chains
            if int(selected_old_layers[selected_old_chains == chain][0]) == source_layer
        ]
        layer_selected.sort(
            key=lambda chain: float(
                selected_coordinates[selected_old_chains == chain, 0].mean()
            )
        )
        for old_chain in layer_selected:
            chain_mapping[old_chain] = next_chain
            next_chain += 1

    new_chains = tuple(chain_mapping[int(chain)] for chain in selected_old_chains)
    new_layers = tuple(layer_mapping[int(layer)] for layer in selected_old_layers)
    model_source = {
        "18-chain": (
            "Cellulose_Modeler_v0.8.0 profile=2,3,4,4,3,2; "
            "Kubicki_et_al_2018 doi:10.1038/s41598-018-32211-w"
        ),
        "24-chain": (
            "profile=3,4,5,5,4,3; Yang_and_Kubicki_2020 "
            "doi:10.1007/s10570-020-02970-9"
        ),
    }.get(model_name, "user_defined_I-beta_profile")
    selected_labels = tuple(
        label for label, keep in zip(full.atom_labels, selected_mask, strict=True) if keep
    )
    selected_residues = tuple(
        residue
        for residue, keep in zip(full.residue_numbers, selected_mask, strict=True)
        if keep
    )
    return BulkCrystal(
        allomorph=full.allomorph,
        space_group=full.space_group,
        residue_name=full.residue_name,
        atom_labels=selected_labels,
        cartesian_coordinates=selected_coordinates,
        basis_vectors=basis,
        replication=full.replication,
        chain_ids=tuple(_PDB_CHAIN_IDENTIFIERS[chain - 1] for chain in new_chains),
        residue_numbers=selected_residues,
        chain_numbers=new_chains,
        layer_numbers=new_layers,
        annotations=(
            f"shape=fibril fibril_model={model_name} total_chains={sum(row_counts)}",
            "fibril_profile_bottom_to_top=" + ",".join(str(value) for value in row_counts),
            f"degree_polymerization={degree_polymerization}",
            f"chain_boundary={chain_boundary}",
            "chains aligned with Cartesian -y",
            "fibril restraint layers numbered from low Cartesian z to high Cartesian z",
            "fibril chains numbered from low Cartesian x to high Cartesian x within each layer",
            f"source={model_source}",
        ),
    )


def build_chain_parallel_slab(  # noqa: PLR0912, PLR0913, PLR0915
    allomorph: AllomorphDefinition,
    face: MillerIndex,
    *,
    layers: int,
    chains_per_layer: int,
    degree_polymerization: int,
    vacuum: float,
    padding_z: float = 0.0,
    chain_boundary: str = "finite",
) -> BulkCrystal:
    """Build a complete-residue slab for a chain-parallel ``(h k 0)`` face.

    This implementation follows the chain-layer construction in the cellulose
    wetting reference systems. Chains in a layer are spaced by ``(a-b)/2``
    and alternate between origin and center sites; successive layers advance
    by ``(a+b)/2``. The chain direction follows ``c`` and the requested face
    normal is oriented along Cartesian -z and the cellulose chains along -y.
    The resulting surface lies in the x-y plane for semi-isotropic pressure
    coupling. A proper rotation is used so molecular chirality is preserved.
    """
    if face.values[allomorph.chain_axis] != 0:
        axis_name = ("h", "k", "l")[allomorph.chain_axis]
        msg = (
            "slab faces must be parallel to the cellulose chain axis "
            f"({axis_name} = 0 for cellulose {allomorph.name})"
        )
        raise ValueError(msg)
    if chain_boundary not in {"finite", "periodic"}:
        msg = "chain boundary must be 'finite' or 'periodic'"
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
    if not np.isfinite(padding_z) or padding_z < 0:
        msg = "padding_z must be a finite, nonnegative angstrom value"
        raise ValueError(msg)
    if chain_boundary == "periodic" and padding_z != 0:
        msg = "padding_z must be zero for a slab with periodic chain boundaries"
        raise ValueError(msg)

    if allomorph.name != "I-beta" or face.values not in {(1, 1, 0), (-1, -1, 0)}:
        return _build_exploratory_chain_parallel_slab(
            allomorph,
            face,
            layers=layers,
            chains_per_layer=chains_per_layer,
            degree_polymerization=degree_polymerization,
            vacuum=vacuum,
            padding_z=padding_z,
            chain_boundary=chain_boundary,
        )
    if chains_per_layer % 2:
        msg = (
            "cellulose I-beta (110) slabs require an even number of chains per "
            "layer so the transverse box is a complete lattice repeat"
        )
        raise ValueError(msg)

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
    axes = _proper_slab_axes(surface_x, chain_axis, normal)
    surface_order_reversed = float(np.dot(axes[0], surface_x)) < 0
    coordinate_blocks: list[FloatArray] = []
    atom_labels: list[str] = []
    chain_ids: list[str] = []
    chain_numbers: list[int] = []
    residue_numbers: list[int] = []
    layer_numbers: list[int] = []
    chain_count = layers * chains_per_layer
    residue_atoms = allomorph.asymmetric_unit.residue_atom_count
    origin_indices = np.r_[0:residue_atoms, 2 * residue_atoms : 3 * residue_atoms]
    center_indices = np.r_[residue_atoms : 2 * residue_atoms, 3 * residue_atoms : 4 * residue_atoms]
    origin_template = unit_cell.cartesian_coordinates[origin_indices]
    center_template = unit_cell.cartesian_coordinates[center_indices].copy()
    center_template[residue_atoms:] += vectors[0] + vectors[1]
    origin_labels = tuple(unit_cell.atom_labels[index] for index in origin_indices)
    center_labels = tuple(unit_cell.atom_labels[index] for index in center_indices)
    layer_step = (vectors[0] + vectors[1]) / 2.0

    for output_layer, crystal_layer in enumerate(reversed(range(layers))):
        for surface_chain in range(chains_per_layer):
            chain_index = output_layer * chains_per_layer + surface_chain
            chain_number = chain_index + 1
            chain_id = _PDB_CHAIN_IDENTIFIERS[chain_index % len(_PDB_CHAIN_IDENTIFIERS)]
            lattice_chain = (
                chains_per_layer - 1 - surface_chain
                if surface_order_reversed
                else surface_chain
            )
            is_center = bool((crystal_layer + lattice_chain) % 2)
            template = center_template if is_center else origin_template
            template_labels = center_labels if is_center else origin_labels
            site = crystal_layer * layer_step + lattice_chain * in_plane_vector
            base_translation = site - layer_step if is_center else site
            for chain_cell in range(chain_repeats):
                translation = base_translation + chain_cell * vectors[2]
                coordinate_blocks.append(template + translation)
                atom_labels.extend(template_labels)
                chain_ids.extend([chain_id] * len(template))
                chain_numbers.extend([chain_number] * len(template))
                layer_numbers.extend([output_layer + 1] * len(template))
                first_residue = 2 * chain_cell + 1
                residue_numbers.extend([first_residue] * residue_atoms)
                residue_numbers.extend([first_residue + 1] * residue_atoms)

    projected = np.vstack(coordinate_blocks) @ axes.T
    projected[:, 0] -= float(np.min(projected[:, 0]))
    projected[:, 1] -= float(np.min(projected[:, 1])) - padding_z / 2.0
    projected[:, 2] -= float(np.min(projected[:, 2])) - vacuum / 2.0
    slab_extent = float(np.ptp(projected[:, 2]))
    chain_length = float(np.linalg.norm(vectors[2])) * chain_repeats
    basis = np.diag(
        [in_plane_length * chains_per_layer, chain_length + padding_z, slab_extent + vacuum]
    )
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
        layer_numbers=tuple(layer_numbers),
        annotations=(
            f"shape=slab face=({face.h},{face.k},{face.l}) layers={layers}",
            f"chains_per_layer={chains_per_layer} total_chains={chain_count}",
            f"degree_polymerization={degree_polymerization}",
            f"chain_boundary={chain_boundary}",
            f"d_hkl={face.spacing(allomorph.unit_cell):.6f} angstrom",
            f"vacuum={vacuum:.6f} angstrom",
            f"padding_z={padding_z:.6f} angstrom total longitudinal padding",
            "surface normal aligned with Cartesian -z; chains aligned with Cartesian -y",
            "orientation_transform=proper_rotation determinant=+1",
            "chain numbering starts at requested low-z face; within each layer low-x to high-x",
        ),
        residue_name=allomorph.residue_name,
    )


def build_slab_by_dimensions(  # noqa: PLR0913
    allomorph: AllomorphDefinition,
    face: MillerIndex,
    *,
    width: float,
    thickness: float,
    length: float,
    vacuum: float,
    padding_z: float = 0.0,
    chain_boundary: str = "finite",
) -> BulkCrystal:
    """Build a slab meeting or exceeding requested physical dimensions.

    Dimensions are in angstrom. Width follows Cartesian x, length follows the
    cellulose chains along Cartesian y, and thickness follows the surface-normal
    Cartesian z direction. Counts are rounded upward so no chain, molecular
    layer, or crystallographic cellobiose repeat is cut.
    """
    return build_slab_with_sizing(
        allomorph,
        face,
        width=width,
        thickness=thickness,
        length=length,
        vacuum=vacuum,
        padding_z=padding_z,
        chain_boundary=chain_boundary,
    )


def build_slab_with_sizing(  # noqa: PLR0913
    allomorph: AllomorphDefinition,
    face: MillerIndex,
    *,
    layers: int | None = None,
    chains_per_layer: int | None = None,
    degree_polymerization: int | None = None,
    width: float | None = None,
    thickness: float | None = None,
    length: float | None = None,
    vacuum: float,
    padding_z: float = 0.0,
    chain_boundary: str = "finite",
) -> BulkCrystal:
    """Build a slab using an independent discrete or physical choice per axis."""
    pairs = (
        ("chains_per_layer", chains_per_layer, "width", width),
        ("layers", layers, "thickness", thickness),
        ("degree_polymerization", degree_polymerization, "length", length),
    )
    for count_name, count, size_name, size in pairs:
        if (count is None) == (size is None):
            msg = f"supply exactly one of {count_name} or {size_name}"
            raise ValueError(msg)
    for name, value in (("width", width), ("thickness", thickness), ("length", length)):
        if value is not None and (not np.isfinite(value) or value <= 0):
            msg = f"slab {name} must be a finite positive angstrom value"
            raise ValueError(msg)
    if not np.isfinite(vacuum) or vacuum < 0:
        msg = "vacuum must be a finite, nonnegative angstrom value"
        raise ValueError(msg)

    chain_repeat_length = float(np.linalg.norm(allomorph.unit_cell.vectors()[allomorph.chain_axis]))
    if degree_polymerization is None:
        assert length is not None
        chain_repeats = max(1, int(np.ceil(length / chain_repeat_length)))
        degree_polymerization = 2 * chain_repeats

    def probe(layers: int, chains: int) -> BulkCrystal:
        return build_chain_parallel_slab(
            allomorph,
            face,
            layers=layers,
            chains_per_layer=chains,
            degree_polymerization=2,
            vacuum=0.0,
            padding_z=0.0,
            chain_boundary=chain_boundary,
        )

    if chains_per_layer is None:
        assert width is not None
        special_i_beta_110 = allomorph.name == "I-beta" and face.values in {
            (1, 1, 0),
            (-1, -1, 0),
        }
        chains_per_layer = 2 if special_i_beta_110 else 1
        while float(probe(1, chains_per_layer).basis_vectors[0, 0]) + _LAYER_TOLERANCE < width:
            chains_per_layer += 2 if special_i_beta_110 else 1
            if chains_per_layer > _MAX_SLAB_DISCRETE_COUNT:
                msg = "requested slab width exceeds the supported discrete-chain search limit"
                raise ValueError(msg)

    if layers is None:
        assert thickness is not None
        layers = 1
        while (
            float(probe(layers, chains_per_layer).basis_vectors[2, 2]) + _LAYER_TOLERANCE
            < thickness
        ):
            layers += 1
            if layers > _MAX_SLAB_DISCRETE_COUNT:
                msg = "requested slab thickness exceeds the supported molecular-layer search limit"
                raise ValueError(msg)

    slab = build_chain_parallel_slab(
        allomorph,
        face,
        layers=layers,
        chains_per_layer=chains_per_layer,
        degree_polymerization=degree_polymerization,
        vacuum=vacuum,
        padding_z=padding_z,
        chain_boundary=chain_boundary,
    )
    realized_width = float(slab.basis_vectors[0, 0])
    realized_thickness = float(slab.basis_vectors[2, 2]) - vacuum
    realized_length = float(slab.basis_vectors[1, 1]) - padding_z
    physical_axes = sum(value is not None for value in (width, thickness, length))
    size_mode = (
        "physical_dimensions"
        if physical_axes == len(pairs)
        else "exact_counts"
        if physical_axes == 0
        else "mixed"
    )
    dimension_annotations = tuple(
        annotation
        for annotation in (
            (
                f"requested_width={width:.6f} angstrom "
                f"realized_width={realized_width:.6f} angstrom"
                if width is not None
                else None
            ),
            (
                f"requested_thickness={thickness:.6f} angstrom "
                f"realized_thickness={realized_thickness:.6f} angstrom"
                if thickness is not None
                else None
            ),
            (
                f"requested_length={length:.6f} angstrom "
                f"realized_length={realized_length:.6f} angstrom"
                if length is not None
                else None
            ),
        )
        if annotation is not None
    )
    return replace(
        slab,
        annotations=(
            *slab.annotations,
            f"size_mode={size_mode}",
            *dimension_annotations,
        ),
    )


def _build_exploratory_chain_parallel_slab(  # noqa: PLR0912, PLR0913, PLR0915
    allomorph: AllomorphDefinition,
    face: MillerIndex,
    *,
    layers: int,
    chains_per_layer: int,
    degree_polymerization: int,
    vacuum: float,
    padding_z: float,
    chain_boundary: str,
) -> BulkCrystal:
    """Build a chain-parallel slab by selecting complete lattice chains.

    Candidate crystallographic chains are generated across the two lattice
    axes transverse to the polymer direction, grouped into planes normal to
    the requested Miller index, and cut as complete chains. This path supports
    all allomorph chain axes and preserves parallel or antiparallel chain
    directions from the crystallographic definition.
    """
    vectors = allomorph.unit_cell.vectors()
    normal = face.normal(allomorph.unit_cell)
    chain_vector = vectors[allomorph.chain_axis]
    chain_axis = chain_vector / np.linalg.norm(chain_vector)
    surface_x = np.cross(chain_axis, normal)
    surface_x /= np.linalg.norm(surface_x)
    axes = _proper_slab_axes(surface_x, chain_axis, normal)

    unit_cell = build_unit_cell(allomorph)
    unit_coordinates = _unwrap_unit_cell_chains(
        unit_cell,
        vectors,
        chain_axis=allomorph.chain_axis,
    )
    unit_chain_order = tuple(dict.fromkeys(unit_cell.chain_numbers))
    templates = tuple(
        unit_coordinates[np.asarray(unit_cell.chain_numbers) == chain_number]
        for chain_number in unit_chain_order
    )
    template_labels = tuple(
        tuple(
            label
            for label, number in zip(unit_cell.atom_labels, unit_cell.chain_numbers, strict=True)
            if number == chain_number
        )
        for chain_number in unit_chain_order
    )
    template_residues = tuple(
        tuple(
            residue
            for residue, number in zip(
                unit_cell.residue_numbers, unit_cell.chain_numbers, strict=True
            )
            if number == chain_number
        )
        for chain_number in unit_chain_order
    )
    direction_values = tuple(
        dict(
            zip(
                allomorph.asymmetric_unit.chain_labels,
                allomorph.asymmetric_unit.chain_directions,
                strict=True,
            )
        ).values()
    )
    if len(direction_values) != len(templates):
        msg = "chain direction metadata does not match crystallographic chain templates"
        raise ValueError(msg)
    transverse_axes = tuple(axis for axis in range(3) if axis != allomorph.chain_axis)

    search_radius = max(8, layers + chains_per_layer + abs(face.h) + abs(face.k))
    candidates: list[tuple[float, float, int, int, int]] = []
    for i in range(-search_radius, search_radius + 1):
        for j in range(-search_radius, search_radius + 1):
            translation = i * vectors[transverse_axes[0]] + j * vectors[transverse_axes[1]]
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
    ideal_start = max(0, min(middle - layers // 2, len(eligible_rows) - layers))
    requested_termination = (
        _I_BETA_TERMINATION_BY_FACE.get(face.values) if allomorph.name == "I-beta" else None
    )
    if requested_termination is None:
        start = ideal_start
    else:
        possible_starts = [
            candidate_start
            for candidate_start in range(len(eligible_rows) - layers + 1)
            if any(
                candidate[2] == requested_termination
                for candidate in eligible_rows[candidate_start + layers - 1]
            )
        ]
        if not possible_starts:
            chain_type = _I_BETA_CHAIN_TYPE_NAMES[requested_termination]
            msg = f"could not expose the required {chain_type}-chain termination for {face.values}"
            raise ValueError(msg)
        start = min(possible_starts, key=lambda value: (abs(value - ideal_start), value))
    selected_rows = eligible_rows[start : start + layers]
    selected: list[tuple[float, float, int, int, int]] = []
    for row in reversed(selected_rows):
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
    layer_numbers: list[int] = []
    for chain_index, (_, _, template_index, i, j) in enumerate(selected):
        chain_number = chain_index + 1
        chain_id = _PDB_CHAIN_IDENTIFIERS[chain_index % len(_PDB_CHAIN_IDENTIFIERS)]
        transverse_translation = i * vectors[transverse_axes[0]] + j * vectors[transverse_axes[1]]
        template = templates[template_index]
        labels = template_labels[template_index]
        residues = template_residues[template_index]
        direction = direction_values[template_index]
        for chain_cell in range(chain_repeats):
            translation = transverse_translation + chain_cell * direction * chain_vector
            coordinate_blocks.append(template + translation)
            atom_labels.extend(labels)
            chain_ids.extend([chain_id] * len(template))
            chain_numbers.extend([chain_number] * len(template))
            layer_numbers.extend(
                [chain_index // chains_per_layer + 1] * len(template)
            )
            residue_numbers.extend(residue + 2 * chain_cell for residue in residues)

    projected = np.vstack(coordinate_blocks) @ axes.T
    projected[:, 0] -= float(np.min(projected[:, 0]))
    projected[:, 1] -= float(np.min(projected[:, 1])) - padding_z / 2.0
    projected[:, 2] -= float(np.min(projected[:, 2])) - vacuum / 2.0
    extents = np.ptp(projected, axis=0)
    chain_length = float(np.linalg.norm(chain_vector)) * chain_repeats
    transverse_indices = (face.values[transverse_axes[0]], face.values[transverse_axes[1]])
    divisor = int(np.gcd(abs(transverse_indices[0]), abs(transverse_indices[1])))
    lattice_surface_vector = (transverse_indices[1] // divisor) * vectors[transverse_axes[0]] - (
        transverse_indices[0] // divisor
    ) * vectors[transverse_axes[1]]
    lattice_surface_repeat = abs(float(np.dot(lattice_surface_vector, surface_x)))
    selected_x_by_row = [
        sorted(candidate[1] for candidate in selected[index : index + chains_per_layer])
        for index in range(0, len(selected), chains_per_layer)
    ]
    spacings = [
        right - left
        for row in selected_x_by_row
        for left, right in pairwise(row)
        if right - left > _LAYER_TOLERANCE
    ]
    surface_repeat_length = float(np.median(spacings)) if spacings else lattice_surface_repeat
    periodic_width = surface_repeat_length * chains_per_layer
    basis = np.diag([periodic_width, chain_length + padding_z, float(extents[2]) + vacuum])
    termination_annotations: tuple[str, ...] = ()
    if requested_termination is not None:
        chain_type = _I_BETA_CHAIN_TYPE_NAMES[requested_termination]
        termination_annotations = (
            f"requested_face_chain_type={chain_type}",
            "requested face is the low-z surface",
        )
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
        layer_numbers=tuple(layer_numbers),
        annotations=(
            f"shape=slab face=({face.h},{face.k},{face.l}) layers={layers}",
            f"chains_per_layer={chains_per_layer} total_chains={len(selected)}",
            f"degree_polymerization={degree_polymerization}",
            f"chain_boundary={chain_boundary}",
            f"d_hkl={face.spacing(allomorph.unit_cell):.6f} angstrom",
            f"vacuum={vacuum:.6f} angstrom",
            f"padding_z={padding_z:.6f} angstrom total longitudinal padding",
            f"lateral_repeat={surface_repeat_length:.6f} angstrom",
            "surface normal aligned with Cartesian -z; chains aligned with Cartesian -y",
            "chain numbering starts at requested low-z face; within each layer low-x to high-x",
            "construction=whole_chain_reciprocal_lattice_cut",
            *termination_annotations,
        ),
        residue_name=allomorph.residue_name,
    )
