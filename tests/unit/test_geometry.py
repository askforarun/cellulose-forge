"""Tests for I-beta unit-cell and legacy-compatible bulk geometry."""

from collections import Counter
from dataclasses import replace

import numpy as np
import pytest
from numpy.testing import assert_allclose

from cellulose_forge.crystallography import MillerIndex, fractional_to_cartesian, load_allomorph
from cellulose_forge.geometry import (
    BulkCrystal,
    build_36_chain_fibril,
    build_bulk_crystal,
    build_chain_parallel_slab,
    build_slab_by_dimensions,
    build_slab_with_sizing,
    build_unit_cell,
    validate_cellulose_chirality,
)


def _c1_signed_volumes(crystal: BulkCrystal) -> np.ndarray:
    """Return C1 stereocentre volumes using O5, C2, and H1 substituents."""
    labels = np.asarray(crystal.atom_labels)
    chains = np.asarray(crystal.chain_numbers)
    residues = np.asarray(crystal.residue_numbers)
    volumes: list[float] = []
    for chain, residue in sorted(set(zip(chains, residues, strict=True))):
        mask = (chains == chain) & (residues == residue)
        lookup = dict(
            zip(labels[mask], crystal.cartesian_coordinates[mask], strict=True)
        )
        if not {"C1", "O5", "C2", "H1"}.issubset(lookup):
            continue
        center = lookup["C1"]
        volumes.append(
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
    return np.asarray(volumes)


@pytest.mark.parametrize(
    ("allomorph_key", "face"),
    [
        ("i_alpha", (0, 1, 0)),
        ("i_alpha", (0, -1, 0)),
        ("i_beta", (1, 1, 0)),
        ("i_beta", (-1, -1, 0)),
        ("ii", (1, 1, 0)),
        ("ii", (-1, -1, 0)),
        ("iii_i", (1, 1, 0)),
        ("iii_i", (-1, -1, 0)),
    ],
)
def test_slab_orientation_preserves_glucose_chirality(
    allomorph_key: str,
    face: tuple[int, int, int],
) -> None:
    definition = load_allomorph(allomorph_key)
    source_volumes = _c1_signed_volumes(build_unit_cell(definition))
    slab_volumes = _c1_signed_volumes(
        build_chain_parallel_slab(
            definition,
            MillerIndex(*face),
            layers=2,
            chains_per_layer=4,
            degree_polymerization=2,
            vacuum=0.0,
        )
    )

    assert source_volumes.size > 0
    assert slab_volumes.size > 0
    assert np.all(np.abs(source_volumes) > 1e-6)
    assert np.all(np.abs(slab_volumes) > 1e-6)
    assert set(np.sign(slab_volumes)) == set(np.sign(source_volumes))


@pytest.mark.parametrize("allomorph_key", ["i_alpha", "i_beta", "ii", "iii_i"])
def test_chirality_validation_reports_d_cellulose(allomorph_key: str) -> None:
    result = validate_cellulose_chirality(build_unit_cell(load_allomorph(allomorph_key)))

    assert result["status"] == "passed"
    assert result["configuration"] == "D-cellulose"
    assert result["expected_signed_volume_sign"] == "negative"
    assert result["residue_count"] > 0


def test_chirality_validation_rejects_mirrored_geometry() -> None:
    crystal = build_unit_cell(load_allomorph("i_beta"))
    mirrored = replace(
        crystal,
        cartesian_coordinates=crystal.cartesian_coordinates * np.asarray([-1.0, 1.0, 1.0]),
        basis_vectors=crystal.basis_vectors * np.asarray([-1.0, 1.0, 1.0]),
    )

    with pytest.raises(ValueError, match="mirrored C1 geometry"):
        validate_cellulose_chirality(mirrored)


@pytest.mark.parametrize(
    ("allomorph_key", "atoms_per_residue"),
    [("i_alpha", 21), ("i_beta", 21), ("ii", 18)],
)
def test_literature_36_chain_fibrils_have_complete_dp_chains(
    allomorph_key: str,
    atoms_per_residue: int,
) -> None:
    fibril = build_36_chain_fibril(
        load_allomorph(allomorph_key),
        degree_polymerization=6,
        vacuum=20.0,
        chain_boundary="periodic",
    )

    assert fibril.atom_count == 36 * 6 * atoms_per_residue
    assert set(fibril.chain_numbers) == set(range(1, 37))
    assert set(fibril.residue_numbers) == set(range(1, 7))
    assert all(
        sum(chain_number == chain for chain_number in fibril.chain_numbers) == 6 * atoms_per_residue
        for chain in range(1, 37)
    )
    assert any("fibril_model=36-chain" in value for value in fibril.annotations)
    assert any("degree_polymerization=6" in value for value in fibril.annotations)
    assert any("legacy_fragment_ids=" in value for value in fibril.annotations)


def test_36_chain_fibril_is_deterministic() -> None:
    definition = load_allomorph("i_beta")
    first = build_36_chain_fibril(definition, degree_polymerization=4, vacuum=12.0)
    second = build_36_chain_fibril(definition, degree_polymerization=4, vacuum=12.0)

    assert first.atom_labels == second.atom_labels
    assert first.chain_numbers == second.chain_numbers
    assert_allclose(first.cartesian_coordinates, second.cartesian_coordinates)


def test_36_chain_fibril_supports_independent_finite_box_padding() -> None:
    definition = load_allomorph("i_beta")
    unpadded = build_36_chain_fibril(
        definition,
        degree_polymerization=4,
        vacuum=0.0,
        chain_boundary="finite",
    )
    padded = build_36_chain_fibril(
        definition,
        degree_polymerization=4,
        vacuum=0.0,
        padding_x=10.0,
        padding_y=20.0,
        padding_z=30.0,
        chain_boundary="finite",
    )

    assert_allclose(
        np.diag(padded.basis_vectors) - np.diag(unpadded.basis_vectors),
        [10.0, 30.0, 20.0],
    )
    assert_allclose(np.min(padded.cartesian_coordinates, axis=0), [5.0, 15.0, 10.0])
    assert any("padding_x=10.000000" in value for value in padded.annotations)
    assert any("padding_y=20.000000" in value for value in padded.annotations)
    assert any("padding_z=30.000000" in value for value in padded.annotations)


def test_36_chain_periodic_fibril_rejects_longitudinal_padding() -> None:
    with pytest.raises(ValueError, match=r"padding_z must be zero.*periodic"):
        build_36_chain_fibril(
            load_allomorph("i_beta"),
            degree_polymerization=4,
            vacuum=10.0,
            padding_z=5.0,
            chain_boundary="periodic",
        )


def test_36_chain_fibril_rejects_unsupported_iii_i() -> None:
    with pytest.raises(ValueError, match=r"not defined.*reference Cellulose-Builder"):
        build_36_chain_fibril(
            load_allomorph("iii_i"),
            degree_polymerization=4,
            vacuum=10.0,
        )


def test_36_chain_fibril_rejects_odd_degree_polymerization() -> None:
    with pytest.raises(ValueError, match=r"requires an even degree of polymerization"):
        build_36_chain_fibril(
            load_allomorph("i_beta"),
            degree_polymerization=5,
            vacuum=10.0,
        )


def test_i_beta_unit_cell_has_two_chains_and_four_residues() -> None:
    definition = load_allomorph("i_beta")

    crystal = build_unit_cell(definition)

    assert crystal.atom_count == 84
    assert Counter(crystal.chain_ids) == {"A": 42, "B": 42}
    assert set(zip(crystal.chain_ids, crystal.residue_numbers, strict=True)) == {
        ("A", 1),
        ("A", 2),
        ("B", 1),
        ("B", 2),
    }


def test_i_alpha_unit_cell_has_one_chain_and_two_residues() -> None:
    crystal = build_unit_cell(load_allomorph("i_alpha"))

    assert crystal.atom_count == 42
    assert set(crystal.chain_numbers) == {1}
    assert set(crystal.residue_numbers) == {1, 2}
    assert Counter(crystal.residue_numbers) == {1: 21, 2: 21}


def test_i_alpha_bulk_replicates_chains_along_crystallographic_a() -> None:
    crystal = build_bulk_crystal(load_allomorph("i_alpha"), (3, 2, 4))

    assert crystal.atom_count == 42 * 3 * 2 * 4
    assert len(set(crystal.chain_numbers)) == 2 * 4
    assert set(crystal.residue_numbers) == set(range(1, 7))


def test_iii_i_unit_cell_and_bulk_follow_p21_chain_geometry() -> None:
    definition = load_allomorph("iii_i")
    unit_cell = build_unit_cell(definition)
    bulk = build_bulk_crystal(definition, (2, 3, 4))

    assert unit_cell.atom_count == 42
    assert set(unit_cell.chain_numbers) == {1}
    assert set(unit_cell.residue_numbers) == {1, 2}
    assert bulk.atom_count == 42 * 2 * 3 * 4
    assert len(set(bulk.chain_numbers)) == 2 * 3
    assert set(bulk.residue_numbers) == set(range(1, 9))


def test_ii_unit_cell_and_bulk_keep_two_antiparallel_chain_identities() -> None:
    definition = load_allomorph("ii")
    unit_cell = build_unit_cell(definition)
    bulk = build_bulk_crystal(definition, (2, 3, 4))

    assert unit_cell.atom_count == 72
    assert set(zip(unit_cell.chain_numbers, unit_cell.residue_numbers, strict=True)) == {
        (1, 1),
        (1, 2),
        (2, 1),
        (2, 2),
    }
    assert bulk.atom_count == 72 * 2 * 3 * 4
    assert len(set(bulk.chain_numbers)) == 2 * 2 * 3
    assert set(bulk.residue_numbers) == set(range(1, 9))


def test_i_beta_unit_cell_coordinates_follow_stored_symmetry() -> None:
    definition = load_allomorph("i_beta")
    asymmetric = definition.asymmetric_unit.fractional_coordinates
    crystal = build_unit_cell(definition)

    expected_identity = fractional_to_cartesian(asymmetric, definition.unit_cell)
    expected_screw = fractional_to_cartesian(
        definition.symmetry_operations[1].apply(asymmetric), definition.unit_cell
    )

    assert_allclose(crystal.cartesian_coordinates[:42], expected_identity, atol=1e-12)
    assert_allclose(crystal.cartesian_coordinates[42:], expected_screw, atol=1e-12)


def test_bulk_replication_retains_every_complete_unit_cell() -> None:
    definition = load_allomorph("i_beta")

    crystal = build_bulk_crystal(definition, (2, 3, 4))

    assert crystal.atom_count == 84 * 2 * 3 * 4
    assert crystal.replication == (2, 3, 4)
    assert_allclose(
        crystal.basis_vectors,
        definition.unit_cell.vectors() * [[2.0], [3.0], [4.0]],
    )
    assert len(set(crystal.chain_ids)) == 2 * 2 * 3
    assert set(crystal.residue_numbers) == set(range(1, 9))


def test_bulk_replication_translates_coordinates_by_lattice_vectors() -> None:
    definition = load_allomorph("i_beta")
    crystal = build_bulk_crystal(definition, (2, 1, 2))

    atoms_per_cell = build_unit_cell(definition).atom_count
    first_cell = crystal.cartesian_coordinates[:atoms_per_cell]
    assert_allclose(
        crystal.cartesian_coordinates[atoms_per_cell : 2 * atoms_per_cell],
        first_cell + definition.unit_cell.vectors()[0],
    )
    assert_allclose(
        crystal.cartesian_coordinates[2 * atoms_per_cell : 3 * atoms_per_cell],
        first_cell + definition.unit_cell.vectors()[2],
    )


def test_bulk_replication_unwraps_center_chain_across_transverse_boundary() -> None:
    crystal = build_bulk_crystal(load_allomorph("i_beta"), (1, 1, 1))
    labels = np.asarray(crystal.atom_labels)
    chains = np.asarray(crystal.chain_numbers)
    residues = np.asarray(crystal.residue_numbers)

    for chain in (1, 2):
        first = crystal.cartesian_coordinates[
            (chains == chain) & (residues == 1) & (labels == "C1")
        ][0]
        second = crystal.cartesian_coordinates[
            (chains == chain) & (residues == 2) & (labels == "C1")
        ][0]
        transverse_displacement = second[:2] - first[:2]
        assert np.linalg.norm(transverse_displacement) < 1.0


def test_bulk_replication_rejects_invalid_counts() -> None:
    definition = load_allomorph("i_beta")

    for replication in ((0, 1, 1), (1, -1, 1), (1, 1, 1.5)):
        try:
            build_bulk_crystal(definition, replication)  # type: ignore[arg-type]
        except ValueError as error:
            assert "three positive integers" in str(error)
        else:
            raise AssertionError("invalid replication was accepted")


def test_i_beta_110_slab_is_surface_aligned_and_keeps_complete_residues() -> None:
    definition = load_allomorph("i_beta")

    slab = build_chain_parallel_slab(
        definition,
        MillerIndex(1, 1, 0),
        layers=4,
        chains_per_layer=20,
        degree_polymerization=12,
        vacuum=30.0,
    )

    assert slab.atom_count == 21 * 4 * 20 * 12
    assert slab.atom_count % definition.asymmetric_unit.residue_atom_count == 0
    assert len(set(slab.chain_numbers)) == 80
    assert set(slab.chain_numbers) == set(range(1, 81))
    assert set(slab.residue_numbers) == set(range(1, 13))
    assert_allclose(slab.basis_vectors - np.diag(np.diag(slab.basis_vectors)), 0.0)
    assert slab.basis_vectors[2, 2] == pytest.approx(
        np.ptp(slab.cartesian_coordinates[:, 2]) + 30.0
    )


def test_i_beta_110_slab_rejects_incomplete_transverse_lattice_repeat() -> None:
    with pytest.raises(ValueError, match="even number of chains per layer"):
        build_chain_parallel_slab(
            load_allomorph("i_beta"),
            MillerIndex(1, 1, 0),
            layers=2,
            chains_per_layer=3,
            degree_polymerization=2,
            vacuum=10.0,
        )


def test_i_beta_110_chain_types_alternate_within_each_layer() -> None:
    definition = load_allomorph("i_beta")
    slab = build_chain_parallel_slab(
        definition,
        MillerIndex(1, 1, 0),
        layers=2,
        chains_per_layer=4,
        degree_polymerization=2,
        vacuum=10.0,
    )
    unit_cell = build_unit_cell(definition)
    residue_atoms = definition.asymmetric_unit.residue_atom_count
    origin_indices = np.r_[0:residue_atoms, 2 * residue_atoms : 3 * residue_atoms]
    center_indices = np.r_[residue_atoms : 2 * residue_atoms, 3 * residue_atoms : 4 * residue_atoms]
    origin = unit_cell.cartesian_coordinates[origin_indices]
    center = unit_cell.cartesian_coordinates[center_indices].copy()
    center[residue_atoms:] += definition.unit_cell.vectors()[0] + definition.unit_cell.vectors()[1]

    def distance_matrix(coordinates: np.ndarray) -> np.ndarray:
        differences = coordinates[:, None, :] - coordinates[None, :, :]
        return np.linalg.norm(differences, axis=2)

    references = (distance_matrix(origin), distance_matrix(center))
    chain_numbers = np.asarray(slab.chain_numbers)
    chain_types: list[int] = []
    for chain in range(1, 9):
        coordinates = slab.cartesian_coordinates[chain_numbers == chain]
        signature = distance_matrix(coordinates)
        errors = [float(np.max(np.abs(signature - reference))) for reference in references]
        chain_types.append(int(np.argmin(errors)))
        assert min(errors) < 1e-8

    assert chain_types[:4] in ([0, 1, 0, 1], [1, 0, 1, 0])
    assert chain_types[4:] in ([0, 1, 0, 1], [1, 0, 1, 0])


def test_i_beta_110_slab_has_no_short_interchain_contacts() -> None:
    slab = build_chain_parallel_slab(
        load_allomorph("i_beta"),
        MillerIndex(1, 1, 0),
        layers=2,
        chains_per_layer=4,
        degree_polymerization=2,
        vacuum=10.0,
        chain_boundary="periodic",
    )
    coordinates = slab.cartesian_coordinates
    chains = np.asarray(slab.chain_numbers)
    box = np.diag(slab.basis_vectors)
    minimum = np.inf
    for first in range(1, 9):
        first_coordinates = coordinates[chains == first]
        for second in range(first + 1, 9):
            differences = first_coordinates[:, None, :] - coordinates[chains == second][None, :, :]
            differences -= np.round(differences / box) * box
            minimum = min(minimum, float(np.min(np.linalg.norm(differences, axis=2))))

    assert minimum >= 1.7


def test_slab_chain_numbering_starts_at_requested_face_and_increases_with_x() -> None:
    definition = load_allomorph("i_beta")
    chains_per_layer = 4
    slab = build_chain_parallel_slab(
        definition,
        MillerIndex(1, 1, 0),
        layers=3,
        chains_per_layer=chains_per_layer,
        degree_polymerization=2,
        vacuum=10.0,
    )
    chain_numbers = np.asarray(slab.chain_numbers)
    centers = {
        chain: np.mean(slab.cartesian_coordinates[chain_numbers == chain], axis=0)
        for chain in range(1, 13)
    }

    first_layer_z = [centers[chain][2] for chain in range(1, 5)]
    second_layer_z = [centers[chain][2] for chain in range(5, 9)]
    assert max(first_layer_z) < min(second_layer_z)
    assert [centers[chain][0] for chain in range(1, 5)] == sorted(
        centers[chain][0] for chain in range(1, 5)
    )


def test_initial_slab_builder_rejects_chain_cutting_face() -> None:
    definition = load_allomorph("i_beta")

    with pytest.raises(ValueError, match=r"parallel.*l = 0"):
        build_chain_parallel_slab(
            definition,
            MillerIndex(1, 1, 1),
            layers=2,
            chains_per_layer=4,
            degree_polymerization=2,
            vacuum=10.0,
        )


@pytest.mark.parametrize(
    "face",
    [
        (1, 0, 0),
        (2, 0, 0),
        (-1, 0, 0),
        (0, 1, 0),
        (0, 2, 0),
        (0, -1, 0),
        (1, -1, 0),
        (-1, 1, 0),
        (-1, -1, 0),
    ],
)
def test_i_beta_chain_parallel_faces_make_complete_aligned_slabs(
    face: tuple[int, int, int],
) -> None:
    definition = load_allomorph("i_beta")
    miller = MillerIndex(*face)
    chains_per_layer = 6 if face == (-1, -1, 0) else 5

    slab = build_chain_parallel_slab(
        definition,
        miller,
        layers=3,
        chains_per_layer=chains_per_layer,
        degree_polymerization=4,
        vacuum=12.0,
    )

    assert slab.atom_count == 21 * 3 * chains_per_layer * 4
    assert len(set(slab.chain_numbers)) == 3 * chains_per_layer
    assert set(slab.residue_numbers) == {1, 2, 3, 4}
    assert_allclose(slab.basis_vectors - np.diag(np.diag(slab.basis_vectors)), 0.0)
    assert slab.basis_vectors[2, 2] == pytest.approx(
        np.ptp(slab.cartesian_coordinates[:, 2]) + 12.0
    )


@pytest.mark.parametrize(
    ("allomorph_key", "face"),
    [
        ("i_alpha", (0, 1, 0)),
        ("ii", (1, 1, 0)),
        ("iii_i", (1, 1, 0)),
    ],
)
def test_allomorph_chain_parallel_slabs_preserve_complete_chains_and_layers(
    allomorph_key: str,
    face: tuple[int, int, int],
) -> None:
    definition = load_allomorph(allomorph_key)
    layers = 3
    chains_per_layer = 5
    degree_polymerization = 4

    slab = build_chain_parallel_slab(
        definition,
        MillerIndex(*face),
        layers=layers,
        chains_per_layer=chains_per_layer,
        degree_polymerization=degree_polymerization,
        vacuum=12.0,
    )

    expected_atoms = (
        definition.asymmetric_unit.residue_atom_count
        * layers
        * chains_per_layer
        * degree_polymerization
    )
    assert slab.atom_count == expected_atoms
    assert set(slab.chain_numbers) == set(range(1, layers * chains_per_layer + 1))
    assert set(slab.residue_numbers) == set(range(1, degree_polymerization + 1))
    assert_allclose(slab.basis_vectors - np.diag(np.diag(slab.basis_vectors)), 0.0)
    chain_numbers = np.asarray(slab.chain_numbers)
    layer_centers = [
        np.mean(
            [
                np.mean(slab.cartesian_coordinates[chain_numbers == chain, 2])
                for chain in range(
                    layer * chains_per_layer + 1,
                    (layer + 1) * chains_per_layer + 1,
                )
            ]
        )
        for layer in range(layers)
    ]
    assert layer_centers == sorted(layer_centers)
    assert "construction=whole_chain_reciprocal_lattice_cut" in slab.annotations


def test_finite_slab_chain_end_padding_expands_and_centers_box() -> None:
    definition = load_allomorph("i_alpha")
    slab = build_chain_parallel_slab(
        definition,
        MillerIndex(0, 1, 0),
        layers=2,
        chains_per_layer=3,
        degree_polymerization=4,
        vacuum=10.0,
        padding_z=12.0,
        chain_boundary="finite",
    )

    assert slab.basis_vectors[1, 1] == pytest.approx(32.8)
    assert np.min(slab.cartesian_coordinates[:, 1]) == pytest.approx(6.0)
    assert any(annotation.startswith("padding_z=12.000000") for annotation in slab.annotations)


def test_periodic_slab_rejects_chain_end_padding() -> None:
    definition = load_allomorph("ii")
    with pytest.raises(ValueError, match="padding_z must be zero"):
        build_chain_parallel_slab(
            definition,
            MillerIndex(1, 1, 0),
            layers=2,
            chains_per_layer=3,
            degree_polymerization=4,
            vacuum=10.0,
            padding_z=2.0,
            chain_boundary="periodic",
        )


def test_parallel_miller_indices_retain_distinct_plane_spacings() -> None:
    definition = load_allomorph("i_beta")
    slab_100 = build_chain_parallel_slab(
        definition,
        MillerIndex(1, 0, 0),
        layers=2,
        chains_per_layer=3,
        degree_polymerization=2,
        vacuum=0.0,
    )
    slab_200 = build_chain_parallel_slab(
        definition,
        MillerIndex(2, 0, 0),
        layers=2,
        chains_per_layer=3,
        degree_polymerization=2,
        vacuum=0.0,
    )

    spacing_100 = MillerIndex(1, 0, 0).spacing(definition.unit_cell)
    spacing_200 = MillerIndex(2, 0, 0).spacing(definition.unit_cell)
    assert spacing_200 == pytest.approx(spacing_100 / 2.0)
    assert f"d_hkl={spacing_100:.6f} angstrom" in slab_100.annotations
    assert f"d_hkl={spacing_200:.6f} angstrom" in slab_200.annotations
    assert "requested_face_chain_type=origin" in slab_100.annotations
    assert "requested_face_chain_type=center" in slab_200.annotations
    assert not np.allclose(slab_100.cartesian_coordinates, slab_200.cartesian_coordinates)


def test_010_and_020_select_origin_and_center_terminations() -> None:
    definition = load_allomorph("i_beta")

    slab_010 = build_chain_parallel_slab(
        definition,
        MillerIndex(0, 1, 0),
        layers=3,
        chains_per_layer=4,
        degree_polymerization=2,
        vacuum=0.0,
    )
    slab_020 = build_chain_parallel_slab(
        definition,
        MillerIndex(0, 2, 0),
        layers=3,
        chains_per_layer=4,
        degree_polymerization=2,
        vacuum=0.0,
    )

    assert "requested_face_chain_type=origin" in slab_010.annotations
    assert "requested_face_chain_type=center" in slab_020.annotations
    assert not np.allclose(slab_010.cartesian_coordinates, slab_020.cartesian_coordinates)


@pytest.mark.parametrize(
    ("face", "vector_index"),
    [((1, 0, 0), 1), ((2, 0, 0), 1), ((0, 1, 0), 0), ((0, 2, 0), 0)],
)
def test_axis_face_lateral_box_has_uniform_periodic_chain_spacing(
    face: tuple[int, int, int], vector_index: int
) -> None:
    definition = load_allomorph("i_beta")
    chains_per_layer = 8

    slab = build_chain_parallel_slab(
        definition,
        MillerIndex(*face),
        layers=3,
        chains_per_layer=chains_per_layer,
        degree_polymerization=2,
        vacuum=0.0,
    )

    expected_spacing = np.linalg.norm(definition.unit_cell.vectors()[vector_index])
    chain_numbers = np.asarray(slab.chain_numbers)
    first_layer_centers = np.asarray(
        [
            np.mean(slab.cartesian_coordinates[chain_numbers == chain], axis=0)
            for chain in range(1, chains_per_layer + 1)
        ]
    )
    x_centers = np.sort(first_layer_centers[:, 0])
    assert_allclose(np.diff(x_centers), expected_spacing, atol=1e-6)
    boundary_spacing = slab.basis_vectors[0, 0] - x_centers[-1] + x_centers[0]
    assert boundary_spacing == pytest.approx(expected_spacing, abs=1e-6)
    assert slab.basis_vectors[0, 0] == pytest.approx(chains_per_layer * expected_spacing, abs=1e-6)


def test_physical_slab_dimensions_round_up_to_complete_components() -> None:
    definition = load_allomorph("i_beta")

    slab = build_slab_by_dimensions(
        definition,
        MillerIndex(1, 1, 0),
        width=40.0,
        thickness=18.0,
        length=25.0,
        vacuum=10.0,
    )

    assert slab.basis_vectors[0, 0] >= 40.0
    assert slab.basis_vectors[2, 2] - 10.0 >= 18.0
    assert slab.basis_vectors[1, 1] >= 25.0
    assert slab.replication == (8, 4, 3)
    assert "size_mode=physical_dimensions" in slab.annotations
    assert any(
        annotation.startswith("requested_width=40.000000") for annotation in slab.annotations
    )


def test_mixed_slab_sizing_preserves_layers_and_derives_width_and_dp() -> None:
    slab = build_slab_with_sizing(
        load_allomorph("i_beta"),
        MillerIndex(1, 1, 0),
        layers=4,
        width=40.0,
        length=40.0,
        vacuum=10.0,
    )

    assert slab.replication == (8, 4, 4)
    assert len(set(slab.chain_numbers)) == 32
    assert set(slab.residue_numbers) == set(range(1, 9))
    assert "size_mode=mixed" in slab.annotations
    assert any(
        annotation.startswith("requested_width=40.000000") for annotation in slab.annotations
    )
    assert any(
        annotation.startswith("requested_length=40.000000") for annotation in slab.annotations
    )
    assert not any(annotation.startswith("requested_thickness=") for annotation in slab.annotations)


def test_i_alpha_physical_slab_length_uses_its_a_chain_axis() -> None:
    definition = load_allomorph("i_alpha")

    slab = build_slab_by_dimensions(
        definition,
        MillerIndex(0, 1, 0),
        width=20.0,
        thickness=12.0,
        length=25.0,
        vacuum=10.0,
        padding_z=12.0,
    )

    assert slab.basis_vectors[1, 1] - 12.0 >= 25.0
    assert slab.basis_vectors[1, 1] == pytest.approx(3 * 10.4 + 12.0)


@pytest.mark.parametrize(
    ("width", "thickness", "length"),
    [(0.0, 10.0, 10.0), (10.0, -1.0, 10.0), (10.0, 10.0, float("nan"))],
)
def test_physical_slab_dimensions_must_be_positive(
    width: float, thickness: float, length: float
) -> None:
    definition = load_allomorph("i_beta")

    with pytest.raises(ValueError, match="finite positive"):
        build_slab_by_dimensions(
            definition,
            MillerIndex(1, 1, 0),
            width=width,
            thickness=thickness,
            length=length,
            vacuum=0.0,
        )


def test_si_wetting_reference_slab_counts() -> None:
    definition = load_allomorph("i_beta")

    slab = build_chain_parallel_slab(
        definition,
        MillerIndex(1, 1, 0),
        layers=4,
        chains_per_layer=100,
        degree_polymerization=12,
        vacuum=0.0,
    )

    assert slab.atom_count == 100_800
    assert len(set(slab.chain_numbers)) == 400
    assert len(set(zip(slab.chain_numbers, slab.residue_numbers, strict=True))) == 4_800
    chain_numbers = np.asarray(slab.chain_numbers)
    chain_centers = [
        np.mean(slab.cartesian_coordinates[chain_numbers == chain], axis=0) for chain in (1, 2, 101)
    ]
    assert_allclose(chain_centers[1][0] - chain_centers[0][0], 5.964480, atol=1e-6)
    assert_allclose(chain_centers[2][0] - chain_centers[0][0], -0.279393, atol=1e-6)
    assert_allclose(chain_centers[2][2] - chain_centers[0][2], 5.316996, atol=1e-6)
    assert np.count_nonzero(chain_numbers == 1) == 252
