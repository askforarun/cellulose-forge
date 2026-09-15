"""Tests for I-beta unit-cell and legacy-compatible bulk geometry."""

from collections import Counter

import numpy as np
import pytest
from numpy.testing import assert_allclose

from cellulose_forge.crystallography import MillerIndex, fractional_to_cartesian, load_allomorph
from cellulose_forge.geometry import build_bulk_crystal, build_chain_parallel_slab, build_unit_cell


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
    unit_cell = build_unit_cell(definition)

    crystal = build_bulk_crystal(definition, (2, 1, 2))

    atoms_per_cell = unit_cell.atom_count
    assert_allclose(
        crystal.cartesian_coordinates[atoms_per_cell : 2 * atoms_per_cell],
        unit_cell.cartesian_coordinates + definition.unit_cell.vectors()[0],
    )
    assert_allclose(
        crystal.cartesian_coordinates[2 * atoms_per_cell : 3 * atoms_per_cell],
        unit_cell.cartesian_coordinates + definition.unit_cell.vectors()[2],
    )


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
    assert set(slab.residue_numbers) == set(range(1, 13))
    assert_allclose(slab.basis_vectors - np.diag(np.diag(slab.basis_vectors)), 0.0)
    assert slab.basis_vectors[2, 2] == pytest.approx(
        np.ptp(slab.cartesian_coordinates[:, 2]) + 30.0
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

    slab = build_chain_parallel_slab(
        definition,
        miller,
        layers=3,
        chains_per_layer=5,
        degree_polymerization=4,
        vacuum=12.0,
    )

    assert slab.atom_count == 21 * 3 * 5 * 4
    assert len(set(slab.chain_numbers)) == 15
    assert set(slab.residue_numbers) == {1, 2, 3, 4}
    assert_allclose(slab.basis_vectors - np.diag(np.diag(slab.basis_vectors)), 0.0)
    assert slab.basis_vectors[2, 2] == pytest.approx(
        np.ptp(slab.cartesian_coordinates[:, 2]) + 12.0
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
        np.mean(slab.cartesian_coordinates[chain_numbers == chain], axis=0) for chain in (0, 1, 100)
    ]
    assert_allclose(chain_centers[1][0] - chain_centers[0][0], 5.964480, atol=1e-6)
    assert_allclose(chain_centers[2][0] - chain_centers[0][0], -0.279393, atol=1e-6)
    assert_allclose(chain_centers[2][2] - chain_centers[0][2], 5.316996, atol=1e-6)
    assert np.count_nonzero(chain_numbers == 0) == 252
