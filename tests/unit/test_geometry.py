"""Tests for I-beta unit-cell and legacy-compatible bulk geometry."""

from collections import Counter

from numpy.testing import assert_allclose

from cellulose_forge.crystallography import fractional_to_cartesian, load_allomorph
from cellulose_forge.geometry import build_bulk_crystal, build_unit_cell


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
