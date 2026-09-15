"""Tests for I-beta unit-cell and legacy-compatible bulk geometry."""

from collections import Counter

from numpy.testing import assert_allclose

from cellulose_forge.crystallography import fractional_to_cartesian, load_allomorph
from cellulose_forge.geometry import build_i_beta_parallelepiped, build_unit_cell


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


def test_legacy_four_by_four_by_five_atom_count() -> None:
    definition = load_allomorph("i_beta")

    crystal = build_i_beta_parallelepiped(definition, 4, 4, 5)

    assert crystal.atom_count == 5_250
    assert crystal.atom_count == 25 * 5 * 42
