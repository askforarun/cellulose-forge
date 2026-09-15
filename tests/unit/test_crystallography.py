"""Tests for crystallographic data loading and coordinate transformations."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from cellulose_forge.crystallography import (
    MillerIndex,
    UnitCell,
    cartesian_to_fractional,
    fractional_to_cartesian,
    load_allomorph,
)


def test_i_beta_definition_matches_published_cell() -> None:
    definition = load_allomorph("i_beta")

    assert definition.name == "I-beta"
    assert definition.space_group == "P2_1"
    assert definition.doi == "10.1021/ja0257319"
    assert definition.asymmetric_unit.fractional_coordinates.shape == (42, 3)
    assert definition.asymmetric_unit.chain_labels == ("origin", "center")
    assert len(definition.symmetry_operations) == 2
    assert definition.unit_cell.volume() == pytest.approx(658.3643009)


def test_i_beta_cell_vectors_match_legacy_convention() -> None:
    cell = load_allomorph("i_beta").unit_cell

    assert_allclose(
        cell.vectors(),
        [
            [7.784, 0.0, 0.0],
            [-0.9283795561, 8.1482827884, 0.0],
            [0.0, 0.0, 10.380],
        ],
        atol=1e-9,
    )


def test_fractional_cartesian_round_trip_for_triclinic_cell() -> None:
    cell = UnitCell(6.0, 7.0, 8.0, 75.0, 82.0, 68.0)
    fractional = np.array([[0.1, 0.2, 0.3], [-0.2, 1.1, 0.7]], dtype=np.float64)

    cartesian = fractional_to_cartesian(fractional, cell)

    assert_allclose(cartesian_to_fractional(cartesian, cell), fractional, atol=1e-12)


def test_p21_operation_is_loaded_from_data() -> None:
    definition = load_allomorph("i_beta")
    coordinate = np.array([[0.25, -0.40, 0.10]], dtype=np.float64)

    transformed = definition.symmetry_operations[1].apply(coordinate)

    assert_allclose(transformed, [[-0.25, 0.40, 0.60]])


def test_i_beta_110_plane_geometry() -> None:
    cell = load_allomorph("i_beta").unit_cell
    face = MillerIndex(1, 1, 0)

    assert face.values == (1, 1, 0)
    assert face.spacing(cell) == pytest.approx(5.316996, abs=1e-6)
    assert np.linalg.norm(face.normal(cell)) == pytest.approx(1.0)


def test_miller_index_preserves_sign_and_magnitude() -> None:
    assert MillerIndex(2, 0, 0).values == (2, 0, 0)
    assert MillerIndex(-1, -1, 0).values == (-1, -1, 0)
    with pytest.raises(ValueError, match="invalid"):
        MillerIndex(0, 0, 0)


@pytest.mark.parametrize("name", ["missing", "../i_beta", "I-beta"])
def test_invalid_allomorph_identifier_is_rejected(name: str) -> None:
    with pytest.raises(ValueError):
        load_allomorph(name)
