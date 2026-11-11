"""Tests for native Amber output helpers."""

import numpy as np
from numpy.testing import assert_allclose

from cellulose_forge.amber import _amber_box


def test_amber_box_preserves_triclinic_lengths_and_angles() -> None:
    basis = np.asarray(
        [
            [10.0, 0.0, 0.0],
            [2.0, 8.0, 0.0],
            [1.0, 1.5, 12.0],
        ]
    )

    box = _amber_box(basis)

    expected_lengths = np.linalg.norm(basis, axis=1)
    expected_angles = []
    for first, second in ((1, 2), (0, 2), (0, 1)):
        expected_angles.append(
            np.degrees(
                np.arccos(
                    np.dot(basis[first], basis[second])
                    / (expected_lengths[first] * expected_lengths[second])
                )
            )
        )
    assert_allclose(box, [*expected_lengths, *expected_angles])
