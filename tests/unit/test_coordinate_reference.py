from pathlib import Path

import numpy as np
import pytest

from cellulose_forge.energy_drivers.coordinates import (
    assert_coordinate_identity,
    load_coordinate_reference,
    read_gro_frame,
    write_coordinate_reference,
)


def test_reference_round_trip_and_periodic_wrapping(tmp_path: Path) -> None:
    path = tmp_path / "coordinates.json"
    box = np.diag([10.0, 11.0, 12.0])
    digest = write_coordinate_reference(path, [[1.0, 2.0, 3.0]], box)
    reference = load_coordinate_reference(path, digest)

    assert_coordinate_identity(
        reference,
        [[11.0, 2.0, 3.0]],
        engine="test",
        box_vectors=box,
        allow_periodic_wrapping=True,
    )


def test_reference_rejects_coordinate_and_box_changes(tmp_path: Path) -> None:
    path = tmp_path / "coordinates.json"
    digest = write_coordinate_reference(path, [[1.0, 2.0, 3.0]], np.diag([10.0] * 3))
    reference = load_coordinate_reference(path, digest)

    with pytest.raises(ValueError, match="coordinates differ"):
        assert_coordinate_identity(
            reference, [[1.1, 2.0, 3.0]], engine="test", box_vectors=np.diag([10.0] * 3)
        )
    with pytest.raises(ValueError, match="box differs"):
        assert_coordinate_identity(
            reference, [[1.0, 2.0, 3.0]], engine="test", box_vectors=np.diag([11.0] * 3)
        )


def test_extended_precision_gro_reader(tmp_path: Path) -> None:
    path = tmp_path / "coordinates.gro"
    path.write_text(
        "extended\n1\n"
        "    1BGLC    C1    1   0.10000000   0.20000000   0.30000000\n"
        "  1.0000000  1.1000000  1.2000000\n",
        encoding="utf-8",
    )
    positions, box = read_gro_frame(path)
    np.testing.assert_allclose(positions, [[1.0, 2.0, 3.0]])
    np.testing.assert_allclose(box, np.diag([10.0, 11.0, 12.0]))
