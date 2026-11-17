"""Engine-neutral coordinates used only to prove cross-engine identity."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]
_BOX_VECTOR_COUNT = 3
_GRO_MINIMUM_LINES = 3
_ORTHORHOMBIC_BOX_VALUES = 3
_POSITION_DIMENSIONS = 2
_TRICLINIC_BOX_VALUES = 9


@dataclass(frozen=True)
class CoordinateReference:
    """Ordered positions and an optional periodic box, both in ångström."""

    positions: FloatArray
    box_vectors: FloatArray | None
    sha256: str


def read_gro_frame(path: Path) -> tuple[FloatArray, FloatArray]:
    """Read positions and box from a standard or extended-precision GRO file."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < _GRO_MINIMUM_LINES:
        raise ValueError("GRO coordinate file is truncated")
    try:
        count = int(lines[1].strip())
    except ValueError as error:
        raise ValueError("GRO atom count is invalid") from error
    if len(lines) != count + 3:
        raise ValueError("GRO atom count does not match the coordinate rows")
    try:
        positions = []
        for line in lines[2 : 2 + count]:
            first_decimal = line.index(".", 20)
            second_decimal = line.index(".", first_decimal + 1)
            width = second_decimal - first_decimal
            positions.append(
                [
                    10.0 * float(line[20 + axis * width : 20 + (axis + 1) * width])
                    for axis in range(3)
                ]
            )
        values = [10.0 * float(value) for value in lines[-1].split()]
    except (ValueError, IndexError) as error:
        raise ValueError("GRO coordinates are invalid") from error
    if len(values) == _ORTHORHOMBIC_BOX_VALUES:
        box = np.diag(values)
    elif len(values) == _TRICLINIC_BOX_VALUES:
        box = np.asarray(
            [
                [values[0], values[3], values[4]],
                [values[5], values[1], values[6]],
                [values[7], values[8], values[2]],
            ]
        )
    else:
        raise ValueError("GRO box must contain three or nine values")
    return np.asarray(positions), box


def load_coordinate_reference(path: Path, expected_sha256: str) -> CoordinateReference:
    """Load a strict, engine-neutral JSON coordinate reference."""
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"coordinate reference does not exist: {path}")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != expected_sha256:
        raise ValueError("coordinate reference does not match the manifest checksum")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise ValueError("coordinate reference must use schema version 1")
    if payload.get("units") != "angstrom":
        raise ValueError("coordinate reference units must be angstrom")
    try:
        positions = np.asarray(payload["positions"], dtype=np.float64)
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("coordinate reference positions are invalid") from error
    if (
        positions.ndim != _POSITION_DIMENSIONS
        or positions.shape[1:] != (_BOX_VECTOR_COUNT,)
        or not np.isfinite(positions).all()
    ):
        raise ValueError("coordinate reference positions must be a finite N x 3 array")
    raw_box = payload.get("box_vectors")
    box = None if raw_box is None else np.asarray(raw_box, dtype=np.float64)
    if box is not None and (box.shape != (3, 3) or not np.isfinite(box).all()):
        raise ValueError("coordinate reference box_vectors must be a finite 3 x 3 array")
    return CoordinateReference(positions, box, digest)


def assert_coordinate_identity(  # noqa: PLR0913 - explicit scientific contract
    reference: CoordinateReference,
    positions: npt.ArrayLike,
    *,
    engine: str,
    box_vectors: npt.ArrayLike | None,
    allow_periodic_wrapping: bool = False,
    tolerance_angstrom: float = 1.0e-8,
) -> None:
    """Require native coordinates and the native periodic box to match the reference."""
    native = np.asarray(positions, dtype=np.float64)
    if native.shape != reference.positions.shape or not np.isfinite(native).all():
        raise ValueError(f"{engine} coordinates differ from the coordinate reference")
    delta = native - reference.positions
    if allow_periodic_wrapping:
        if reference.box_vectors is None:
            raise ValueError(f"{engine} periodic wrapping requires a reference box")
        fractional = delta @ np.linalg.inv(reference.box_vectors)
        delta -= np.rint(fractional) @ reference.box_vectors
    if not np.allclose(delta, 0.0, rtol=0.0, atol=tolerance_angstrom):
        raise ValueError(f"{engine} coordinates differ from the coordinate reference")
    if reference.box_vectors is None:
        if box_vectors is not None:
            raise ValueError(f"{engine} has a periodic box but the coordinate reference does not")
        return
    if box_vectors is None:
        raise ValueError(f"{engine} is missing the coordinate-reference box")
    native_box = np.asarray(box_vectors, dtype=np.float64)
    if native_box.shape != (3, 3) or not np.allclose(
        native_box, reference.box_vectors, rtol=0.0, atol=tolerance_angstrom
    ):
        raise ValueError(f"{engine} box differs from the coordinate reference")


def write_coordinate_reference(
    path: Path, positions: npt.ArrayLike, box_vectors: npt.ArrayLike | None
) -> str:
    """Write a deterministic reference JSON and return its SHA-256 digest."""
    payload = {
        "schema_version": 1,
        "units": "angstrom",
        "positions": np.asarray(positions, dtype=np.float64).tolist(),
        "box_vectors": None
        if box_vectors is None
        else np.asarray(box_vectors, dtype=np.float64).tolist(),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return hashlib.sha256(path.read_bytes()).hexdigest()
