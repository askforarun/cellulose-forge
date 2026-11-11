"""Validated crystallographic definitions and coordinate transformations."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from importlib import resources
from math import acos, cos, degrees, radians, sin, sqrt
from typing import Any, cast

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]
_ARRAY_DIMENSIONS = 2
_COORDINATE_DIMENSIONS = 3
_STRAIGHT_ANGLE_DEGREES = 180.0
_SINGULAR_TOLERANCE = 1e-12


def _readonly_array(value: object, shape: tuple[int, ...], name: str) -> FloatArray:
    array = np.asarray(value, dtype=np.float64)
    if array.shape != shape:
        msg = f"{name} must have shape {shape}, got {array.shape}"
        raise ValueError(msg)
    if not np.all(np.isfinite(array)):
        msg = f"{name} must contain only finite values"
        raise ValueError(msg)
    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class UnitCell:
    """Triclinic unit-cell lengths (angstrom) and angles (degrees)."""

    a: float
    b: float
    c: float
    alpha_deg: float
    beta_deg: float
    gamma_deg: float

    def __post_init__(self) -> None:
        if not all(np.isfinite(value) and value > 0 for value in (self.a, self.b, self.c)):
            msg = "unit-cell lengths must be finite and greater than zero"
            raise ValueError(msg)
        if not all(
            np.isfinite(value) and 0 < value < _STRAIGHT_ANGLE_DEGREES
            for value in (self.alpha_deg, self.beta_deg, self.gamma_deg)
        ):
            msg = "unit-cell angles must be finite and between 0 and 180 degrees"
            raise ValueError(msg)
        self.vectors()

    @property
    def gamma_rad(self) -> float:
        """Return gamma in radians for compatibility with legacy geometry code."""
        return radians(self.gamma_deg)

    def vectors(self) -> FloatArray:
        """Return row-wise Cartesian cell vectors for a general triclinic cell."""
        alpha = radians(self.alpha_deg)
        beta = radians(self.beta_deg)
        gamma = radians(self.gamma_deg)
        sin_gamma = sin(gamma)
        if abs(sin_gamma) < _SINGULAR_TOLERANCE:
            msg = "gamma produces a singular unit cell"
            raise ValueError(msg)

        vector_a = np.array([self.a, 0.0, 0.0], dtype=np.float64)
        vector_b = np.array([self.b * cos(gamma), self.b * sin_gamma, 0.0], dtype=np.float64)
        c_x = self.c * cos(beta)
        c_y = self.c * (cos(alpha) - cos(beta) * cos(gamma)) / sin_gamma
        c_z_squared = self.c**2 - c_x**2 - c_y**2
        if c_z_squared <= 0:
            msg = "unit-cell angles do not define a positive-volume cell"
            raise ValueError(msg)
        vector_c = np.array([c_x, c_y, sqrt(c_z_squared)], dtype=np.float64)
        vectors = np.vstack([vector_a, vector_b, vector_c])
        vectors.setflags(write=False)
        return vectors

    def volume(self) -> float:
        """Return the unit-cell volume in cubic angstrom."""
        return float(abs(np.linalg.det(self.vectors())))


@dataclass(frozen=True)
class SpaceGroupOperation:
    """A crystallographic transformation ``x' = R x + t``."""

    name: str
    rotation: FloatArray
    translation: FloatArray

    def __post_init__(self) -> None:
        object.__setattr__(self, "rotation", _readonly_array(self.rotation, (3, 3), "rotation"))
        object.__setattr__(
            self, "translation", _readonly_array(self.translation, (3,), "translation")
        )

    def apply(self, fractional: FloatArray) -> FloatArray:
        """Apply the operation to an ``(n, 3)`` fractional-coordinate array."""
        coordinates = np.asarray(fractional, dtype=np.float64)
        if coordinates.ndim != _ARRAY_DIMENSIONS or coordinates.shape[1] != _COORDINATE_DIMENSIONS:
            msg = f"fractional coordinates must have shape (n, 3), got {coordinates.shape}"
            raise ValueError(msg)
        return coordinates @ self.rotation.T + self.translation


@dataclass(frozen=True)
class MillerIndex:
    """An exact, signed crystallographic plane index."""

    h: int
    k: int
    l: int  # noqa: E741 - conventional third Miller index

    def __post_init__(self) -> None:
        values = (self.h, self.k, self.l)
        if any(not isinstance(value, int) or isinstance(value, bool) for value in values):
            msg = "Miller indices must be integers"
            raise TypeError(msg)
        if values == (0, 0, 0):
            msg = "Miller index (0, 0, 0) is invalid"
            raise ValueError(msg)

    @property
    def values(self) -> tuple[int, int, int]:
        """Return ``(h, k, l)`` without reducing magnitude or changing sign."""
        return self.h, self.k, self.l

    def reciprocal_vector(self, unit_cell: UnitCell) -> FloatArray:
        """Return the Cartesian reciprocal vector without the ``2 pi`` factor."""
        return np.asarray(
            np.linalg.solve(unit_cell.vectors(), np.asarray(self.values, dtype=np.float64)),
            dtype=np.float64,
        )

    def normal(self, unit_cell: UnitCell) -> FloatArray:
        """Return the oriented Cartesian unit normal."""
        reciprocal = self.reciprocal_vector(unit_cell)
        return reciprocal / np.linalg.norm(reciprocal)

    def spacing(self, unit_cell: UnitCell) -> float:
        """Return the interplanar spacing in angstrom."""
        return float(1.0 / np.linalg.norm(self.reciprocal_vector(unit_cell)))


@dataclass(frozen=True)
class AsymmetricUnit:
    """Ordered atom labels and fractional coordinates for one asymmetric unit."""

    atom_labels: tuple[str, ...]
    fractional_coordinates: FloatArray
    residue_atom_count: int
    chain_labels: tuple[str, ...]
    chain_directions: tuple[int, ...]

    def __post_init__(self) -> None:
        coordinates = np.asarray(self.fractional_coordinates, dtype=np.float64)
        if coordinates.ndim != _ARRAY_DIMENSIONS or coordinates.shape[1] != _COORDINATE_DIMENSIONS:
            msg = f"fractional_coordinates must have shape (n, 3), got {coordinates.shape}"
            raise ValueError(msg)
        if len(self.atom_labels) != len(coordinates):
            msg = (
                f"atom_labels has {len(self.atom_labels)} entries but "
                f"fractional_coordinates has {len(coordinates)} rows"
            )
            raise ValueError(msg)
        if not self.atom_labels or any(not label for label in self.atom_labels):
            msg = "atom labels must be nonempty"
            raise ValueError(msg)
        if not np.all(np.isfinite(coordinates)):
            msg = "fractional coordinates must contain only finite values"
            raise ValueError(msg)
        if self.residue_atom_count <= 0 or len(coordinates) % self.residue_atom_count:
            msg = "residue_atom_count must divide the asymmetric-unit atom count"
            raise ValueError(msg)
        residue_count = len(coordinates) // self.residue_atom_count
        if len(self.chain_labels) != residue_count:
            msg = f"chain_labels must contain {residue_count} entries"
            raise ValueError(msg)
        if len(self.chain_directions) != residue_count or any(
            direction not in (-1, 1) for direction in self.chain_directions
        ):
            msg = f"chain_directions must contain {residue_count} entries equal to -1 or 1"
            raise ValueError(msg)
        directions_by_chain: dict[str, int] = {}
        for chain, direction in zip(self.chain_labels, self.chain_directions, strict=True):
            if chain in directions_by_chain and directions_by_chain[chain] != direction:
                msg = f"chain direction is inconsistent for chain label {chain!r}"
                raise ValueError(msg)
            directions_by_chain[chain] = direction
        coordinates.setflags(write=False)
        object.__setattr__(self, "fractional_coordinates", coordinates)


@dataclass(frozen=True)
class AllomorphDefinition:
    """A validated, traceable crystallographic definition of a cellulose allomorph."""

    name: str
    space_group: str
    chain_axis: int
    unit_cell: UnitCell
    asymmetric_unit: AsymmetricUnit
    symmetry_operations: tuple[SpaceGroupOperation, ...]
    citation: str
    doi: str
    polymer: str = "cellulose"
    residue_name: str = "BGC"
    temperature_k: float | None = None

    def __post_init__(self) -> None:
        if not self.name or not self.space_group or not self.symmetry_operations:
            msg = "allomorph name, space group, and symmetry operations are required"
            raise ValueError(msg)
        if self.chain_axis not in (0, 1, 2):
            msg = "chain_axis must identify crystallographic a, b, or c"
            raise ValueError(msg)


SUPPORTED_POLYMERS = ("cellulose", "chitin")


def _load_payload(name: str, polymer: str = "cellulose") -> dict[str, Any]:
    if re.fullmatch(r"[a-z][a-z0-9_]*", name) is None:
        msg = f"invalid allomorph identifier: {name!r}"
        raise ValueError(msg)
    if polymer not in SUPPORTED_POLYMERS:
        msg = f"unknown polymer: {polymer!r}"
        raise ValueError(msg)
    directory = resources.files("cellulose_forge.data.crystallography")
    if polymer != "cellulose":
        directory = directory.joinpath(polymer)
    resource = directory.joinpath(f"{name}.json")
    if not resource.is_file():
        msg = f"unknown allomorph definition: {name!r}"
        raise ValueError(msg)
    with resource.open(encoding="utf-8") as handle:
        return cast(dict[str, Any], json.load(handle))


def load_allomorph(name: str, polymer: str = "cellulose") -> AllomorphDefinition:
    """Load a validated, package-owned crystallographic definition."""
    payload = _load_payload(name, polymer)
    cell_payload = cast(dict[str, Any], payload["unit_cell"])
    unit_cell = UnitCell(
        a=float(cell_payload["a"]),
        b=float(cell_payload["b"]),
        c=float(cell_payload["c"]),
        alpha_deg=float(cell_payload["alpha_deg"]),
        beta_deg=float(cell_payload["beta_deg"]),
        gamma_deg=float(cell_payload["gamma_deg"]),
    )

    asy_payload = cast(dict[str, Any], payload["asymmetric_unit"])
    asymmetric_unit = AsymmetricUnit(
        atom_labels=tuple(cast(list[str], asy_payload["atom_labels"])),
        fractional_coordinates=np.asarray(asy_payload["fractional_coordinates"], dtype=np.float64),
        residue_atom_count=int(asy_payload["residue_atom_count"]),
        chain_labels=tuple(cast(list[str], asy_payload["chain_labels"])),
        chain_directions=tuple(cast(list[int], asy_payload["chain_directions"])),
    )

    operation_payloads = cast(list[dict[str, Any]], payload["symmetry_operations"])
    symmetry_operations = tuple(
        SpaceGroupOperation(
            name=str(operation["name"]),
            rotation=np.asarray(operation["rotation"], dtype=np.float64),
            translation=np.asarray(operation["translation"], dtype=np.float64),
        )
        for operation in operation_payloads
    )
    provenance = cast(dict[str, Any], payload["provenance"])

    return AllomorphDefinition(
        name=str(payload["allomorph"]),
        space_group=str(payload["space_group"]),
        chain_axis={"a": 0, "b": 1, "c": 2}[str(payload["chain_axis"])],
        unit_cell=unit_cell,
        asymmetric_unit=asymmetric_unit,
        symmetry_operations=symmetry_operations,
        citation=str(provenance["citation"]),
        doi=str(provenance["doi"]),
        polymer=str(payload.get("polymer", "cellulose")),
        residue_name=str(payload.get("residue_name", "BGC")),
        temperature_k=(
            float(provenance["temperature_k"]) if "temperature_k" in provenance else None
        ),
    )


def fractional_to_cartesian(fractional: FloatArray, unit_cell: UnitCell) -> FloatArray:
    """Convert an ``(n, 3)`` fractional-coordinate array to Cartesian coordinates."""
    coordinates = np.asarray(fractional, dtype=np.float64)
    if coordinates.ndim != _ARRAY_DIMENSIONS or coordinates.shape[1] != _COORDINATE_DIMENSIONS:
        msg = f"fractional coordinates must have shape (n, 3), got {coordinates.shape}"
        raise ValueError(msg)
    return coordinates @ unit_cell.vectors()


def cartesian_to_fractional(cartesian: FloatArray, unit_cell: UnitCell) -> FloatArray:
    """Convert an ``(n, 3)`` Cartesian-coordinate array to fractional coordinates."""
    coordinates = np.asarray(cartesian, dtype=np.float64)
    if coordinates.ndim != _ARRAY_DIMENSIONS or coordinates.shape[1] != _COORDINATE_DIMENSIONS:
        msg = f"Cartesian coordinates must have shape (n, 3), got {coordinates.shape}"
        raise ValueError(msg)
    result = np.linalg.solve(unit_cell.vectors().T, coordinates.T).T
    return np.asarray(result, dtype=np.float64)


def basis_lengths_and_angles(basis_vectors: FloatArray) -> tuple[float, ...]:
    """Return ``a, b, c, alpha, beta, gamma`` for three row-wise basis vectors."""
    vectors = _readonly_array(basis_vectors, (3, 3), "basis_vectors")
    lengths = np.linalg.norm(vectors, axis=1)

    def angle(first: FloatArray, second: FloatArray, denominator: float) -> float:
        cosine = float(np.dot(first, second) / denominator)
        return degrees(acos(float(np.clip(cosine, -1.0, 1.0))))

    alpha = angle(vectors[1], vectors[2], float(lengths[1] * lengths[2]))
    beta = angle(vectors[0], vectors[2], float(lengths[0] * lengths[2]))
    gamma = angle(vectors[0], vectors[1], float(lengths[0] * lengths[1]))
    return float(lengths[0]), float(lengths[1]), float(lengths[2]), alpha, beta, gamma
