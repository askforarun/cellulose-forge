"""Geometry helpers for hydroxyl hydrogens absent from diffraction data."""

from __future__ import annotations

import numpy as np

from cellulose_forge.crystallography import FloatArray

_DEGENERACY_TOLERANCE = 1.0e-12
_PARALLEL_AXIS_THRESHOLD = 0.9
_COORDINATE_MATCH_TOLERANCE = 1.0e-6


def tetrahedral_hydroxyl_hydrogen(  # noqa: PLR0913
    carbon: FloatArray,
    oxygen: FloatArray,
    preferred_hydrogen: FloatArray,
    *,
    bond_length: float = 0.96,
    angle_degrees: float = 108.5,
    avoid_coordinates: FloatArray | None = None,
    periodic_box: FloatArray | None = None,
) -> FloatArray:
    """Project a preferred O--H torsion onto an exact C--O--H cone."""
    carbon_direction = np.asarray(carbon, dtype=np.float64) - oxygen
    carbon_direction /= np.linalg.norm(carbon_direction)
    preferred = np.asarray(preferred_hydrogen, dtype=np.float64) - oxygen
    perpendicular = preferred - np.dot(preferred, carbon_direction) * carbon_direction
    norm = float(np.linalg.norm(perpendicular))
    if norm <= _DEGENERACY_TOLERANCE:
        trial = np.asarray((1.0, 0.0, 0.0))
        if abs(float(np.dot(trial, carbon_direction))) > _PARALLEL_AXIS_THRESHOLD:
            trial = np.asarray((0.0, 1.0, 0.0))
        perpendicular = trial - np.dot(trial, carbon_direction) * carbon_direction
        norm = float(np.linalg.norm(perpendicular))
    perpendicular /= norm
    angle = np.radians(angle_degrees)
    cross_direction = np.cross(carbon_direction, perpendicular)
    candidates: list[FloatArray] = []
    for torsion in np.linspace(0.0, 2.0 * np.pi, 360, endpoint=False):
        cone = np.cos(torsion) * perpendicular + np.sin(torsion) * cross_direction
        direction = np.cos(angle) * carbon_direction + np.sin(angle) * cone
        candidates.append(np.asarray(oxygen, dtype=np.float64) + bond_length * direction)
    if avoid_coordinates is None:
        coordinate = candidates[0]
    else:
        avoid = np.asarray(avoid_coordinates, dtype=np.float64)
        keep = np.logical_and(
            np.linalg.norm(avoid - oxygen, axis=1) > _COORDINATE_MATCH_TOLERANCE,
            np.linalg.norm(avoid - carbon, axis=1) > _COORDINATE_MATCH_TOLERANCE,
        )
        avoid = avoid[keep]
        inverse_box = np.linalg.inv(periodic_box) if periodic_box is not None else None

        def score(candidate: FloatArray) -> tuple[float, float]:
            delta = avoid - candidate
            if periodic_box is not None and inverse_box is not None:
                fractional = delta @ inverse_box
                fractional -= np.round(fractional)
                delta = fractional @ periodic_box
            closest = float(np.min(np.linalg.norm(delta, axis=1)))
            preferred_alignment = float(
                np.dot((candidate - oxygen) / bond_length, preferred / np.linalg.norm(preferred))
            )
            return (closest, preferred_alignment)

        coordinate = max(candidates, key=score)
    coordinate.setflags(write=False)
    return coordinate
