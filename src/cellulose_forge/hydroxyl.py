"""Geometry helpers for hydroxyl hydrogens absent from diffraction data."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

from cellulose_forge.crystallography import FloatArray

_DEGENERACY_TOLERANCE = 1.0e-12
_PARALLEL_AXIS_THRESHOLD = 0.9
_COORDINATE_MATCH_TOLERANCE = 1.0e-6
_TORSION_SAMPLES = 360
_ENVIRONMENT_RADIUS_ANGSTROM = 4.5
# Minimum separations accepted for a modelled hydroxyl hydrogen.  Oxygen is
# allowed closer than carbon because it is the hydrogen-bond acceptor.
MINIMUM_HYDROGEN_SEPARATION_ANGSTROM = 2.0
MINIMUM_OXYGEN_SEPARATION_ANGSTROM = 1.6
MINIMUM_OTHER_HEAVY_SEPARATION_ANGSTROM = 2.0
# Geometric hydrogen-bond criteria used to choose among clash-free torsions.
MAXIMUM_DONOR_ACCEPTOR_ANGSTROM = 3.3
MAXIMUM_HYDROGEN_ACCEPTOR_ANGSTROM = 2.6
# H...O distance at or below which a hydrogen bond scores full closeness.
IDEAL_HYDROGEN_ACCEPTOR_ANGSTROM = 1.8
# Weight of reference-torsion alignment (range -1..1) relative to bond quality.
_ALIGNMENT_TIE_BREAK = 0.01
# Hydroxyl classes are placed in this order; see HydroxylHydrogenPlacer.resolve.
PLACEMENT_ORDER = ("O3", "O2", "O6")
MINIMUM_HYDROGEN_BOND_ANGLE_DEGREES = 120.0


def _cone_candidates(
    carbon: FloatArray,
    oxygen: FloatArray,
    preferred_hydrogen: FloatArray,
    bond_length: float,
    angle_degrees: float,
) -> tuple[FloatArray, FloatArray]:
    """Return sampled hydrogen positions on the exact C--O--H cone.

    The torsion grid starts from a fixed global direction, not from the
    preferred hydrogen, so every builder evaluates identical trial positions;
    the preferred hydrogen only affects the returned alignment scores.
    """
    carbon_direction = np.asarray(carbon, dtype=np.float64) - oxygen
    carbon_direction /= np.linalg.norm(carbon_direction)
    preferred = np.asarray(preferred_hydrogen, dtype=np.float64) - oxygen
    trial = np.asarray((1.0, 0.0, 0.0))
    if abs(float(np.dot(trial, carbon_direction))) > _PARALLEL_AXIS_THRESHOLD:
        trial = np.asarray((0.0, 1.0, 0.0))
    perpendicular = trial - np.dot(trial, carbon_direction) * carbon_direction
    perpendicular /= np.linalg.norm(perpendicular)
    cross_direction = np.cross(carbon_direction, perpendicular)
    torsions = np.linspace(0.0, 2.0 * np.pi, _TORSION_SAMPLES, endpoint=False)
    cone = np.cos(torsions)[:, None] * perpendicular + np.sin(torsions)[:, None] * cross_direction
    angle = np.radians(angle_degrees)
    directions = np.cos(angle) * carbon_direction + np.sin(angle) * cone
    candidates = np.asarray(oxygen, dtype=np.float64) + bond_length * directions
    preferred_norm = float(np.linalg.norm(preferred))
    alignment = (
        directions @ (preferred / preferred_norm)
        if preferred_norm > _DEGENERACY_TOLERANCE
        else np.zeros(len(directions))
    )
    return candidates, alignment


class HydroxylHydrogenPlacer:
    """Place missing hydroxyl hydrogens on a tetrahedral cone, one at a time.

    Among clash-free torsions, the placer chooses the one with the highest
    geometric hydrogen-bond quality, which rewards both a short H...O distance
    and a near-linear O-H...O angle; the reference torsion breaks near-ties and
    is used whenever no torsion forms a hydrogen bond.  Each placed hydrogen
    joins the environment so later hydrogens avoid it.

    Neighbours are found with a periodic cell list whose bins are at least the
    environment radius wide perpendicular to every lattice face, so each
    placement examines a bounded neighbourhood and the total cost grows
    linearly with system size.
    """

    def __init__(
        self,
        labels: Sequence[str],
        coordinates: FloatArray,
        periodic_box: FloatArray | None = None,
        *,
        bond_length: float = 0.96,
        angle_degrees: float = 108.5,
    ) -> None:
        self._coordinates = np.asarray(coordinates, dtype=np.float64)
        self._labels = np.asarray([str(label).upper() for label in labels])
        elements = np.asarray([label[:1] for label in self._labels])
        self._hydrogen = elements == "H"
        self._oxygen = elements == "O"
        self._box = None if periodic_box is None else np.asarray(periodic_box, dtype=np.float64)
        self._inverse_box = None if self._box is None else np.linalg.inv(self._box)
        self._bond_length = bond_length
        self._angle_degrees = angle_degrees
        self._shape = self._bin_shape()
        self._crystal_bins: dict[tuple[int, int, int], list[int]] = {}
        for index, key in enumerate(self._bin_keys(self._coordinates)):
            self._crystal_bins.setdefault(key, []).append(index)
        self._placed_bins: dict[tuple[int, int, int], list[FloatArray]] = {}
        self._requests: list[
            tuple[int, str, FloatArray, FloatArray, FloatArray, str | None, FloatArray | None]
        ] = []

    def _bin_shape(self) -> tuple[int, int, int]:
        if self._box is None:
            return (1, 1, 1)
        volume = abs(float(np.linalg.det(self._box)))
        heights = [
            volume / float(np.linalg.norm(np.cross(self._box[(i + 1) % 3], self._box[(i + 2) % 3])))
            for i in range(3)
        ]
        counts = [max(1, int(height // _ENVIRONMENT_RADIUS_ANGSTROM)) for height in heights]
        return (counts[0], counts[1], counts[2])

    def _bin_keys(self, points: FloatArray) -> list[tuple[int, int, int]]:
        if self._inverse_box is None:
            return [(0, 0, 0)] * len(points)
        fractional = np.asarray(points, dtype=np.float64).reshape(-1, 3) @ self._inverse_box
        fractional -= np.floor(fractional)
        indices = np.floor(fractional * np.asarray(self._shape)).astype(int) % np.asarray(
            self._shape
        )
        return [(int(i), int(j), int(k)) for i, j, k in indices]

    def _neighbour_keys(self, point: FloatArray) -> list[tuple[int, int, int]]:
        centre = self._bin_keys(point)[0]
        axes = [
            sorted({(centre[axis] + offset) % count for offset in (-1, 0, 1)})
            if count >= 3  # noqa: PLR2004
            else list(range(count))
            for axis, count in enumerate(self._shape)
        ]
        return [(i, j, k) for i in axes[0] for j in axes[1] for k in axes[2]]

    def _minimum_image(self, delta: FloatArray) -> FloatArray:
        if self._box is None or self._inverse_box is None:
            return delta
        fractional = delta @ self._inverse_box
        fractional -= np.round(fractional)
        return np.asarray(fractional @ self._box, dtype=np.float64)

    def place(  # noqa: PLR0913, PLR0915 - geometric selection kept in one place
        self,
        carbon: FloatArray,
        oxygen: FloatArray,
        preferred_hydrogen: FloatArray,
        *,
        preferred_acceptor: str | None = None,
        preferred_acceptor_coordinate: FloatArray | None = None,
        record: bool = True,
    ) -> FloatArray:
        """Return one hydrogen position and add it to the placed environment.

        ``preferred_acceptor`` names an acceptor oxygen label whose hydrogen
        bond, when any clash-free torsion can form one, takes priority over
        the purely geometric choice.  Cellulose builders pass ``"O5"`` for
        HO3, encoding the intrachain O3-H...O5 bond present in every cellulose
        allomorph. ``preferred_acceptor_coordinate`` identifies one exact
        periodic acceptor site and takes precedence over a label. It is used
        for the chain-specific cellulose-II network determined by neutron
        diffraction.
        """
        oxygen = np.asarray(oxygen, dtype=np.float64)
        candidates, alignment = _cone_candidates(
            carbon, oxygen, preferred_hydrogen, self._bond_length, self._angle_degrees
        )
        keys = self._neighbour_keys(oxygen)
        indices = np.asarray(
            [index for key in keys for index in self._crystal_bins.get(key, [])], dtype=int
        )
        placed = [point for key in keys for point in self._placed_bins.get(key, [])]
        crystal_delta = self._minimum_image(self._coordinates[indices] - oxygen)
        carbon_delta = self._minimum_image(self._coordinates[indices] - np.asarray(carbon))
        crystal_distance = np.linalg.norm(crystal_delta, axis=1)
        near = (
            (crystal_distance < _ENVIRONMENT_RADIUS_ANGSTROM)
            & (crystal_distance > _COORDINATE_MATCH_TOLERANCE)
            & (np.linalg.norm(carbon_delta, axis=1) > _COORDINATE_MATCH_TOLERANCE)
        )
        placed_delta = self._minimum_image(
            np.asarray(placed, dtype=np.float64).reshape(-1, 3) - oxygen
        )
        placed_near = np.linalg.norm(placed_delta, axis=1) < _ENVIRONMENT_RADIUS_ANGSTROM
        environment = oxygen + np.vstack((crystal_delta[near], placed_delta[placed_near]))
        is_hydrogen = np.concatenate(
            (self._hydrogen[indices][near], np.ones(int(placed_near.sum()), bool))
        )
        is_oxygen = np.concatenate(
            (self._oxygen[indices][near], np.zeros(int(placed_near.sum()), bool))
        )
        labels = np.concatenate((self._labels[indices][near], np.full(int(placed_near.sum()), "H")))
        distances = np.linalg.norm(candidates[:, None, :] - environment[None, :, :], axis=2)

        def closest(mask: np.ndarray[Any, np.dtype[np.bool_]]) -> FloatArray:
            if not mask.any():
                return np.full(len(candidates), np.inf)
            return np.asarray(distances[:, mask].min(axis=1))

        hydrogen_clearance = closest(is_hydrogen)
        oxygen_clearance = closest(is_oxygen)
        other_clearance = closest(~is_hydrogen & ~is_oxygen)
        feasible = (
            (hydrogen_clearance >= MINIMUM_HYDROGEN_SEPARATION_ANGSTROM)
            & (oxygen_clearance >= MINIMUM_OXYGEN_SEPARATION_ANGSTROM)
            & (other_clearance >= MINIMUM_OTHER_HEAVY_SEPARATION_ANGSTROM)
        )

        acceptors = is_oxygen & (
            np.linalg.norm(environment - oxygen, axis=1) <= MAXIMUM_DONOR_ACCEPTOR_ANGSTROM
        )
        quality = np.zeros(len(candidates))
        preferred_quality = np.zeros(len(candidates))
        if acceptors.any():
            to_donor = oxygen - candidates
            to_acceptor = environment[acceptors][None, :, :] - candidates[:, None, :]
            acceptor_distance = distances[:, acceptors]
            cosine = np.einsum("ij,ikj->ik", to_donor, to_acceptor) / (
                np.linalg.norm(to_donor, axis=1)[:, None] * acceptor_distance
            )
            angle = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
            # Quality is 1 for a linear (180 degree) bond at or below the ideal
            # distance and falls to 0 at the 120 degree and 2.6 A limits, so a
            # short but strongly bent contact cannot beat a near-linear one.
            linearity = np.clip(
                (angle - MINIMUM_HYDROGEN_BOND_ANGLE_DEGREES)
                / (180.0 - MINIMUM_HYDROGEN_BOND_ANGLE_DEGREES),
                0.0,
                1.0,
            )
            closeness = np.clip(
                (MAXIMUM_HYDROGEN_ACCEPTOR_ANGSTROM - acceptor_distance)
                / (MAXIMUM_HYDROGEN_ACCEPTOR_ANGSTROM - IDEAL_HYDROGEN_ACCEPTOR_ANGSTROM),
                0.0,
                1.0,
            )
            per_acceptor = linearity * closeness
            quality = per_acceptor.max(axis=1)
            if preferred_acceptor is not None:
                preferred = labels[acceptors] == preferred_acceptor.upper()
                if preferred.any():
                    preferred_quality = per_acceptor[:, preferred].max(axis=1)
        if preferred_acceptor_coordinate is not None:
            preferred_delta = self._minimum_image(
                np.asarray(preferred_acceptor_coordinate, dtype=np.float64) - oxygen
            )
            preferred_position = oxygen + preferred_delta
            to_donor = oxygen - candidates
            to_acceptor = preferred_position - candidates
            acceptor_distance = np.linalg.norm(to_acceptor, axis=1)
            cosine = np.einsum("ij,ij->i", to_donor, to_acceptor) / (
                np.linalg.norm(to_donor, axis=1) * acceptor_distance
            )
            angle = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
            linearity = np.clip(
                (angle - MINIMUM_HYDROGEN_BOND_ANGLE_DEGREES)
                / (180.0 - MINIMUM_HYDROGEN_BOND_ANGLE_DEGREES),
                0.0,
                1.0,
            )
            closeness = np.clip(
                (MAXIMUM_HYDROGEN_ACCEPTOR_ANGSTROM - acceptor_distance)
                / (MAXIMUM_HYDROGEN_ACCEPTOR_ANGSTROM - IDEAL_HYDROGEN_ACCEPTOR_ANGSTROM),
                0.0,
                1.0,
            )
            preferred_quality = linearity * closeness
        bonded = feasible & (quality > 0.0)
        preferred_bonded = feasible & (preferred_quality > 0.0)

        if preferred_bonded.any():
            score = np.where(
                preferred_bonded,
                preferred_quality + _ALIGNMENT_TIE_BREAK * alignment,
                -np.inf,
            )
            choice = int(np.argmax(score))
        elif bonded.any():
            # Alignment with the reference torsion only breaks near-ties.
            score = np.where(bonded, quality + _ALIGNMENT_TIE_BREAK * alignment, -np.inf)
            choice = int(np.argmax(score))
        elif feasible.any():
            choice = int(np.argmax(np.where(feasible, alignment, -np.inf)))
        else:
            clearance = np.minimum.reduce(
                (
                    hydrogen_clearance / MINIMUM_HYDROGEN_SEPARATION_ANGSTROM,
                    oxygen_clearance / MINIMUM_OXYGEN_SEPARATION_ANGSTROM,
                    other_clearance / MINIMUM_OTHER_HEAVY_SEPARATION_ANGSTROM,
                )
            )
            choice = int(np.argmax(clearance))
        coordinate = np.asarray(candidates[choice], dtype=np.float64)
        if record:
            self._record(coordinate)
        coordinate.setflags(write=False)
        return coordinate

    def _record(self, coordinate: FloatArray) -> None:
        self._placed_bins.setdefault(self._bin_keys(coordinate)[0], []).append(
            np.array(coordinate, dtype=np.float64)
        )

    def request(  # noqa: PLR0913
        self,
        key: int,
        carbon: FloatArray,
        oxygen: FloatArray,
        preferred_hydrogen: FloatArray,
        *,
        donor: str,
        preferred_acceptor: str | None = None,
        preferred_acceptor_coordinate: FloatArray | None = None,
    ) -> None:
        """Queue one hydroxyl hydrogen for order-independent placement."""
        self._requests.append(
            (
                key,
                donor.upper(),
                carbon,
                oxygen,
                preferred_hydrogen,
                preferred_acceptor,
                preferred_acceptor_coordinate,
            )
        )

    def resolve(self) -> dict[int, FloatArray]:
        """Place every queued hydrogen class by class and return them by key.

        Classes are placed in ``PLACEMENT_ORDER`` (O3, then O2, then O6). A
        hydrogen sees the crystal and every hydrogen of earlier classes, but
        not other members of its own class, so symmetry-equivalent hydroxyls
        always see identical environments and the result does not depend on
        the order in which a builder visits atoms.
        """
        rank = {name: index for index, name in enumerate(PLACEMENT_ORDER)}
        donors = sorted(
            {request[1] for request in self._requests}, key=lambda d: (rank.get(d, len(rank)), d)
        )
        results: dict[int, FloatArray] = {}
        for donor in donors:
            batch = [
                (
                    request[0],
                    self.place(
                        request[2],
                        request[3],
                        request[4],
                        preferred_acceptor=request[5],
                        preferred_acceptor_coordinate=request[6],
                        record=False,
                    ),
                )
                for request in self._requests
                if request[1] == donor
            ]
            for key, coordinate in batch:
                self._record(coordinate)
                results[key] = coordinate
        self._requests.clear()
        return results
