"""Modelled hydroxyl-hydrogen placement: order independence and scaling."""

from __future__ import annotations

import time

import numpy as np

from cellulose_forge.crystallography import load_allomorph
from cellulose_forge.geometry import build_bulk_crystal
from cellulose_forge.hydroxyl import HydroxylHydrogenPlacer


def _place_all(replication: tuple[int, int, int], *, reverse: bool = False) -> tuple[
    dict[tuple[int, int, str], np.ndarray], float
]:
    crystal = build_bulk_crystal(load_allomorph("ii"), replication)
    labels = list(crystal.atom_labels)
    xyz = np.asarray(crystal.cartesian_coordinates)
    lookup = {
        (chain, residue, label): index
        for index, (chain, residue, label) in enumerate(
            zip(crystal.chain_numbers, crystal.residue_numbers, labels, strict=True)
        )
    }
    donors = [key for key in lookup if key[2] in {"O2", "O3", "O6"}]
    if reverse:
        donors.reverse()
    start = time.perf_counter()
    placer = HydroxylHydrogenPlacer(labels, xyz, crystal.basis_vectors)
    for number, key in enumerate(donors):
        carbon = xyz[lookup[(key[0], key[1], "C" + key[2][1])]]
        oxygen = xyz[lookup[key]]
        placer.request(
            number,
            carbon,
            oxygen,
            2 * oxygen - carbon,
            donor=key[2],
            preferred_acceptor="O5" if key[2] == "O3" else None,
        )
    placed = placer.resolve()
    elapsed = time.perf_counter() - start
    return {key: placed[number] for number, key in enumerate(donors)}, elapsed


def test_placement_does_not_depend_on_request_order() -> None:
    forward, _ = _place_all((2, 2, 2))
    backward, _ = _place_all((2, 2, 2), reverse=True)
    assert forward.keys() == backward.keys()
    for key, coordinate in forward.items():
        np.testing.assert_allclose(coordinate, backward[key], atol=1e-9)


def test_placement_cost_grows_linearly() -> None:
    small_placed, small = _place_all((6, 5, 2))
    large_placed, large = _place_all((6, 5, 8))
    per_small = small / len(small_placed)
    per_large = large / len(large_placed)
    # Four times the hydroxyls; quadratic placement would cost ~4x more per hydroxyl.
    assert per_large < 2.0 * per_small


class _BruteForcePlacer(HydroxylHydrogenPlacer):
    """Reference placer whose single bin makes every atom a neighbour."""

    def _bin_shape(self) -> tuple[int, int, int]:
        return (1, 1, 1)


def _requests(crystal: object) -> list[tuple[int, np.ndarray, np.ndarray, str]]:
    labels = list(crystal.atom_labels)  # type: ignore[attr-defined]
    xyz = np.asarray(crystal.cartesian_coordinates)  # type: ignore[attr-defined]
    lookup = {
        (chain, residue, label): index
        for index, (chain, residue, label) in enumerate(
            zip(crystal.chain_numbers, crystal.residue_numbers, labels, strict=True)  # type: ignore[attr-defined]
        )
    }
    return [
        (number, xyz[lookup[(key[0], key[1], "C" + key[2][1])]], xyz[lookup[key]], key[2])
        for number, key in enumerate(k for k in lookup if k[2] in {"O2", "O3", "O6"})
    ]


def test_cell_list_matches_brute_force_on_triclinic_periodic_boxes() -> None:
    # I-alpha is triclinic (P1); a 3 x 3 x 3 box gives 3 bins per axis, so the
    # cell list wraps neighbours across every periodic face.
    crystal = build_bulk_crystal(load_allomorph("i_alpha"), (3, 3, 3))
    box = np.asarray(crystal.basis_vectors)
    labels = list(crystal.atom_labels)
    xyz = np.asarray(crystal.cartesian_coordinates)
    # Shift the whole crystal by a non-lattice vector so atoms straddle faces.
    shift = 0.37 * box[0] - 0.61 * box[1] + 0.83 * box[2]
    for coordinates in (xyz, xyz + shift):
        results = []
        for placer_type in (HydroxylHydrogenPlacer, _BruteForcePlacer):
            placer = placer_type(labels, coordinates, box)
            for number, carbon, oxygen, donor in _requests(crystal):
                offset = coordinates[0] - xyz[0]
                placer.request(
                    number,
                    carbon + offset,
                    oxygen + offset,
                    2 * (oxygen + offset) - (carbon + offset),
                    donor=donor,
                    preferred_acceptor="O5" if donor == "O3" else None,
                )
            results.append(placer.resolve())
        cell_list, brute_force = results
        assert cell_list.keys() == brute_force.keys()
        for key, coordinate in cell_list.items():
            np.testing.assert_allclose(coordinate, brute_force[key], atol=1e-9)
