"""Cellulose-II chain roles and neutron-derived hydroxyl network."""

from __future__ import annotations

import numpy as np
import pytest

from cellulose_forge.cellulose_ii import (
    cellulose_ii_chain_roles,
    cellulose_ii_preferred_acceptor,
)
from cellulose_forge.crystallography import load_allomorph
from cellulose_forge.geometry import build_bulk_crystal


def _minimum_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    fractional = delta @ np.linalg.inv(box)
    fractional -= np.round(fractional)
    return fractional @ box


def test_chain_roles_follow_crystallographic_origin_and_center_labels() -> None:
    crystal = build_bulk_crystal(load_allomorph("ii"), (2, 2, 2))
    assert cellulose_ii_chain_roles(crystal) == {
        1: "origin",
        2: "center",
        3: "origin",
        4: "center",
        5: "origin",
        6: "center",
        7: "origin",
        8: "center",
    }


@pytest.mark.parametrize(
    ("chain", "donor", "acceptor", "acceptor_role", "distance"),
    (
        (1, "O2", "O6", "origin", 2.776),
        (1, "O6", "O6", "center", 2.796),
        (2, "O2", "O2", "origin", 2.916),
        (2, "O6", "O2", "center", 2.526),
    ),
)
def test_preferred_acceptors_follow_langan_table_3(
    chain: int,
    donor: str,
    acceptor: str,
    acceptor_role: str,
    distance: float,
) -> None:
    crystal = build_bulk_crystal(load_allomorph("ii"), (2, 2, 2))
    roles = cellulose_ii_chain_roles(crystal)
    labels = np.asarray(crystal.atom_labels)
    chains = np.asarray(crystal.chain_numbers)
    residues = np.asarray(crystal.residue_numbers)
    xyz = np.asarray(crystal.cartesian_coordinates)
    residue = int(residues[chains == chain].min())
    donor_coordinate = xyz[(chains == chain) & (residues == residue) & (labels == donor)][0]
    preferred = cellulose_ii_preferred_acceptor(
        crystal,
        chain_number=chain,
        residue_number=residue,
        donor=donor,
    )
    assert np.linalg.norm(
        _minimum_image(preferred - donor_coordinate, crystal.basis_vectors)
    ) == pytest.approx(distance, abs=0.001)
    candidates = np.flatnonzero(labels == acceptor)
    matched = min(
        candidates,
        key=lambda index: np.linalg.norm(
            _minimum_image(xyz[index] - preferred, crystal.basis_vectors)
        ),
    )
    assert roles[int(chains[matched])] == acceptor_role
