"""Chain-specific hydroxyl acceptors in cellulose II."""

from __future__ import annotations

import numpy as np

from cellulose_forge.crystallography import FloatArray
from cellulose_forge.geometry import BulkCrystal


def _minimum_image(delta: FloatArray, box: FloatArray) -> FloatArray:
    fractional = np.asarray(delta, dtype=np.float64) @ np.linalg.inv(box)
    fractional -= np.round(fractional)
    return np.asarray(fractional @ box, dtype=np.float64)


def _positive_crystallographic_c(crystal: BulkCrystal) -> FloatArray:
    if any("chains aligned with Cartesian -y" in note for note in crystal.annotations):
        return np.asarray((0.0, -1.0, 0.0))
    vector = np.asarray(crystal.basis_vectors[2], dtype=np.float64)
    return vector / np.linalg.norm(vector)


def cellulose_ii_chain_roles(crystal: BulkCrystal) -> dict[int, str]:
    """Classify cellulose-II chains using Langan's O5/C6 convention.

    Langan, Nishiyama and Chanzy define the origin (up) chain by O5 lying
    above C6 along crystallographic +c; the centre (down) chain has the
    opposite relation. This remains valid after CelluloseForge's rigid slab
    and fibril orientations.
    """
    if crystal.allomorph != "II":
        msg = "cellulose-II chain roles require allomorph II"
        raise ValueError(msg)
    labels = np.asarray(crystal.atom_labels)
    chains = np.asarray(crystal.chain_numbers)
    residues = np.asarray(crystal.residue_numbers)
    coordinates = np.asarray(crystal.cartesian_coordinates)
    c_axis = _positive_crystallographic_c(crystal)
    roles: dict[int, str] = {}
    for chain in sorted(set(chains)):
        residue = min(residues[chains == chain])
        mask = (chains == chain) & (residues == residue)
        lookup = {
            name: coordinates[index]
            for index, name in zip(np.flatnonzero(mask), labels[mask], strict=True)
        }
        direction = _minimum_image(lookup["O5"] - lookup["C6"], crystal.basis_vectors)
        roles[int(chain)] = "origin" if float(np.dot(direction, c_axis)) > 0.0 else "center"
    return roles


def cellulose_ii_preferred_acceptor(
    crystal: BulkCrystal,
    *,
    chain_number: int,
    residue_number: int,
    donor: str,
) -> FloatArray:
    """Return the Table-3 acceptor nearest one cellulose-II hydroxyl donor.

    The major components from Langan et al., JACS 1999, 121, 9940 are:
    O2o-H...O6o, O6o-H...O6c, O2c-H...O2o, O6c-H...O2c, and
    O3-H...O5 on the same chain.
    """
    roles = cellulose_ii_chain_roles(crystal)
    role = roles[chain_number]
    target_by_donor = {
        ("origin", "O2"): ("origin", "O6"),
        ("origin", "O6"): ("center", "O6"),
        ("center", "O2"): ("origin", "O2"),
        ("center", "O6"): ("center", "O2"),
    }
    labels = np.asarray(crystal.atom_labels)
    chains = np.asarray(crystal.chain_numbers)
    residues = np.asarray(crystal.residue_numbers)
    coordinates = np.asarray(crystal.cartesian_coordinates)
    donor_mask = (chains == chain_number) & (residues == residue_number) & (labels == donor)
    if int(donor_mask.sum()) != 1:
        msg = f"expected one {donor} in cellulose-II chain {chain_number} residue {residue_number}"
        raise ValueError(msg)
    donor_coordinate = coordinates[np.flatnonzero(donor_mask)[0]]
    if donor == "O3":
        target_mask = (chains == chain_number) & (labels == "O5") & (residues != residue_number)
    else:
        try:
            target_role, target_label = target_by_donor[role, donor]
        except KeyError as error:
            msg = f"unsupported cellulose-II hydroxyl donor {role} {donor}"
            raise ValueError(msg) from error
        target_chains = np.asarray([roles[int(chain)] == target_role for chain in chains])
        target_mask = target_chains & (labels == target_label)
    candidates = coordinates[target_mask]
    if not len(candidates):
        msg = f"no cellulose-II acceptor found for {role} {donor}"
        raise ValueError(msg)
    deltas = np.vstack(
        [
            _minimum_image(candidate - donor_coordinate, crystal.basis_vectors)
            for candidate in candidates
        ]
    )
    return np.asarray(donor_coordinate + deltas[np.argmin(np.linalg.norm(deltas, axis=1))])
