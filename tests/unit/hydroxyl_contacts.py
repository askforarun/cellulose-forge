"""Shared contact metrics for modelled cellulose II hydroxyl hydrogens."""

from __future__ import annotations

import numpy as np


def modelled_hydroxyl_contacts(mapped: object, donors: dict[str, str]) -> dict[str, float]:
    """Closest contacts and geometric hydrogen bonding of modelled hydroxyl H."""
    names = np.asarray(mapped.atom_names)  # type: ignore[attr-defined]
    xyz = np.asarray(mapped.cartesian_coordinates)  # type: ignore[attr-defined]
    box = np.asarray(mapped.basis_vectors)  # type: ignore[attr-defined]
    chains = np.asarray(mapped.chain_numbers)  # type: ignore[attr-defined]
    residues = np.asarray(mapped.residue_numbers)  # type: ignore[attr-defined]
    inverse = np.linalg.inv(box)
    elements = np.asarray([name[:1] for name in names])
    closest = {"hydrogen": np.inf, "oxygen": np.inf, "other_heavy": np.inf}
    bonded = 0
    o3_donors = 0
    o3_to_intrachain_o5 = 0
    hydroxyls = [index for index, name in enumerate(names) if name in donors]
    for index in hydroxyls:
        same_residue = (chains == chains[index]) & (residues == residues[index])
        donor_name = donors[names[index]]
        donor = int(np.flatnonzero((names == donor_name) & same_residue)[0])
        # The donor's own carbon is a fixed 1-3 neighbour (~1.95 A at 108.5 degrees).
        carbon = int(np.flatnonzero((names == "C" + donor_name[1:]) & same_residue)[0])
        fractional = (xyz - xyz[index]) @ inverse
        fractional -= np.round(fractional)
        delta = fractional @ box
        distance = np.linalg.norm(delta, axis=1)
        distance[[index, donor, carbon]] = np.inf
        closest["hydrogen"] = min(closest["hydrogen"], float(distance[elements == "H"].min()))
        closest["oxygen"] = min(closest["oxygen"], float(distance[elements == "O"].min()))
        closest["other_heavy"] = min(
            closest["other_heavy"], float(distance[(elements != "H") & (elements != "O")].min())
        )
        to_donor = xyz[donor] - xyz[index]
        partners = []
        for acceptor in np.flatnonzero((elements == "O") & (distance < 2.6)):
            cosine = np.dot(to_donor, delta[acceptor]) / (
                np.linalg.norm(to_donor) * distance[acceptor]
            )
            if cosine <= np.cos(np.radians(120.0)):
                partners.append(int(acceptor))
        if partners:
            bonded += 1
        if donor_name == "O3":
            o3_donors += 1
            o3_to_intrachain_o5 += any(
                names[acceptor] == "O5" and chains[acceptor] == chains[index]
                for acceptor in partners
            )
    return {
        **closest,
        "hydrogen_bonded_fraction": bonded / len(hydroxyls),
        "o3_to_intrachain_o5_fraction": o3_to_intrachain_o5 / o3_donors,
    }
