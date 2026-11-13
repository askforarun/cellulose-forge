"""Record the hydrogen-bond network formed by modelled cellulose II hydroxyls.

The cellulose II diffraction model (Langan, Nishiyama and Chanzy 2001) has no
hydroxyl-hydrogen coordinates, so CelluloseForge models HO2, HO3 and HO6.
This script reports, separately for the origin and center chains, which
acceptor each modelled hydroxyl hydrogen-bonds to and the distance and angle
distributions, so the result can be compared with the published network:

* intrachain O3-H...O5 to the next residue in both chains, the only
  intrachain hydrogen bond in cellulose II; and
* intermolecular O2-H...O6, O6-H...O6 and O2-H...O2 bonds, including bonds
  between sheets (Langan, Nishiyama and Chanzy, J. Am. Chem. Soc. 1999, 121,
  9940; the OH2 and OH6 directions are reported to be the least certain).

It uses the committed periodic GLYCAM reference chain and needs no external
programs.  It is a measurement, not a pass/fail test of the literature
assignment, which the regression tests check separately.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "validation"))

from _provenance import provenance  # noqa: E402

from cellulose_forge.crystallography import load_allomorph  # noqa: E402
from cellulose_forge.geometry import build_bulk_crystal  # noqa: E402
from cellulose_forge.glycam_mapping import map_crystal_to_periodic_glycam  # noqa: E402
from cellulose_forge.hydroxyl import HydroxylHydrogenPlacer  # noqa: E402

PERIODIC_REFERENCE_PDB = (
    ROOT / "examples/glycam/dp12_periodic_reference/cellulose_periodic_chain.pdb"
)
REPLICATION = (2, 2, 6)
DONORS = {"H2O": "O2", "H3O": "O3", "H6O": "O6"}
MAXIMUM_HYDROGEN_ACCEPTOR_ANGSTROM = 2.6
MINIMUM_ANGLE_DEGREES = 120.0
# Origin chains sit at fractional (0, 0) of the unit cell, center chains at (1/2, 1/2).
ORIGIN_TOLERANCE = 0.25


def _chain_roles(mapped: object, replication: tuple[int, int, int]) -> dict[int, str]:
    """Classify chains as origin (fractional a,b near 0) or center (near 1/2)."""
    inverse = np.linalg.inv(np.asarray(mapped.basis_vectors))  # type: ignore[attr-defined]
    roles = {}
    chains = np.asarray(mapped.chain_numbers)  # type: ignore[attr-defined]
    xyz = np.asarray(mapped.cartesian_coordinates)  # type: ignore[attr-defined]
    names = np.asarray(mapped.atom_names)  # type: ignore[attr-defined]
    for chain in sorted(set(chains.tolist())):
        ring = xyz[(chains == chain) & (names == "C1")]
        cell = (ring[0] @ inverse) * np.asarray(replication)
        offset = cell[:2] - np.floor(cell[:2])
        is_origin = np.all(np.minimum(offset, 1 - offset) < ORIGIN_TOLERANCE)
        roles[chain] = "origin" if is_origin else "center"
    return roles


def _closest_contacts(
    names: np.ndarray,
    chains: np.ndarray,
    residues: np.ndarray,
    xyz: np.ndarray,
    box: np.ndarray,
) -> dict[str, float]:
    """Closest H, O and C to any modelled hydrogen, excluding its own O and C."""
    inverse = np.linalg.inv(box)
    oxygen = np.flatnonzero(np.char.startswith(names.astype(str), "O"))
    hydrogens = np.flatnonzero(np.char.startswith(names.astype(str), "H"))
    contacts = {"hydrogen": np.inf, "oxygen": np.inf, "carbon": np.inf}
    for index in np.flatnonzero(np.isin(names, list(DONORS))):
        donor_name = DONORS[str(names[index])]
        same = (chains == chains[index]) & (residues == residues[index])
        own = np.flatnonzero(same & np.isin(names, [donor_name, "C" + donor_name[1:]]))
        fractional = (xyz - xyz[index]) @ inverse
        fractional -= np.round(fractional)
        distance = np.linalg.norm(fractional @ box, axis=1)
        distance[np.concatenate(([index], own))] = np.inf
        contacts["hydrogen"] = min(contacts["hydrogen"], float(distance[hydrogens].min()))
        contacts["oxygen"] = min(contacts["oxygen"], float(distance[oxygen].min()))
        carbons = np.flatnonzero(np.char.startswith(names.astype(str), "C"))
        contacts["carbon"] = min(contacts["carbon"], float(distance[carbons].min()))

    return contacts


def measure() -> dict[str, object]:
    crystal = build_bulk_crystal(load_allomorph("ii"), REPLICATION)
    mapped = map_crystal_to_periodic_glycam(crystal, PERIODIC_REFERENCE_PDB)
    box = np.asarray(mapped.basis_vectors)
    inverse = np.linalg.inv(box)
    names = np.asarray(mapped.atom_names)
    chains = np.asarray(mapped.chain_numbers)
    residues = np.asarray(mapped.residue_numbers)
    xyz = np.asarray(mapped.cartesian_coordinates)
    degree = int(residues.max())
    roles = _chain_roles(mapped, REPLICATION)
    oxygen = np.flatnonzero(np.char.startswith(names.astype(str), "O"))

    records: dict[str, dict[str, list[dict[str, object]]]] = {}
    for index in np.flatnonzero(np.isin(names, list(DONORS))):
        donor_name = DONORS[str(names[index])]
        same = (chains == chains[index]) & (residues == residues[index])
        donor = int(np.flatnonzero(same & (names == donor_name))[0])
        fractional = (xyz[oxygen] - xyz[index]) @ inverse
        fractional -= np.round(fractional)
        delta = fractional @ box
        distance = np.linalg.norm(delta, axis=1)
        distance[oxygen == donor] = np.inf
        to_donor = xyz[donor] - xyz[index]
        cosine = (delta @ to_donor) / (distance * np.linalg.norm(to_donor))
        angle = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
        valid = (distance <= MAXIMUM_HYDROGEN_ACCEPTOR_ANGSTROM) & (
            angle >= MINIMUM_ANGLE_DEGREES
        )
        role = roles[int(chains[index])]
        key = f"{donor_name}-H"
        entry: dict[str, object]
        if valid.any():
            best = int(np.argmin(np.where(valid, distance, np.inf)))
            acceptor = int(oxygen[best])
            same_chain = bool(chains[acceptor] == chains[index])
            offset = int(
                (residues[acceptor] - residues[index] + degree // 2) % degree - degree // 2
            )
            partner = (
                f"{names[acceptor]} (same chain, residue {offset:+d})"
                if same_chain
                else f"{names[acceptor]} ({roles[int(chains[acceptor])]} chain)"
            )
            entry = {
                "acceptor": partner,
                "h_acceptor_angstrom": float(distance[best]),
                "angle_degrees": float(angle[best]),
            }
        else:
            entry = {"acceptor": "none"}
        records.setdefault(role, {}).setdefault(key, []).append(entry)

    contacts = _closest_contacts(names, chains, residues, xyz, box)

    summary: dict[str, object] = {}
    for role, donors in sorted(records.items()):
        summary[role] = {}
        for donor, entries in sorted(donors.items()):
            bonded = [entry for entry in entries if entry["acceptor"] != "none"]
            summary[role][donor] = {  # type: ignore[index]
                "hydroxyls": len(entries),
                "acceptors": dict(Counter(str(entry["acceptor"]) for entry in entries)),
                "h_acceptor_angstrom_range": (
                    [
                        min(float(entry["h_acceptor_angstrom"]) for entry in bonded),
                        max(float(entry["h_acceptor_angstrom"]) for entry in bonded),
                    ]
                    if bonded
                    else None
                ),
                "angle_degrees_range": (
                    [
                        min(float(entry["angle_degrees"]) for entry in bonded),
                        max(float(entry["angle_degrees"]) for entry in bonded),
                    ]
                    if bonded
                    else None
                ),
            }
    return {
        "system": f"cellulose II periodic bulk {REPLICATION}, GLYCAM mapping",
        "provenance": provenance(("numpy",)),
        "criteria": {
            "maximum_h_acceptor_angstrom": MAXIMUM_HYDROGEN_ACCEPTOR_ANGSTROM,
            "minimum_o_h_o_angle_degrees": MINIMUM_ANGLE_DEGREES,
        },
        "literature_network": {
            "source": "Langan, Nishiyama and Chanzy, J. Am. Chem. Soc. 1999, 121, 9940",
            "intrachain": "O3-H...O5 (next residue) in both chains",
            "intermolecular": "O2-H...O6, O6-H...O6 and O2-H...O2, including intersheet",
            "least_certain": "OH2 and OH6 directions",
        },
        "chains": summary,
        "closest_modelled_hydrogen_contacts_angstrom": contacts,
    }


SCALING_REPLICATIONS = ((6, 5, 2), (6, 5, 4), (6, 5, 8), (6, 5, 16))


def scaling() -> list[dict[str, float | int]]:
    """Time placement of every cellulose II hydroxyl at increasing sizes."""
    rows: list[dict[str, float | int]] = []
    for replication in SCALING_REPLICATIONS:
        crystal = build_bulk_crystal(load_allomorph("ii"), replication)
        labels = list(crystal.atom_labels)
        xyz = np.asarray(crystal.cartesian_coordinates)
        lookup = {
            (chain, residue, label): index
            for index, (chain, residue, label) in enumerate(
                zip(crystal.chain_numbers, crystal.residue_numbers, labels, strict=True)
            )
        }
        start = time.perf_counter()
        placer = HydroxylHydrogenPlacer(labels, xyz, crystal.basis_vectors)
        count = 0
        for (chain, residue, label), index in lookup.items():
            if label in {"O2", "O3", "O6"}:
                carbon = xyz[lookup[(chain, residue, "C" + label[1])]]
                placer.request(
                    count,
                    carbon,
                    xyz[index],
                    2 * xyz[index] - carbon,
                    donor=label,
                    preferred_acceptor="O5" if label == "O3" else None,
                )
                count += 1
        placer.resolve()
        seconds = time.perf_counter() - start
        rows.append(
            {
                "atoms": len(xyz),
                "hydroxyls": count,
                "seconds": seconds,
                "microseconds_per_hydroxyl": 1.0e6 * seconds / count,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation/results/cellulose_ii_hydroxyl_network.json",
    )
    args = parser.parse_args()
    result = measure()
    result["placement_scaling"] = scaling()
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["chains"], indent=2))


if __name__ == "__main__":
    main()
