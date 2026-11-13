"""Prototype a covalently periodic DP6 OpenFF cellulose repeat."""

from __future__ import annotations

import hashlib
import json
from importlib import metadata
from pathlib import Path

from openff.toolkit import ForceField, Molecule
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
FINITE = ROOT / "validation/results/openff_cellodextrin_dp6_nagl_charged.json"
OUTPUT = ROOT / "validation/results/openff_cellodextrin_dp6_periodic_nagl_charged.json"
FIRST_RESIDUE = 2
LAST_RESIDUE = 7
REFERENCE_INTERNAL_RESIDUE = 3
TERMINAL_ATOM_COUNT = 3


def periodic_dp6_template() -> Molecule:
    """Close the finite DP6 graph into a six-residue periodic repeat."""
    finite = Molecule.from_json(FINITE.read_text(encoding="utf-8"))
    rdkit = Chem.RWMol(finite.to_rdkit())
    first_c1 = next(
        atom.molecule_atom_index
        for atom in finite.atoms
        if atom.name == "C1" and int(atom.metadata["residue_number"]) == FIRST_RESIDUE
    )
    last_o4 = next(
        atom.molecule_atom_index
        for atom in finite.atoms
        if atom.name == "O4" and int(atom.metadata["residue_number"]) == LAST_RESIDUE
    )
    terminal_indices = [
        atom.molecule_atom_index
        for atom in finite.atoms
        if (atom.name, int(atom.metadata["residue_number"]))
        in {("HO1", 1), ("O1", 1), ("H4O", 7)}
    ]
    if len(terminal_indices) != TERMINAL_ATOM_COUNT:
        raise RuntimeError(f"Expected three finite terminal atoms: {terminal_indices}")
    rdkit.AddBond(last_o4, first_c1, Chem.BondType.SINGLE)
    for index in sorted(terminal_indices, reverse=True):
        rdkit.RemoveAtom(index)

    molecule = Molecule.from_rdkit(
        rdkit.GetMol(),
        allow_undefined_stereo=False,
        hydrogens_are_explicit=True,
    )
    internal_names = [
        atom.name
        for atom in molecule.atoms
        if int(atom.metadata["residue_number"]) == FIRST_RESIDUE
    ]
    atom_order = [
        next(
            atom.molecule_atom_index
            for atom in molecule.atoms
            if int(atom.metadata["residue_number"]) == residue and atom.name == name
        )
        for residue in range(FIRST_RESIDUE, 8)
        for name in internal_names
    ]
    molecule = Molecule.from_rdkit(
        Chem.RenumberAtoms(molecule.to_rdkit(), atom_order),
        allow_undefined_stereo=False,
        hydrogens_are_explicit=True,
    )
    reference_stereo = {
        atom.name: atom.stereochemistry
        for atom in finite.atoms
        if int(atom.metadata["residue_number"]) == REFERENCE_INTERNAL_RESIDUE
        and atom.name in {"C1", "C2", "C3", "C4", "C5"}
    }
    if reference_stereo != {"C1": "S", "C2": "R", "C3": "R", "C4": "S", "C5": "R"}:
        raise RuntimeError(f"Unexpected PubChem internal-residue stereo: {reference_stereo}")
    for atom in molecule.atoms:
        if atom.name in reference_stereo:
            atom._stereochemistry = reference_stereo[atom.name]
    molecule.assign_partial_charges("openff-gnn-am1bcc-1.0.0.pt")
    return molecule


def main() -> None:
    molecule = periodic_dp6_template()
    force_field = ForceField("openff-2.3.0.offxml")
    labels = force_field.label_molecules(molecule.to_topology())[0]
    force_field.create_openmm_system(molecule.to_topology())
    OUTPUT.write_text(molecule.to_json(indent=2) + "\n", encoding="utf-8")
    cache_hash = hashlib.sha256(OUTPUT.read_bytes()).hexdigest()
    print(
        json.dumps(
            {
                "atoms": molecule.n_atoms,
                "bonds": molecule.n_bonds,
                "charge_e": float(molecule.total_charge.m),
                "stereochemistry_source": "PubChem-derived finite DP6 internal residue 3",
                "stereochemistry_by_atom_name": {
                    name: next(
                        atom.stereochemistry for atom in molecule.atoms if atom.name == name
                    )
                    for name in ("C1", "C2", "C3", "C4", "C5")
                },
                "charge_provenance": {
                    "method": "openff-gnn-am1bcc-1.0.0.pt",
                    "openff_toolkit": metadata.version("openff-toolkit"),
                    "openff_nagl": metadata.version("openff-nagl"),
                    "openff_nagl_models": metadata.version("openff-nagl-models"),
                    "template_sha256": cache_hash,
                },
                "parameter_counts": {name: len(values) for name, values in labels.items()},
                "output": str(OUTPUT.relative_to(ROOT)),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
