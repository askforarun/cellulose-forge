"""Compare native finite CHARMM36 topology with an independent CHARMM-GUI ITP."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from cellulose_forge.charmm36 import (
    Charmm36ChainTemplate,
    FiniteCharmm36Topology,
    native_finite_charmm36_topology,
    read_charmm36_chain_template,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _reference_interactions(
    topology: Charmm36ChainTemplate,
) -> set[tuple[str, tuple[tuple[int, str], ...], tuple[str, ...]]]:
    result = set()
    for interaction in topology.interactions:
        atoms = tuple(
            (
                topology.atoms[index].residue_number,
                topology.atoms[index].atom_name,
            )
            for index in interaction.atom_indices
        )
        result.add((interaction.section, min(atoms, atoms[::-1]), interaction.parameters))
    return result


def _generated_interactions(
    topology: FiniteCharmm36Topology,
) -> set[tuple[str, tuple[tuple[int, str], ...], tuple[str, ...]]]:
    result = set()
    for interaction in topology.interactions:
        atoms = tuple(
            (
                topology.topology_atoms[index].residue_number,
                topology.topology_atoms[index].atom_name,
            )
            for index in interaction.atom_indices
        )
        result.add((interaction.section, min(atoms, atoms[::-1]), interaction.parameters))
    return result


def validate(force_field_directory: Path, itp_path: Path) -> dict[str, object]:
    """Return exact atom and bonded-graph comparison results."""
    reference = read_charmm36_chain_template(
        itp_path, require_uniform_residues=False
    )
    generated = native_finite_charmm36_topology(
        force_field_directory, degree_polymerization=reference.degree_polymerization
    )
    reference_atoms = {
        (atom.residue_number, atom.atom_name): (atom.atom_type, atom.charge)
        for atom in reference.atoms
    }
    generated_atoms = {
        (atom.residue_number, atom.atom_name): (atom.atom_type, atom.charge)
        for atom in generated.topology_atoms
    }
    if reference_atoms != generated_atoms:
        common_atoms = reference_atoms.keys() & generated_atoms.keys()
        changed_atoms = sum(
            reference_atoms[key] != generated_atoms[key] for key in common_atoms
        )
        raise ValueError(
            "generated finite atoms differ from CHARMM-GUI: "
            f"missing={len(reference_atoms.keys() - generated_atoms.keys())}, "
            f"extra={len(generated_atoms.keys() - reference_atoms.keys())}, "
            f"changed={changed_atoms}"
        )
    reference_graph = _reference_interactions(reference)
    generated_graph = _generated_interactions(generated)
    if reference_graph != generated_graph:
        raise ValueError(
            "generated finite bonded graph differs from CHARMM-GUI: "
            f"missing={len(reference_graph - generated_graph)}, "
            f"extra={len(generated_graph - reference_graph)}"
        )
    counts = {
        section: sum(item.section == section for item in reference.interactions)
        for section in ("bonds", "pairs", "angles", "dihedrals")
    }
    return {
        "status": "passed",
        "reference": "doGlycans cellulose example processed independently by CHARMM-GUI v1.7",
        "reference_itp_sha256": _sha256(itp_path),
        "carb_rtp_sha256": _sha256(force_field_directory / "carb.rtp"),
        "reference_dp": reference.degree_polymerization,
        "atom_count": len(reference.atoms),
        "chain_charge_e": reference.total_charge,
        "atom_names_types_charges_match": True,
        "interaction_counts": counts,
        "bonded_graph_matches": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force-field-directory", type=Path, required=True)
    parser.add_argument("--itp", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = validate(args.force_field_directory, args.itp)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
