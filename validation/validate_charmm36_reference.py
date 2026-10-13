"""Validate a periodic CHARMM36 cellulose ITP against BGLC and 14bb.

External CHARMM36 files are supplied by the user and are never copied into the
repository.  This script checks the scientifically important transformation:
the BGLC residue after it has participated on both sides of a cyclic 14bb
beta(1->4) linkage.  It also confirms that CelluloseForge regenerates the
reference topology's bonded interaction graph at the same DP.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
from dataclasses import dataclass
from pathlib import Path

from cellulose_forge.charmm36 import (
    Charmm36ChainTemplate,
    derive_periodic_charmm36_topology,
    read_charmm36_chain_template,
    render_periodic_charmm36_itp,
)


@dataclass(frozen=True)
class RtfAtom:
    atom_type: str
    charge: float


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rtf_blocks(path: Path) -> tuple[dict[str, RtfAtom], list[list[str]]]:
    residue: dict[str, RtfAtom] = {}
    patch_rows: list[list[str]] = []
    block = ""
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = raw_line.split("!", maxsplit=1)[0].split()
        if not fields:
            continue
        if fields[:2] == ["RESI", "BGLC"]:
            block = "BGLC"
            continue
        if fields[:2] == ["PRES", "14bb"]:
            block = "14bb"
            continue
        if fields[0] in {"RESI", "PRES"}:
            block = ""
            continue
        if block == "BGLC" and fields[0] == "ATOM":
            residue[fields[1]] = RtfAtom(fields[2], float(fields[3]))
        elif block == "14bb" and fields[0].lower() in {"dele", "atom"}:
            patch_rows.append(fields)
    if not residue or not patch_rows:
        raise ValueError("could not locate BGLC and PRES 14bb in the CHARMM RTF")
    return residue, patch_rows


def _cyclic_repeat_from_rtf(path: Path) -> dict[str, RtfAtom]:
    atoms, rows = _rtf_blocks(path)
    expected = dict(atoms)
    for fields in rows:
        if fields[0].lower() == "dele" and fields[1].lower() == "atom":
            expected.pop(fields[2][1:], None)
        elif fields[0] == "ATOM":
            expected[fields[1][1:]] = RtfAtom(fields[2], float(fields[3]))
    return expected


def _canonical_interactions(
    topology: Charmm36ChainTemplate,
) -> set[tuple[str, tuple[tuple[int, str], ...], tuple[str, ...]]]:
    canonical: set[tuple[str, tuple[tuple[int, str], ...], tuple[str, ...]]] = set()
    for interaction in topology.interactions:
        atoms = tuple(
            (
                topology.atoms[index].residue_number - 1,
                topology.atoms[index].atom_name,
            )
            for index in interaction.atom_indices
        )
        canonical.add(
            (
                interaction.section,
                min(atoms, tuple(reversed(atoms))),
                interaction.parameters,
            )
        )
    return canonical


def validate(rtf_path: Path, itp_path: Path) -> dict[str, object]:
    """Return a JSON-serializable CHARMM36 reference-validation summary."""
    expected = _cyclic_repeat_from_rtf(rtf_path)
    source = read_charmm36_chain_template(itp_path)
    middle = source.degree_polymerization // 2
    actual = {
        atom.atom_name: RtfAtom(atom.atom_type, atom.charge)
        for atom in source.atoms
        if atom.residue_number == middle
    }
    if actual != expected:
        missing = sorted(expected.keys() - actual.keys())
        extra = sorted(actual.keys() - expected.keys())
        changed = sorted(
            name
            for name in expected.keys() & actual.keys()
            if expected[name] != actual[name]
        )
        raise ValueError(
            "periodic ITP repeat does not match BGLC+14bb: "
            f"missing={missing}, extra={extra}, changed={changed}"
        )

    generated = derive_periodic_charmm36_topology(
        itp_path, degree_polymerization=source.degree_polymerization
    )
    with tempfile.TemporaryDirectory() as temporary_directory:
        generated_path = Path(temporary_directory) / "generated.itp"
        generated_path.write_text(
            render_periodic_charmm36_itp(generated), encoding="utf-8"
        )
        regenerated = read_charmm36_chain_template(generated_path)
    source_graph = _canonical_interactions(source)
    generated_graph = _canonical_interactions(regenerated)
    if source_graph != generated_graph:
        raise ValueError(
            "generated periodic bonded graph differs from the reference: "
            f"missing={len(source_graph - generated_graph)}, "
            f"extra={len(generated_graph - source_graph)}"
        )

    counts = {
        section: sum(item.section == section for item in source.interactions)
        for section in ("bonds", "pairs", "angles", "dihedrals")
    }
    return {
        "status": "passed",
        "rtf_sha256": _sha256(rtf_path),
        "itp_sha256": _sha256(itp_path),
        "reference_dp": source.degree_polymerization,
        "atoms_per_repeat": len(actual),
        "repeat_charge_e": sum(atom.charge for atom in actual.values()),
        "deleted_by_cyclic_14bb": sorted({"O1", "HO1", "HO4"}),
        "modified_atoms": {
            name: {"type": actual[name].atom_type, "charge_e": actual[name].charge}
            for name in ("C1", "C4", "O4")
        },
        "interaction_counts": counts,
        "regenerated_graph_matches_reference": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rtf", type=Path, required=True)
    parser.add_argument("--itp", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = validate(args.rtf, args.itp)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
