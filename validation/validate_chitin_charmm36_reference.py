"""Compare CelluloseForge's finite chitin CHARMM36 topology with a psfgen reference.

The reference PSF is built outside CelluloseForge from the canonical CHARMM
``top_all36_carb.rtf`` (``RESI BGLCNA`` and ``PRES 14bb``) with VMD psfgen:

    package require psfgen
    topology top_all36_carb.rtf
    segment CHT { first none; last none; auto none;
                  for {set i 1} {$i <= 6} {incr i} { residue $i BGLCNA } }
    for {set i 1} {$i < 6} {incr i} { patch 14bb CHT:$i CHT:[expr {$i + 1}] }
    regenerate angles dihedrals
    writepsf charmm glcnac_dp6.psf

psfgen is used only to produce this validation reference. CelluloseForge does
not import or call it. The CHARMM archive is not redistributed; its filename
and SHA-256 go into the result JSON.

Usage:
    python validation/validate_chitin_charmm36_reference.py \
        --psf glcnac_dp6.psf --archive toppar_c36_feb26.tgz \
        --force-field-directory .../charmm36-feb2026_cgenff-5.0.ff \
        --output chitin_charmm36_reference.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from cellulose_forge.charmm36 import native_finite_charmm36_topology

_CHARGE_TOLERANCE = 1e-6
_MASS_TOLERANCE = 1e-3
_PSF_SECTIONS = {"!NATOM": 0, "!NBOND": 2, "!NTHETA": 3, "!NPHI": 4, "!NIMPHI": 4}

Key = tuple[int, str]


def _read_psf(path: Path) -> dict[str, list]:
    lines = path.read_text(encoding="utf-8").splitlines()
    data: dict[str, list] = {}
    index = 0
    while index < len(lines):
        fields = lines[index].split()
        tag = next((f.split(":")[0] for f in fields if f.startswith("!N")), None)
        if tag not in _PSF_SECTIONS:
            index += 1
            continue
        count = int(fields[0])
        index += 1
        if tag == "!NATOM":
            data["atoms"] = [lines[index + i].split() for i in range(count)]
            index += count
            continue
        arity = _PSF_SECTIONS[tag]
        numbers: list[int] = []
        while len(numbers) < count * arity:
            numbers.extend(int(value) for value in lines[index].split())
            index += 1
        data[tag] = [tuple(numbers[i : i + arity]) for i in range(0, len(numbers), arity)]
    return data


def _canonical(row: tuple[Key, ...]) -> tuple[Key, ...]:
    return min(row, tuple(reversed(row)))


def _graph_pairs(bonds: set[tuple[Key, ...]], dihedrals: set[tuple[Key, ...]]) -> set:
    neighbours: dict[Key, set[Key]] = defaultdict(set)
    for first, second in bonds:
        neighbours[first].add(second)
        neighbours[second].add(first)
    pairs = set()
    for row in dihedrals:
        first, last = row[0], row[-1]
        close = neighbours[first] | {
            atom for neighbour in neighbours[first] for atom in neighbours[neighbour]
        }
        if last not in close and last != first:
            pairs.add(_canonical((first, last)))
    return pairs


def _compare(name: str, reference: set, generated: set) -> dict[str, object]:
    return {
        "term": name,
        "reference": len(reference),
        "generated": len(generated),
        "missing": sorted(map(str, reference - generated))[:10],
        "extra": sorted(map(str, generated - reference))[:10],
        "match": reference == generated,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--psf", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True, help="archive path or name")
    parser.add_argument(
        "--archive-sha256", help="checksum to record when the archive is no longer on disk"
    )
    parser.add_argument("--force-field-directory", type=Path, required=True)
    parser.add_argument("--degree-polymerization", type=int, default=6)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    psf = _read_psf(args.psf)
    ref_atoms = {
        (int(fields[2]), fields[4]): (fields[5], float(fields[6]), float(fields[7]))
        for fields in psf["atoms"]
    }
    key_by_serial = {int(fields[0]): (int(fields[2]), fields[4]) for fields in psf["atoms"]}

    def keys(rows: list[tuple[int, ...]]) -> list[tuple[Key, ...]]:
        return [tuple(key_by_serial[serial] for serial in row) for row in rows]

    topology = native_finite_charmm36_topology(
        args.force_field_directory,
        degree_polymerization=args.degree_polymerization,
        residue="BGLCNA",
    )
    gen_atoms = {
        (atom.residue_number, atom.atom_name): (atom.atom_type, atom.charge, atom.mass)
        for atom in topology.topology_atoms
    }
    key_by_index = [(atom.residue_number, atom.atom_name) for atom in topology.topology_atoms]
    sections: dict[str, set] = defaultdict(set)
    for interaction in topology.interactions:
        row = tuple(key_by_index[index] for index in interaction.atom_indices)
        if interaction.section == "dihedrals" and interaction.parameters == ("2",):
            sections["impropers"].add(row)
        else:
            sections[interaction.section].add(_canonical(row))

    ref_bonds = {_canonical(row) for row in keys(psf["!NBOND"])}
    ref_dihedrals = {_canonical(row) for row in keys(psf["!NPHI"])}
    results = [
        _compare("atoms", set(ref_atoms), set(gen_atoms)),
        _compare("bonds", ref_bonds, sections["bonds"]),
        _compare("angles", {_canonical(r) for r in keys(psf["!NTHETA"])}, sections["angles"]),
        _compare("proper_dihedrals", ref_dihedrals, sections["dihedrals"]),
        _compare("impropers", set(keys(psf["!NIMPHI"])), sections["impropers"]),
        _compare("pairs_1_4", _graph_pairs(ref_bonds, ref_dihedrals), sections["pairs"]),
    ]

    mismatched = [
        str(key)
        for key in sorted(set(ref_atoms) & set(gen_atoms))
        if ref_atoms[key][0] != gen_atoms[key][0]
        or abs(ref_atoms[key][1] - gen_atoms[key][1]) > _CHARGE_TOLERANCE
        or abs(ref_atoms[key][2] - gen_atoms[key][2]) > _MASS_TOLERANCE
    ]
    results.append(
        {"term": "types_charges_masses", "mismatched": mismatched[:10], "match": not mismatched}
    )

    def residue_summary(atoms: dict[Key, tuple[str, float, float]]) -> dict[int, list]:
        summary: dict[int, list] = defaultdict(lambda: [0, 0.0])
        for (residue, _), (_, charge, _) in atoms.items():
            summary[residue][0] += 1
            summary[residue][1] += charge
        return {residue: [n, round(q, 6)] for residue, (n, q) in sorted(summary.items())}

    ref_residues = residue_summary(ref_atoms)
    gen_residues = residue_summary(gen_atoms)
    results.append(
        {
            "term": "per_residue_atom_count_and_charge",
            "reference": ref_residues,
            "generated": gen_residues,
            "total_charge": round(sum(q for _, q, _ in gen_atoms.values()), 6),
            "match": ref_residues == gen_residues,
        }
    )
    termini = {
        "reducing_end_names": sorted(name for r, name in gen_atoms if r == 1),
        "nonreducing_end_names": sorted(
            name for r, name in gen_atoms if r == args.degree_polymerization
        ),
    }
    ref_termini = {
        "reducing_end_names": sorted(name for r, name in ref_atoms if r == 1),
        "nonreducing_end_names": sorted(
            name for r, name in ref_atoms if r == args.degree_polymerization
        ),
    }
    results.append({"term": "termini", "generated": termini, "match": termini == ref_termini})

    passed = all(result["match"] for result in results)
    report = {
        "test": "finite BGLCNA beta(1->4) chain vs psfgen reference",
        "degree_polymerization": args.degree_polymerization,
        "reference_generator": (
            "VMD 1.9.4a57 psfgen 2.0 (validation only, not a runtime dependency)"
        ),
        "reference_topology": "top_all36_carb.rtf: RESI BGLCNA, PRES 14bb",
        "charmm_archive": args.archive.name,
        "charmm_archive_sha256": (
            hashlib.sha256(args.archive.read_bytes()).hexdigest()
            if args.archive.is_file()
            else args.archive_sha256
        ),
        "gromacs_force_field": args.force_field_directory.name,
        "exclusions": "nrexcl=3 from the bond graph; equal bond sets imply equal exclusions",
        "results": results,
        "passed": passed,
    }
    text = json.dumps(report, indent=2)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
