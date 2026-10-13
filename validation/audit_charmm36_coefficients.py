"""Audit CHARMM carbohydrate coefficients against a converted GROMACS set."""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

from parmed.charmm import CharmmParameterSet

from cellulose_forge.charmm36 import read_charmm36_chain_template

_MINIMUM_NONBONDED_FIELDS = 4
_MINIMUM_ANGLETYPE_FIELDS = 8


def _gromacs_sections(path: Path) -> dict[str, list[list[str]]]:
    result: dict[str, list[list[str]]] = defaultdict(list)
    section = ""
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.split(";", 1)[0].strip()
        if line.startswith("["):
            section = line.strip("[] ")
        elif line and not line.startswith("#"):
            result[section].append(line.split())
    return result


def _lookup(rows: list[list[str]], types: tuple[str, ...]) -> list[list[str]]:
    width = len(types)
    direct = [row for row in rows if tuple(row[:width]) == types]
    reverse = [row for row in rows if tuple(row[:width]) == tuple(reversed(types))]
    return direct or reverse


def _close(left: float, right: float) -> bool:
    # Converted tables round some energies to five decimal places.
    return math.isclose(left, right, rel_tol=5e-5, abs_tol=5e-6)


def audit(  # noqa: PLR0912, PLR0915
    prm: Path, itp: Path, ff_bonded: Path, ff_nonbonded: Path
) -> dict[str, object]:
    source = CharmmParameterSet(str(prm))
    topology = read_charmm36_chain_template(itp)
    bonded = _gromacs_sections(ff_bonded)
    nonbonded = _gromacs_sections(ff_nonbonded)
    atom_types = sorted({atom.atom_type for atom in topology.atoms})
    failures: list[str] = []

    nb_source: dict[str, tuple[float, float]] = {}
    in_nonbonded = False
    for raw in prm.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = raw.split("!", 1)[0].split()
        if fields and fields[0] == "NONBONDED":
            in_nonbonded = True
            continue
        if in_nonbonded and fields and fields[0] == "END":
            break
        if in_nonbonded and len(fields) >= _MINIMUM_NONBONDED_FIELDS:
            try:
                nb_source[fields[0]] = (abs(float(fields[2])) * 4.184, float(fields[3]))
            except ValueError:
                pass
    nb_gmx = {row[0]: (float(row[-2]), float(row[-1])) for row in nonbonded["atomtypes"]}
    for kind in atom_types:
        epsilon, rmin_half = nb_source[kind]
        sigma_expected = 2.0 * rmin_half * 0.1 / (2.0 ** (1.0 / 6.0))
        sigma, epsilon_gmx = nb_gmx[kind]
        if not (_close(sigma, sigma_expected) and _close(epsilon_gmx, epsilon)):
            failures.append(f"nonbonded {kind}")

    used: dict[str, set[tuple[str, ...]]] = defaultdict(set)
    for interaction in topology.interactions:
        used[interaction.section].add(
            tuple(topology.atoms[index].atom_type for index in interaction.atom_indices)
        )

    for types in used["bonds"]:
        value = source.bond_types.get(types) or source.bond_types.get(tuple(reversed(types)))
        rows = _lookup(bonded["bondtypes"], types)
        if value is None or not rows:
            failures.append(f"missing bond {types}")
            continue
        row = rows[0]
        if not (_close(float(row[3]), value.req / 10) and _close(float(row[4]), value.k * 836.8)):
            failures.append(f"bond {types}")

    for types in used["angles"]:
        value = source.angle_types.get(types) or source.angle_types.get(tuple(reversed(types)))
        rows = _lookup(bonded["angletypes"], types)
        if value is None or not rows:
            failures.append(f"missing angle {types}")
            continue
        row = rows[0]
        urey_bradley = source.urey_bradley_types.get(types) or source.urey_bradley_types.get(
            tuple(reversed(types))
        )
        if urey_bradley is None:
            failures.append(f"missing Urey-Bradley {types}")
            continue
        expected = (
            value.theteq,
            value.k * 8.368,
            urey_bradley.req / 10,
            urey_bradley.k * 836.8,
        )
        converted = zip(row[4:8], expected, strict=True)
        if len(row) < _MINIMUM_ANGLETYPE_FIELDS or not all(
            _close(float(actual), target) for actual, target in converted
        ):
            failures.append(f"angle {types}")

    for types in used["dihedrals"]:
        value = source.dihedral_types.get(types) or source.dihedral_types.get(
            tuple(reversed(types))
        )
        rows = _lookup(bonded["dihedraltypes"], types)
        if value is None or not rows:
            failures.append(f"missing dihedral {types}")
            continue
        expected = sorted((term.phase, term.phi_k * 4.184, term.per) for term in value)
        actual = sorted(
            (float(row[5]), float(row[6]), int(row[7]))
            for row in rows
            if row[4] == "9"
        )
        if len(actual) != len(expected) or any(
            not (_close(a[0], e[0]) and _close(a[1], e[1]) and a[2] == e[2])
            for a, e in zip(actual, expected, strict=True)
        ):
            failures.append(f"dihedral {types}")

    if failures:
        raise ValueError(f"CHARMM36 coefficient audit failed ({len(failures)}): {failures[:12]}")
    return {
        "status": "passed",
        "atom_types": len(atom_types),
        "bond_types": len(used["bonds"]),
        "angle_types": len(used["angles"]),
        "proper_dihedral_types": len(used["dihedrals"]),
        "improper_types": 0,
        "conversions": {
            "distance": "angstrom to nm",
            "energy": "kcal/mol to kJ/mol",
            "harmonic_factor": 2,
            "lj": "CHARMM Rmin/2 and epsilon to GROMACS sigma and epsilon",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prm", type=Path, required=True)
    parser.add_argument("--itp", type=Path, required=True)
    parser.add_argument("--ffbonded", type=Path, required=True)
    parser.add_argument("--ffnonbonded", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.prm, args.itp, args.ffbonded, args.ffnonbonded)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
