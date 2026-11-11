"""Convert published chitin CIF files into CelluloseForge crystal definitions.

The CIF files are the supporting information of the cited papers and are not
redistributed with CelluloseForge. This script records exactly how the package
data were derived from them:

1. Read the cell, symmetry operations, and atom sites, keeping one
   conformation of any disordered group.
2. Name every atom from its covalent connectivity (C1-C8, N2, O3-O7 and
   their hydrogens), independently of the labels used in the CIF.
3. Represent the glycosidic oxygen as O4 of the residue whose C4 it bonds,
   matching the cellulose definitions, by taking the symmetry image of the
   glycosidic oxygen that lies next to C4.
4. Express the space group as a two-fold screw along the chain axis acting on
   one residue per chain. For P2_1 2_1 2_1 the second, antiparallel chain is
   the image of the first under the operation (-x, y+1/2, -z+1/2).

Usage::

    python tools/convert_chitin_cif.py alpha path/to/300K.cif --disorder-group 2
    python tools/convert_chitin_cif.py alpha path/to/300K.cif --disorder-group 1 \
        --name alpha_o6_minor
    python tools/convert_chitin_cif.py beta path/to/ma102240r_si_001.cif
"""

from __future__ import annotations

import argparse
import json
import re
import shlex
from fractions import Fraction
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1] / "src/cellulose_forge/data/crystallography/chitin"
ATOM_ORDER = (
    "C1", "H1", "C2", "H2", "C3", "H3", "C4", "H4", "C5", "H5", "C6", "H61", "H62",
    "N2", "HN2", "C7", "O7", "C8", "H81", "H82", "H83",
    "O3", "O4", "O5", "O6", "HO3", "HO6",
)
SOURCES = {
    "alpha": {
        "allomorph": "alpha",
        "space_group": "P2_12_12_1",
        "citation": (
            "Sikorski, P.; Hori, R.; Wada, M. Biomacromolecules 2009, 10, 1100-1105."
        ),
        "doi": "10.1021/bm801251e",
        "file": "300K.cif (supporting information)",
        "chain_generator": "-x+1/2, -y, z+1/2",
        "second_chain": "-x, y+1/2, -z+1/2",
    },
    "beta": {
        "allomorph": "beta",
        "space_group": "P2_1",
        "citation": (
            "Nishiyama, Y.; Noishiki, Y.; Wada, M. Macromolecules 2011, 44, 950-957."
        ),
        "doi": "10.1021/ma102240r",
        "file": "ma102240r_si_001.cif (supporting information)",
        "chain_generator": "-x, -y, z+1/2",
        "second_chain": None,
    },
}


def _number(text: str) -> float:
    return float(re.sub(r"\(\d+\)$", "", text))


def _read_cif(path: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    lines = path.read_bytes().decode("latin-1").splitlines()
    scalars: dict[str, str] = {}
    atoms: list[dict[str, str]] = []
    index = 0
    while index < len(lines):
        line = lines[index].strip()
        if line == "loop_":
            columns: list[str] = []
            index += 1
            while index < len(lines) and lines[index].strip().startswith("_"):
                columns.append(lines[index].strip())
                index += 1
            rows: list[dict[str, str]] = []
            while index < len(lines):
                text = lines[index].strip()
                if not text or text.startswith(("_", "loop_", "#", ";")):
                    break
                values = shlex.split(text)
                if len(values) == len(columns):
                    rows.append(dict(zip(columns, values, strict=True)))
                elif len(values) == 1 and len(columns) == 1:
                    rows.append({columns[0]: values[0]})
                index += 1
            if "_atom_site_fract_x" in columns:
                atoms = rows
            if columns == ["_symmetry_equiv_pos_as_xyz"]:
                scalars["symmetry"] = "|".join(row[columns[0]] for row in rows)
            continue
        if line.startswith("_") and " " in line:
            key, value = line.split(None, 1)
            scalars[key] = value.strip()
        index += 1
    return scalars, atoms


def _operation(text: str) -> tuple[np.ndarray, np.ndarray]:
    rotation = np.zeros((3, 3))
    translation = np.zeros(3)
    for row, term in enumerate(text.replace(" ", "").split(",")):
        for sign, symbol in re.findall(r"([+-]?)([xyz]|\d+/\d+|\d+)", term):
            factor = -1.0 if sign == "-" else 1.0
            if symbol in "xyz":
                rotation[row, "xyz".index(symbol)] = factor
            else:
                translation[row] += factor * float(Fraction(symbol))
    return rotation, translation


def _cell_matrix(cell: list[float]) -> np.ndarray:
    a, b, c, alpha, beta, gamma = cell
    al, be, ga = np.radians([alpha, beta, gamma])
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    cz = np.sqrt(c * c - cx * cx - cy * cy)
    return np.asarray([[a, 0.0, 0.0], [b * np.cos(ga), b * np.sin(ga), 0.0], [cx, cy, cz]])


def _name_atoms(elements: list[str], cartesian: np.ndarray) -> list[str]:
    """Assign residue atom names from covalent connectivity."""
    count = len(elements)
    neighbors: list[set[int]] = [set() for _ in range(count)]
    for first in range(count):
        for second in range(first + 1, count):
            limit = 1.25 if "H" in (elements[first], elements[second]) else 1.75
            if np.linalg.norm(cartesian[first] - cartesian[second]) < limit:
                neighbors[first].add(second)
                neighbors[second].add(first)
    names: dict[int, str] = {}

    def heavy(atom: int, element: str | None = None) -> list[int]:
        return [n for n in neighbors[atom] if elements[n] != "H"
                and (element is None or elements[n] == element)]

    nitrogen = next(i for i, e in enumerate(elements) if e == "N")
    names[nitrogen] = "N2"
    carbonyl = next(n for n in heavy(nitrogen, "C") if heavy(n, "O"))
    c2 = next(n for n in heavy(nitrogen, "C") if n != carbonyl)
    names[carbonyl], names[c2] = "C7", "C2"
    names[next(iter(heavy(carbonyl, "O")))] = "O7"
    names[next(n for n in heavy(carbonyl, "C"))] = "C8"
    ring_oxygen = next(
        i for i, e in enumerate(elements) if e == "O" and len(heavy(i, "C")) == 2  # noqa: PLR2004
    )
    names[ring_oxygen] = "O5"
    c1 = next(n for n in heavy(ring_oxygen, "C") if c2 in neighbors[n])
    c5 = next(n for n in heavy(ring_oxygen, "C") if n != c1)
    names[c1], names[c5] = "C1", "C5"
    c3 = next(n for n in heavy(c2, "C") if n != c1)
    names[c3] = "C3"
    c4 = next(n for n in heavy(c3, "C") if n != c2)
    names[c4] = "C4"
    names[next(n for n in heavy(c5, "C") if n != c4)] = "C6"
    c6 = next(i for i, name in names.items() if name == "C6")
    names[next(iter(heavy(c6, "O")))] = "O6"
    names[next(iter(heavy(c3, "O")))] = "O3"
    glycosidic = [n for atom in (c1, c4) for n in heavy(atom, "O") if n != ring_oxygen]
    if len(glycosidic) != 1:
        raise SystemExit(f"expected one glycosidic oxygen, found {len(glycosidic)}")
    names[glycosidic[0]] = "OGLY"
    hydrogen_names = {
        "C1": ["H1"], "C2": ["H2"], "C3": ["H3"], "C4": ["H4"], "C5": ["H5"],
        "C6": ["H61", "H62"], "N2": ["HN2"], "C8": ["H81", "H82", "H83"],
        "O3": ["HO3"], "O6": ["HO6"],
    }
    for atom, name in list(names.items()):
        hydrogens = sorted(n for n in neighbors[atom] if elements[n] == "H")
        expected = hydrogen_names.get(name, [])
        if len(hydrogens) != len(expected):
            raise SystemExit(f"{name}: expected {len(expected)} hydrogens, found {len(hydrogens)}")
        names.update(zip(hydrogens, expected, strict=True))
    if len(names) != count:
        missing = [elements[i] for i in range(count) if i not in names]
        raise SystemExit(f"unnamed atoms remain: {missing}")
    return [names[i] for i in range(count)]


def convert(allomorph: str, cif: Path, disorder_group: str | None) -> dict[str, object]:
    """Return the CelluloseForge definition derived from one CIF."""
    source = SOURCES[allomorph]
    scalars, rows = _read_cif(cif)
    cell = [_number(scalars[f"_cell_{key}"]) for key in (
        "length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"
    )]
    matrix = _cell_matrix(cell)
    # Keep every fully occupied site and the partially occupied sites of the
    # requested disorder group. Selecting by group alone would drop H5 of the
    # alpha-chitin CIF, which is labelled group 1 although its occupancy is 1.
    selected = [
        row for row in rows
        if _number(row.get("_atom_site_occupancy", "1")) >= 1.0 - 1e-6
        or row.get("_atom_site_disorder_group", ".") == disorder_group
    ]
    elements = [row["_atom_site_type_symbol"] for row in selected]
    fractional = np.asarray([
        [_number(row[f"_atom_site_fract_{axis}"]) for axis in "xyz"] for row in selected
    ])
    names = _name_atoms(elements, fractional @ matrix)
    if sorted(names) != sorted([*ATOM_ORDER[:22], "OGLY", *ATOM_ORDER[23:]]):
        raise SystemExit(f"unexpected atom set: {sorted(names)}")
    position = dict(zip(names, fractional, strict=True))

    rotation, translation = _operation(str(source["chain_generator"]))
    candidates = []
    for direction in (-1, 1):
        image = position["OGLY"]
        for _ in range(abs(direction)):
            image = rotation @ image + translation * direction if direction > 0 else (
                np.linalg.solve(rotation, image - translation)
            )
        for shift in range(-2, 3):
            for sx in (-1, 0, 1):
                for sy in (-1, 0, 1):
                    shifted = image + np.asarray([sx, sy, shift])
                    distance = float(np.linalg.norm((shifted - position["C4"]) @ matrix))
                    candidates.append((distance, shifted))
    distance, o4 = min(candidates, key=lambda item: item[0])
    if not 1.3 < distance < 1.6:  # noqa: PLR2004
        raise SystemExit(f"no symmetry image of the glycosidic oxygen bonds C4 ({distance:.2f} A)")
    position["O4"] = o4
    del position["OGLY"]
    residue = np.asarray([position[name] for name in ATOM_ORDER])

    residues = [residue]
    chain_labels = ["origin"]
    chain_directions = [1]
    if source["second_chain"] is not None:
        second_rotation, second_translation = _operation(str(source["second_chain"]))
        residues.append(residue @ second_rotation.T + second_translation)
        chain_labels.append("center")
        chain_directions.append(-1)

    screw_rotation, screw_translation = _operation(str(source["chain_generator"]))
    return {
        "polymer": "chitin",
        "allomorph": source["allomorph"],
        "space_group": source["space_group"],
        "chain_axis": "c",
        "provenance": {
            "citation": source["citation"],
            "doi": source["doi"],
            "transcribed_from": (
                f"{source['file']}; converted by tools/convert_chitin_cif.py"
                + (f" keeping disorder group {disorder_group}" if disorder_group else "")
            ),
            "temperature_k": _number(scalars.get("_diffrn_ambient_temperature", "nan")),
            "notes": (
                "Atom names were assigned from covalent connectivity. The glycosidic "
                "oxygen is stored as O4 of the residue whose C4 it bonds, using its "
                "symmetry image along the chain."
                + (
                    " The antiparallel second chain is the image of the first under "
                    f"({source['second_chain']}); together with the chain screw this "
                    "generates all four P2_12_12_1 positions."
                    if source["second_chain"] else ""
                )
            ),
        },
        "unit_cell": {
            "a": cell[0], "b": cell[1], "c": cell[2],
            "alpha_deg": cell[3], "beta_deg": cell[4], "gamma_deg": cell[5],
            "length_unit": "angstrom",
        },
        "symmetry_operations": [
            {"name": "identity", "rotation": np.eye(3, dtype=int).tolist(),
             "translation": [0.0, 0.0, 0.0]},
            {"name": "chain_screw", "description": (
                f"Two-fold screw along c, ({source['chain_generator']}), generating the "
                "second residue of every chain."),
             "rotation": screw_rotation.astype(int).tolist(),
             "translation": screw_translation.tolist()},
        ],
        "asymmetric_unit": {
            "description": (
                f"{len(residues)} N-acetyl-beta-D-glucosamine residue(s), one per chain, "
                f"{len(ATOM_ORDER)} atoms each."
            ),
            "residue_atom_count": len(ATOM_ORDER),
            "chain_labels": chain_labels,
            "chain_directions": chain_directions,
            "atom_labels": list(ATOM_ORDER) * len(residues),
            "fractional_coordinates": np.round(np.vstack(residues), 6).tolist(),
        },
    }


def main() -> int:
    """Convert one CIF and write the package data file."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("allomorph", choices=sorted(SOURCES))
    parser.add_argument("cif", type=Path)
    parser.add_argument("--disorder-group", help="CIF disorder group to keep")
    parser.add_argument("--name", help="data file stem (default: the allomorph)")
    args = parser.parse_args()
    definition = convert(args.allomorph, args.cif, args.disorder_group)
    DATA.mkdir(parents=True, exist_ok=True)
    output = DATA / f"{args.name or args.allomorph}.json"
    output.write_text(json.dumps(definition, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {output.relative_to(DATA.parents[3])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
