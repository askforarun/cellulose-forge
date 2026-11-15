"""Independently compare a generated NAMD bundle with its GROMACS topology."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

_CHARGE_TOLERANCE = 5.0e-7
_MASS_TOLERANCE = 5.0e-5


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("--output", type=Path)
    return parser


def _canonical(values: tuple[int, ...]) -> tuple[int, ...]:
    return min(values, values[::-1])


def _parse_psf(path: Path) -> dict[str, object]:
    lines = path.read_text(encoding="utf-8").splitlines()
    atoms: list[dict[str, object]] = []
    interactions: dict[str, list[tuple[int, ...]]] = {}
    arities = {"!NBOND": 2, "!NTHETA": 3, "!NPHI": 4, "!NIMPHI": 4}
    row = 0
    while row < len(lines):
        line = lines[row]
        if "!NATOM" in line:
            count = int(line.split()[0])
            for atom_line in lines[row + 1 : row + count + 1]:
                fields = atom_line.split()
                atoms.append(
                    {
                        "segment": fields[1],
                        "residue": int(fields[2]),
                        "residue_name": fields[3],
                        "name": fields[4],
                        "type": fields[5],
                        "charge": float(fields[6]),
                        "mass": float(fields[7]),
                    }
                )
            row += count + 1
            continue
        marker = next((key for key in arities if key in line), None)
        if marker is None:
            row += 1
            continue
        count = int(line.split()[0])
        arity = arities[marker]
        values: list[int] = []
        row += 1
        while len(values) < count * arity:
            values.extend(int(value) for value in lines[row].split())
            row += 1
        interactions[marker] = [
            tuple(value - 1 for value in values[start : start + arity])
            for start in range(0, len(values), arity)
        ]
    return {"atoms": atoms, **interactions}


def _checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compare(build: Path) -> dict[str, object]:
    """Return an independent PSF/GROMACS comparison for one build folder."""
    try:
        import parmed  # noqa: PLC0415 - optional validation dependency
    except ImportError as error:
        raise RuntimeError("compare_bundle.py requires ParmEd") from error

    psf = _parse_psf(build / "system.psf")
    atoms = psf["atoms"]
    assert isinstance(atoms, list)
    topology = parmed.gromacs.GromacsTopologyFile(
        str(build / "system.top"), xyz=str(build / "structure.gro")
    )
    mismatches: Counter[str] = Counter()
    for gromacs_atom, psf_atom in zip(topology.atoms, atoms, strict=True):
        assert isinstance(psf_atom, dict)
        mismatches["name"] += gromacs_atom.name != psf_atom["name"]
        mismatches["type"] += gromacs_atom.type != psf_atom["type"]
        mismatches["charge"] += abs(gromacs_atom.charge - psf_atom["charge"]) > _CHARGE_TOLERANCE
        mismatches["mass"] += abs(gromacs_atom.mass - psf_atom["mass"]) > _MASS_TOLERANCE

    gromacs_bonds = Counter(_canonical((bond.atom1.idx, bond.atom2.idx)) for bond in topology.bonds)
    gromacs_angles = Counter(
        _canonical((angle.atom1.idx, angle.atom2.idx, angle.atom3.idx)) for angle in topology.angles
    )
    gromacs_propers = Counter(
        _canonical((term.atom1.idx, term.atom2.idx, term.atom3.idx, term.atom4.idx))
        for term in topology.dihedrals
        if not term.improper
    )
    improper_dihedrals = Counter(
        (term.atom1.idx, term.atom2.idx, term.atom3.idx, term.atom4.idx)
        for term in topology.dihedrals
        if term.improper
    )
    explicit_impropers = Counter(
        (term.atom1.idx, term.atom2.idx, term.atom3.idx, term.atom4.idx)
        for term in topology.impropers
    )
    # ParmEd may expose the same ordered CHARMM improper through both
    # collections. Counter union keeps its true multiplicity without adding
    # the duplicate view; atom order and the central atom remain significant.
    gromacs_impropers = improper_dihedrals | explicit_impropers
    psf_bonds = Counter(_canonical(term) for term in psf["!NBOND"])
    psf_angles = Counter(_canonical(term) for term in psf["!NTHETA"])
    psf_propers = Counter(_canonical(term) for term in psf["!NPHI"])
    psf_impropers = Counter(psf["!NIMPHI"])

    configuration = (build / "run_namd.conf").read_text(encoding="utf-8").splitlines()
    vectors = np.asarray(
        [
            [float(value) for value in line.split()[1:]]
            for line in configuration
            if line.startswith("cellBasisVector")
        ]
    )
    volume = abs(float(np.linalg.det(vectors)))
    widths = [
        volume / float(np.linalg.norm(np.cross(vectors[(axis + 1) % 3], vectors[(axis + 2) % 3])))
        for axis in range(3)
    ]
    return {
        "build": build.name,
        "atom_count": len(atoms),
        "atom_row_mismatches": {key: value for key, value in mismatches.items() if value},
        "net_charge_psf": round(sum(float(atom["charge"]) for atom in atoms), 8),
        "sets_equal": {
            "bonds": gromacs_bonds == psf_bonds,
            "angles": gromacs_angles == psf_angles,
            "proper_dihedrals": gromacs_propers == psf_propers,
            "impropers": gromacs_impropers == psf_impropers,
        },
        "cell_widths_angstrom": widths,
        "parameter_sha256": _checksum(build / "par_all36_carb.prm"),
    }


def main() -> None:
    args = _parser().parse_args()
    result = compare(args.build)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
