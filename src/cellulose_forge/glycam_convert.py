"""ParmEd-backed AMBER-to-GROMACS conversion with explicit mixed 1-4 pairs.

ParmEd refuses to write a GROMACS topology when an AMBER topology contains
more than one 1-4 scaling convention.  GLYCAM GlcNAc has exactly that case:
most proper torsions use unity scaling, while the amide H-N-C-O term uses
SCEE=1.2 and SCNB=2.0.  GROMACS pair function 2 can represent the exception
explicitly, so this module normalizes the temporary ParmEd structure for the
ordinary conversion and then replaces the affected pair rows with explicit
parameters from the original AMBER topology.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

KCAL_TO_KJ = 4.184
ANGSTROM_TO_NM = 0.1
_MINIMUM_PAIR_FIELDS = 3


@dataclass(frozen=True)
class _PairOverride:
    first: int
    second: int
    fudge_qq: float
    charge_first: float
    charge_second: float
    sigma_nm: float
    epsilon_kj: float


def _mixed_pair_overrides(structure: Any) -> dict[tuple[int, int], _PairOverride]:
    """Extract non-unity proper-dihedral pairs before normalizing ParmEd types."""
    overrides: dict[tuple[int, int], _PairOverride] = {}
    for dihedral in structure.dihedrals:
        kind = dihedral.type
        if dihedral.improper or dihedral.ignore_end or kind is None:
            continue
        scee, scnb = float(kind.scee), float(kind.scnb)
        if math.isclose(scee, 1.0) and math.isclose(scnb, 1.0):
            continue
        if scee <= 0.0 or scnb <= 0.0:
            raise ValueError("proper GLYCAM dihedral has non-positive 1-4 scaling")
        first, second = dihedral.atom1, dihedral.atom4
        first_index, second_index = sorted((first.idx + 1, second.idx + 1))
        indices: tuple[int, int] = (first_index, second_index)
        sigma_angstrom = (float(first.rmin) + float(second.rmin)) / (2.0 ** (1.0 / 6.0))
        epsilon_kcal = math.sqrt(float(first.epsilon) * float(second.epsilon)) / scnb
        override = _PairOverride(
            first=indices[0],
            second=indices[1],
            fudge_qq=1.0 / scee,
            charge_first=float(structure.atoms[indices[0] - 1].charge),
            charge_second=float(structure.atoms[indices[1] - 1].charge),
            sigma_nm=sigma_angstrom * ANGSTROM_TO_NM,
            epsilon_kj=epsilon_kcal * KCAL_TO_KJ,
        )
        previous = overrides.get(indices)
        if previous is not None and previous != override:
            raise ValueError(f"conflicting GLYCAM 1-4 parameters for pair {indices}")
        overrides[indices] = override
    return overrides


def _normalize_proper_scaling(structure: Any) -> None:
    """Set proper-dihedral scaling to unity in the temporary conversion object."""
    for dihedral in structure.dihedrals:
        if not dihedral.improper and dihedral.type is not None:
            dihedral.type.scee = 1.0
            dihedral.type.scnb = 1.0


def _write_explicit_pairs(path: Path, overrides: dict[tuple[int, int], _PairOverride]) -> None:
    """Replace ordinary pair rows with explicit GROMACS function-2 rows."""
    lines = path.read_text(encoding="utf-8").splitlines()
    section = ""
    occurrences: dict[tuple[int, int], int] = dict.fromkeys(overrides, 0)
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped[1:-1].strip().lower()
            continue
        if section != "pairs" or not stripped or stripped.startswith((";", "#")):
            continue
        fields = stripped.split()
        if len(fields) < _MINIMUM_PAIR_FIELDS:
            continue
        pair_indices = sorted((int(fields[0]), int(fields[1])))
        pair: tuple[int, int] = (pair_indices[0], pair_indices[1])
        override = overrides.get(pair)
        if override is None:
            continue
        occurrences[pair] += 1
        if occurrences[pair] > 1:
            raise ValueError(f"mixed-scaling pair appears more than once: {pair}")
        lines[index] = (
            f"{override.first:7d} {override.second:7d} 2 "
            f"{override.fudge_qq:.10f} {override.charge_first:.10f} "
            f"{override.charge_second:.10f} {override.sigma_nm:.10f} "
            f"{override.epsilon_kj:.10f}"
        )
    missing = sorted(pair for pair, count in occurrences.items() if count == 0)
    if missing:
        raise ValueError(f"converted GROMACS topology omitted mixed-scaling pairs: {missing}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _replace_molecule_name(path: Path, molecule_name: str) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    section = ""
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped[1:-1].strip().lower()
            continue
        if section in {"moleculetype", "molecules"} and stripped and not stripped.startswith(";"):
            fields = stripped.split()
            fields[0] = molecule_name
            lines[index] = "  ".join(fields)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def convert(
    prmtop: Path,
    inpcrd: Path,
    topology: Path,
    coordinates: Path,
    molecule_name: str,
) -> None:
    """Convert one AMBER topology while preserving every 1-4 scaling exception."""
    import parmed as pmd  # type: ignore[import-not-found]  # noqa: PLC0415

    structure = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    overrides = _mixed_pair_overrides(structure)
    structure.save(str(coordinates), overwrite=True)
    _normalize_proper_scaling(structure)
    structure.save(str(topology), overwrite=True)
    _write_explicit_pairs(topology, overrides)
    _replace_molecule_name(topology, molecule_name)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prmtop", type=Path)
    parser.add_argument("inpcrd", type=Path)
    parser.add_argument("topology", type=Path)
    parser.add_argument("coordinates", type=Path)
    parser.add_argument("molecule_name")
    args = parser.parse_args()
    convert(args.prmtop, args.inpcrd, args.topology, args.coordinates, args.molecule_name)


if __name__ == "__main__":
    main()
