"""Prepare GROMACS explicit pair overrides for OpenMM-based readers.

OpenMM's :class:`GromacsTopFile` accepts ordinary ``[ pairs ]`` function-1
rows, but not the function-2 rows used by CelluloseForge to preserve mixed
GLYCAM 1-4 scaling.  This module makes a temporary parser-compatible copy and
expands every local pair override over the molecules in ``[ molecules ]``.
The caller must then install the recorded parameters in the target engine.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

_EXPLICIT_PAIR_FIELD_COUNT = 8

@dataclass(frozen=True)
class ExplicitPairOverride:
    """One global, zero-based nonbonded exception."""

    first: int
    second: int
    charge_product_e2: float
    sigma_nm: float
    epsilon_kj_mol: float


def prepare_explicit_pairs(  # noqa: PLR0912
    source: str | Path,
    destination: str | Path,
) -> tuple[Path, tuple[ExplicitPairOverride, ...]]:
    """Write a function-1 copy and return expanded function-2 parameters."""
    source_path = Path(source)
    destination_path = Path(destination)
    lines = source_path.read_text(encoding="utf-8").splitlines()
    section = ""
    molecule_type: str | None = None
    awaiting_molecule_type = False
    atom_counts: dict[str, int] = {}
    local_overrides: dict[str, list[tuple[int, int, float, float, float]]] = {}
    molecules: list[tuple[str, int]] = []
    rewritten: list[str] = []

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped[1:-1].strip().lower()
            awaiting_molecule_type = section == "moleculetype"
            rewritten.append(line)
            continue
        if not stripped or stripped.startswith((";", "#")):
            rewritten.append(line)
            continue
        fields = stripped.split()
        if awaiting_molecule_type:
            molecule_type = fields[0]
            atom_counts.setdefault(molecule_type, 0)
            local_overrides.setdefault(molecule_type, [])
            awaiting_molecule_type = False
        elif section == "atoms" and molecule_type is not None:
            atom_counts[molecule_type] = max(atom_counts[molecule_type], int(fields[0]))
        elif (
            section == "pairs"
            and molecule_type is not None
            and len(fields) >= _EXPLICIT_PAIR_FIELD_COUNT
        ):
            if fields[2] == "2":
                first, second = int(fields[0]) - 1, int(fields[1]) - 1
                fudge_qq = float(fields[3])
                charge_product = fudge_qq * float(fields[4]) * float(fields[5])
                local_overrides[molecule_type].append(
                    (first, second, charge_product, float(fields[6]), float(fields[7]))
                )
                rewritten.append(f"{fields[0]:>7} {fields[1]:>7} 1")
                continue
        elif section == "molecules":
            molecules.append((fields[0], int(fields[1])))
        rewritten.append(line)

    if not any(local_overrides.values()):
        return source_path, ()
    if not molecules:
        raise ValueError("explicit-pair topology has no [ molecules ] entries")

    expanded: list[ExplicitPairOverride] = []
    offset = 0
    for name, count in molecules:
        if name not in atom_counts:
            raise ValueError(f"[ molecules ] references unknown molecule type {name}")
        atom_count = atom_counts[name]
        for _ in range(count):
            for first, second, charge, sigma, epsilon in local_overrides[name]:
                expanded.append(
                    ExplicitPairOverride(
                        offset + first,
                        offset + second,
                        charge,
                        sigma,
                        epsilon,
                    )
                )
            offset += atom_count

    destination_path.write_text("\n".join(rewritten) + "\n", encoding="utf-8")
    return destination_path, tuple(expanded)
