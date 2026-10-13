"""Chain selection and GROMACS position-restraint support."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from cellulose_forge.geometry import BulkCrystal
from cellulose_forge.glycam_mapping import GlycamMappedSystem


@dataclass(frozen=True)
class ChainRestraints:
    """Resolved cellulose-chain position restraints."""

    chain_ids: tuple[int, ...]
    force_constants: tuple[float, float, float]


def resolve_chain_restraints(
    crystal: BulkCrystal,
    *,
    chain_ids: list[int] | None,
    force_constants: list[float] | None,
) -> ChainRestraints | None:
    """Resolve universal selectors to deterministic chain IDs."""
    if not chain_ids:
        if force_constants is not None:
            raise ValueError("--restraint-force requires a chain-restraint selector")
        return None
    if force_constants is None:
        raise ValueError("--restraint-force FX FY FZ is required with chain restraints")
    if any(value < 0 for value in force_constants):
        raise ValueError("position-restraint force constants must be nonnegative")

    available = set(crystal.chain_numbers)
    selected = set(chain_ids or [])
    unknown = selected.difference(available)
    if unknown:
        raise ValueError(f"unknown chain IDs: {', '.join(map(str, sorted(unknown)))}")

    return ChainRestraints(
        chain_ids=tuple(sorted(selected)),
        force_constants=tuple(force_constants),  # type: ignore[arg-type]
    )


def group_mapped_chains(
    system: GlycamMappedSystem,
    restrained_chain_ids: tuple[int, ...],
    *,
    free_molecule_type: str = "CELLULOSE_CHAIN",
    restrained_molecule_type: str = "CELLULOSE_CHAIN_RESTRAINED",
) -> tuple[GlycamMappedSystem, list[dict[str, int | str]]]:
    """Place free chains before restrained chains and report new molecule indices."""
    restrained = set(restrained_chain_ids)
    chain_order = tuple(dict.fromkeys(system.chain_numbers))
    ordered = tuple(chain for chain in chain_order if chain not in restrained) + tuple(
        chain for chain in chain_order if chain in restrained
    )
    atom_indices = [
        index
        for chain in ordered
        for index, atom_chain in enumerate(system.chain_numbers)
        if atom_chain == chain
    ]
    reordered = GlycamMappedSystem(
        atom_names=tuple(system.atom_names[index] for index in atom_indices),
        residue_names=tuple(system.residue_names[index] for index in atom_indices),
        residue_numbers=tuple(system.residue_numbers[index] for index in atom_indices),
        chain_numbers=tuple(system.chain_numbers[index] for index in atom_indices),
        cartesian_coordinates=np.asarray(system.cartesian_coordinates)[atom_indices],
        basis_vectors=system.basis_vectors,
    )
    mapping: list[dict[str, int | str]] = []
    first_atom = 1
    for molecule_index, chain in enumerate(ordered, start=1):
        atom_count = sum(value == chain for value in reordered.chain_numbers)
        is_restrained = chain in restrained
        mapping.append(
            {
                "original_chain_id": chain,
                "gromacs_molecule_index": molecule_index,
                "molecule_type": restrained_molecule_type if is_restrained else free_molecule_type,
                "global_atom_start": first_atom,
                "global_atom_end": first_atom + atom_count - 1,
                "restraint_status": "restrained" if is_restrained else "free",
            }
        )
        first_atom += atom_count
    return reordered, mapping


def position_restraints_section(
    atom_count: int,
    force_constants: tuple[float, float, float],
) -> str:
    """Render conditional GROMACS position restraints for every atom in a chain."""
    fx, fy, fz = force_constants
    lines = [
        "#ifdef POSRES_CELLULOSE",
        "[ position_restraints ]",
        "; atom  funct       fx       fy       fz",
    ]
    lines.extend(
        f"{atom:6d}      1 {fx:10.3f} {fy:10.3f} {fz:10.3f}"
        for atom in range(1, atom_count + 1)
    )
    lines.extend(("#endif", ""))
    return "\n".join(lines)
