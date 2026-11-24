"""Chain selection and GROMACS position-restraint support."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from cellulose_forge.geometry import BulkCrystal
from cellulose_forge.glycam_mapping import GlycamMappedSystem, glycam_to_crystal_atom_name


@dataclass(frozen=True)
class ChainRestraints:
    """Resolved cellulose-chain position restraints."""

    chain_ids: tuple[int, ...]
    force_constants: tuple[float, float, float]
    layer_ids: tuple[int, ...] = ()
    layer_mapping: tuple[tuple[int, tuple[int, ...]], ...] = ()
    atom_names: tuple[str, ...] = ()
    selection_mode: str = "all_heavy_atoms"


@dataclass(frozen=True)
class ResolvedRestraintAtoms:
    """One engine-neutral real-heavy-atom restraint selection."""

    atom_indices: tuple[int, ...]
    atoms_per_chain: tuple[tuple[int, int], ...]


_POLYMER_HEAVY_ATOM_NAMES = {
    "cellulose": frozenset({*(f"C{i}" for i in range(1, 7)), *(f"O{i}" for i in range(1, 7))}),
    "mannan": frozenset({*(f"C{i}" for i in range(1, 7)), *(f"O{i}" for i in range(1, 7))}),
    "chitin": frozenset(
        {
            *(f"C{i}" for i in range(1, 9)),
            *(f"O{i}" for i in range(1, 8)),
            "N2",
        }
    ),
}

_CHARMM_CHITIN_CANONICAL_NAMES = {
    "N": "N2",
    "C": "C7",
    "O": "O7",
    "CT": "C8",
    "HN": "HN2",
    "HT1": "H81",
    "HT2": "H82",
    "HT3": "H83",
}


def canonical_atom_name(name: str, *, polymer: str, naming: str = "canonical") -> str:
    """Translate one backend atom name to the mapped-template namespace."""
    normalized = name.upper()
    if naming == "charmm" and polymer == "chitin":
        return _CHARMM_CHITIN_CANONICAL_NAMES.get(normalized, normalized)
    if naming == "glycam":
        return glycam_to_crystal_atom_name(normalized)
    if naming not in {"canonical", "charmm", "glycam"}:
        raise ValueError(f"unknown restraint atom naming convention {naming!r}")
    return normalized


def atom_element(name: str) -> str:
    """Return the element encoded by a canonical mapped-template atom name."""
    element = next((character for character in name.upper() if character.isalpha()), "")
    if element not in {"C", "H", "N", "O", "P", "S"}:
        raise ValueError(f"cannot determine element for restraint atom {name!r}")
    return element


def resolve_restraint_atom_indices(  # noqa: PLR0913
    *,
    canonical_atom_names: tuple[str, ...],
    elements: tuple[str, ...],
    site_kinds: tuple[str, ...],
    chain_numbers: tuple[int, ...],
    restrained_chain_ids: tuple[int, ...],
    requested_atom_names: tuple[str, ...] = (),
) -> ResolvedRestraintAtoms:
    """Resolve canonical names to real-heavy-particle indices exactly once."""
    size = len(canonical_atom_names)
    if any(len(values) != size for values in (elements, site_kinds, chain_numbers)):
        raise ValueError("restraint atom metadata arrays must have equal lengths")
    selected_chains = set(restrained_chain_ids)
    requested = {name.upper() for name in requested_atom_names}
    counts = dict.fromkeys(restrained_chain_ids, 0)
    names_by_chain: dict[int, set[str]] = {chain_id: set() for chain_id in restrained_chain_ids}
    indices: list[int] = []
    for index, (name, element, kind, chain_id) in enumerate(
        zip(canonical_atom_names, elements, site_kinds, chain_numbers, strict=True)
    ):
        if chain_id not in selected_chains or kind != "atom" or element.upper() == "H":
            continue
        if requested and name.upper() not in requested:
            continue
        indices.append(index)
        counts[chain_id] += 1
        names_by_chain[chain_id].add(name.upper())
    missing_by_chain = {
        chain_id: sorted(requested.difference(names))
        for chain_id, names in names_by_chain.items()
        if requested.difference(names)
    }
    if missing_by_chain:
        detail = "; ".join(
            f"chain {chain_id}: {', '.join(names)}"
            for chain_id, names in sorted(missing_by_chain.items())
        )
        raise ValueError(
            "position restraint atom names are absent from the written topology: " + detail
        )
    empty_chains = [chain_id for chain_id, count in counts.items() if count == 0]
    if empty_chains:
        raise ValueError(
            "position restraints selected no real heavy atoms for chain IDs: "
            + ", ".join(map(str, empty_chains))
        )
    return ResolvedRestraintAtoms(tuple(indices), tuple(sorted(counts.items())))


def resolve_mapped_restraint_atoms(
    system: GlycamMappedSystem,
    restrained_chain_ids: tuple[int, ...],
    *,
    requested_atom_names: tuple[str, ...] = (),
    polymer: str = "cellulose",
    naming: str = "canonical",
) -> ResolvedRestraintAtoms:
    """Resolve real-heavy atoms in a mapped all-atom system."""
    canonical = tuple(
        canonical_atom_name(name, polymer=polymer, naming=naming) for name in system.atom_names
    )
    return resolve_restraint_atom_indices(
        canonical_atom_names=canonical,
        elements=tuple(atom_element(name) for name in canonical),
        site_kinds=("atom",) * len(canonical),
        chain_numbers=system.chain_numbers,
        restrained_chain_ids=restrained_chain_ids,
        requested_atom_names=requested_atom_names,
    )


def transverse_chain_layers(crystal: BulkCrystal) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Group chains into deterministic transverse molecular layers.

    Slabs and fibrils carry explicit layer membership assigned by their geometry
    builders. The coordinate-based fallback groups chains along Cartesian z.
    """
    chain_numbers = np.asarray(crystal.chain_numbers)
    coordinates = np.asarray(crystal.cartesian_coordinates)
    if crystal.layer_numbers:
        return tuple(
            (
                layer_id,
                tuple(
                    sorted(
                        {
                            chain_id
                            for chain_id, atom_layer in zip(
                                crystal.chain_numbers, crystal.layer_numbers, strict=True
                            )
                            if atom_layer == layer_id
                        }
                    )
                ),
            )
            for layer_id in sorted(set(crystal.layer_numbers))
        )
    centers = [
        (chain_id, float(coordinates[chain_numbers == chain_id, 2].mean()))
        for chain_id in sorted(set(crystal.chain_numbers))
    ]
    grouped: list[list[tuple[int, float]]] = []
    for chain_id, center in sorted(centers, key=lambda item: (item[1], item[0])):
        if not grouped or not np.isclose(center, grouped[-1][0][1], atol=1e-5, rtol=0.0):
            grouped.append([])
        grouped[-1].append((chain_id, center))
    return tuple(
        (layer_id, tuple(sorted(chain_id for chain_id, _ in members)))
        for layer_id, members in enumerate(grouped, start=1)
    )


def resolve_chain_restraints(  # noqa: PLR0913
    crystal: BulkCrystal,
    *,
    chain_ids: list[int] | None,
    layer_ids: list[int] | None = None,
    force_constants: list[float] | None = None,
    atom_names: list[str] | None = None,
    polymer: str = "cellulose",
) -> ChainRestraints | None:
    """Resolve universal selectors to deterministic chain IDs."""
    if chain_ids and layer_ids:
        raise ValueError("--restrain-chains and --restrain-layers are mutually exclusive")
    if not chain_ids and not layer_ids:
        if force_constants is not None or atom_names:
            raise ValueError(
                "--restraint-force and --restrain-atom-names require a chain-restraint selector"
            )
        return None
    if force_constants is None:
        raise ValueError("--restraint-force FX FY FZ is required with chain restraints")
    if not all(np.isfinite(value) for value in force_constants):
        raise ValueError("position-restraint force constants must be finite")
    if any(value < 0 for value in force_constants):
        raise ValueError("position-restraint force constants must be nonnegative")

    if polymer not in _POLYMER_HEAVY_ATOM_NAMES:
        raise ValueError(f"position restraints are not defined for polymer {polymer!r}")
    requested_names = tuple(dict.fromkeys(name.upper() for name in (atom_names or [])))
    unknown_names = set(requested_names).difference(_POLYMER_HEAVY_ATOM_NAMES[polymer])
    if unknown_names:
        raise ValueError(
            f"unknown {polymer} restraint atom names: " + ", ".join(sorted(unknown_names))
        )

    available = set(crystal.chain_numbers)
    selected = set(chain_ids or [])
    unknown = selected.difference(available)
    if unknown:
        raise ValueError(f"unknown chain IDs: {', '.join(map(str, sorted(unknown)))}")

    layer_mapping = transverse_chain_layers(crystal)
    layers_by_id = dict(layer_mapping)
    selected_layers = set(layer_ids or [])
    unknown_layers = selected_layers.difference(layers_by_id)
    if unknown_layers:
        raise ValueError(
            f"unknown layer IDs: {', '.join(map(str, sorted(unknown_layers)))}; "
            f"available layers are 1-{len(layer_mapping)}"
        )
    for layer_id in selected_layers:
        selected.update(layers_by_id[layer_id])

    return ChainRestraints(
        chain_ids=tuple(sorted(selected)),
        force_constants=tuple(force_constants),  # type: ignore[arg-type]
        layer_ids=tuple(sorted(selected_layers)),
        layer_mapping=layer_mapping,
        atom_names=requested_names,
        selection_mode="named_heavy_atoms" if requested_names else "all_heavy_atoms",
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
    atom_names: tuple[str, ...] | list[str],
    force_constants: tuple[float, float, float],
    define: str = "POSRES_CELLULOSE",
    *,
    canonical_atom_names: tuple[str, ...] | list[str] | None = None,
    selected_atom_names: tuple[str, ...] = (),
) -> str:
    """Render conditional GROMACS position restraints for chain heavy atoms."""
    fx, fy, fz = force_constants
    lines = [
        f"#ifdef {define}",
        "[ position_restraints ]",
        "; atom  funct       fx       fy       fz",
    ]
    canonical = tuple(canonical_atom_names or atom_names)
    if len(canonical) != len(atom_names):
        raise ValueError("canonical restraint names must match the molecule atom count")
    selected_indices = set(
        resolve_restraint_atom_indices(
            canonical_atom_names=canonical,
            elements=tuple(atom_element(name) for name in canonical),
            site_kinds=("atom",) * len(canonical),
            chain_numbers=(1,) * len(canonical),
            restrained_chain_ids=(1,),
            requested_atom_names=selected_atom_names,
        ).atom_indices
    )
    lines.extend(
        f"{atom:6d}      1 {fx:10.3f} {fy:10.3f} {fz:10.3f}"
        for atom in range(1, len(atom_names) + 1)
        if atom - 1 in selected_indices
    )
    lines.extend(("#endif", ""))
    return "\n".join(lines)
