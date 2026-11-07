"""Periodic-image safety checks for force-field-mapped cellulose systems."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from cellulose_forge.crystallography import FloatArray

if TYPE_CHECKING:
    from cellulose_forge.glycam_mapping import GlycamMappedSystem

_BOX_DETERMINANT_TOLERANCE = 1.0e-12
_MAXIMUM_GLYCOSIDIC_BOND_ANGSTROM = 2.0
_GROMACS_ATOM_FIELD_COUNT = 5
_GROMACS_BOND_FIELD_COUNT = 2

_TopologyBlock = tuple[str, dict[int, tuple[int, str]], list[tuple[int, int]]]


@dataclass(frozen=True)
class PeriodicImageContact:
    """Closest contact between a finite chain and its longitudinal image."""

    distance_angstrom: float
    chain_id: int
    atom_index: int
    image_atom_index: int
    atom_name: str
    image_atom_name: str

    def as_dict(self) -> dict[str, float | int | str]:
        """Return JSON-compatible contact metadata."""
        return {
            "distance_angstrom": self.distance_angstrom,
            "chain_id": self.chain_id,
            "atom_index": self.atom_index,
            "image_atom_index": self.image_atom_index,
            "atom_name": self.atom_name,
            "image_atom_name": self.image_atom_name,
        }


def closest_longitudinal_periodic_image_contact(
    system: GlycamMappedSystem,
) -> PeriodicImageContact:
    """Find the closest same-chain contact across the chain-periodic boundary.

    Oriented slab and fibril outputs store the chain-periodic translation in
    the second box vector. Lateral components are reduced by the two remaining
    box vectors so the result remains valid for triclinic cells.
    """
    basis = np.asarray(system.basis_vectors, dtype=float)
    if basis.shape != (3, 3) or abs(np.linalg.det(basis)) < _BOX_DETERMINANT_TOLERANCE:
        msg = "periodic-image clash validation requires an invertible 3x3 box"
        raise ValueError(msg)
    inverse = np.linalg.inv(basis)
    coordinates = np.asarray(system.cartesian_coordinates, dtype=float)
    best: PeriodicImageContact | None = None
    chain_numbers = np.asarray(system.chain_numbers)
    for chain_id in dict.fromkeys(system.chain_numbers):
        global_indices = np.flatnonzero(chain_numbers == chain_id)
        chain_coordinates = coordinates[global_indices]
        for start in range(0, len(global_indices), 256):
            first = chain_coordinates[start : start + 256]
            delta = (
                chain_coordinates[np.newaxis, :, :]
                + basis[1]
                - first[:, np.newaxis, :]
            )
            fractional = delta @ inverse
            fractional[..., 0] -= np.rint(fractional[..., 0])
            fractional[..., 2] -= np.rint(fractional[..., 2])
            delta = fractional @ basis
            distances = np.linalg.norm(delta, axis=2)
            local_flat = int(np.argmin(distances))
            local_i, local_j = np.unravel_index(local_flat, distances.shape)
            distance = float(distances[local_i, local_j])
            atom_index = int(global_indices[start + local_i])
            image_atom_index = int(global_indices[local_j])
            candidate = PeriodicImageContact(
                distance_angstrom=distance,
                chain_id=int(chain_id),
                atom_index=atom_index + 1,
                image_atom_index=image_atom_index + 1,
                atom_name=system.atom_names[atom_index],
                image_atom_name=system.atom_names[image_atom_index],
            )
            if best is None or candidate.distance_angstrom < best.distance_angstrom:
                best = candidate
    if best is None:
        msg = "periodic-image clash validation requires at least one mapped chain"
        raise ValueError(msg)
    return best


def validate_finite_periodic_image_contact(
    system: GlycamMappedSystem,
) -> PeriodicImageContact:
    """Reject unphysical finite-chain contacts with longitudinal images."""
    contact = closest_longitudinal_periodic_image_contact(system)
    hydrogen_contact = contact.atom_name.startswith("H") or contact.image_atom_name.startswith(
        "H"
    )
    threshold = 1.0 if hydrogen_contact else 1.5
    if contact.distance_angstrom < threshold:
        msg = (
            "finite-chain periodic-image clash: "
            f"chain {contact.chain_id} atoms {contact.atom_index} "
            f"({contact.atom_name}) and image {contact.image_atom_index} "
            f"({contact.image_atom_name}) are {contact.distance_angstrom:.3f} angstrom "
            f"apart (minimum {threshold:.1f} angstrom)"
        )
        raise ValueError(msg)
    return contact


def validate_periodic_glycosidic_bonds(
    system: GlycamMappedSystem,
) -> dict[str, float | int]:
    """Validate every periodic O4--C1 linkage using minimum-image distances."""
    basis = np.asarray(system.basis_vectors, dtype=float)
    if basis.shape != (3, 3) or abs(np.linalg.det(basis)) < _BOX_DETERMINANT_TOLERANCE:
        msg = "periodic bond validation requires an invertible 3x3 box"
        raise ValueError(msg)
    inverse = np.linalg.inv(basis)
    coordinates = np.asarray(system.cartesian_coordinates, dtype=float)
    atom_lookup = {
        (chain, residue, name): index
        for index, (name, residue, chain) in enumerate(
            zip(
                system.atom_names,
                system.residue_numbers,
                system.chain_numbers,
                strict=True,
            )
        )
    }
    maximum = 0.0
    count = 0
    for chain_id in dict.fromkeys(system.chain_numbers):
        residue_order = tuple(
            dict.fromkeys(
                residue
                for residue, chain in zip(
                    system.residue_numbers,
                    system.chain_numbers,
                    strict=True,
                )
                if chain == chain_id
            )
        )
        if not residue_order:
            continue
        for residue, following in zip(
            residue_order,
            residue_order[1:] + residue_order[:1],
            strict=True,
        ):
            o4 = atom_lookup[(chain_id, residue, "O4")]
            c1 = atom_lookup[(chain_id, following, "C1")]
            displacement = coordinates[c1] - coordinates[o4]
            fractional = displacement @ inverse
            fractional -= np.rint(fractional)
            distance = float(np.linalg.norm(fractional @ basis))
            if distance > _MAXIMUM_GLYCOSIDIC_BOND_ANGSTROM:
                msg = (
                    "periodic glycosidic bond is too long: "
                    f"chain {chain_id}, residues {residue}->{following}, "
                    f"O4--C1 = {distance:.3f} angstrom "
                    f"(maximum {_MAXIMUM_GLYCOSIDIC_BOND_ANGSTROM:.1f})"
                )
                raise ValueError(msg)
            maximum = max(maximum, distance)
            count += 1
    if count == 0:
        msg = "periodic glycosidic bond validation found no cellulose linkages"
        raise ValueError(msg)
    return {"bond_count": count, "maximum_distance_angstrom": maximum}


def _read_gromacs_molecule_blocks(topology_path: str | Path) -> list[_TopologyBlock]:
    """Read local atoms and bonds from every molecule type in a topology."""
    blocks: list[_TopologyBlock] = []
    atoms: dict[int, tuple[int, str]] | None = None
    bonds: list[tuple[int, int]] = []
    molecule_name: str | None = None
    section = ""
    for raw_line in Path(topology_path).read_text(encoding="utf-8").splitlines():
        stripped = raw_line.split(";", 1)[0].strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("[") and stripped.endswith("]"):
            new_section = stripped[1:-1].strip().lower()
            if new_section == "moleculetype":
                if atoms is not None and molecule_name is not None:
                    blocks.append((molecule_name, atoms, bonds))
                atoms = {}
                bonds = []
                molecule_name = None
            section = new_section
            continue
        if atoms is None:
            continue
        fields = stripped.split()
        if section == "moleculetype" and molecule_name is None:
            molecule_name = fields[0]
        elif section == "atoms" and len(fields) >= _GROMACS_ATOM_FIELD_COUNT:
            atoms[int(fields[0])] = (int(fields[2]), fields[4])
        elif section == "bonds" and len(fields) >= _GROMACS_BOND_FIELD_COUNT:
            bonds.append((int(fields[0]), int(fields[1])))
    if atoms is not None and molecule_name is not None:
        blocks.append((molecule_name, atoms, bonds))
    if not blocks:
        msg = f"no [ moleculetype ] definition found in {topology_path}"
        raise ValueError(msg)
    return blocks


def _read_gromacs_molecule_counts(topology_path: str | Path) -> list[tuple[str, int]]:
    """Read the ordered molecule inventory from a complete system topology."""
    rows: list[tuple[str, int]] = []
    section = ""
    for raw_line in Path(topology_path).read_text(encoding="utf-8").splitlines():
        stripped = raw_line.split(";", 1)[0].strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped[1:-1].strip().lower()
            continue
        if section == "molecules":
            fields = stripped.split()
            rows.append((fields[0], int(fields[1])))
    if not rows:
        msg = f"no [ molecules ] inventory found in {topology_path}"
        raise ValueError(msg)
    return rows


def _validate_gromacs_molecule_block(
    block: _TopologyBlock,
    coordinates: FloatArray,
    basis: FloatArray,
    inverse: FloatArray,
) -> tuple[int, float]:
    """Validate actual inter-residue glycosidic bonds in one molecule type."""
    _, block_atoms, block_bonds = block
    if not block_atoms:
        return 0, 0.0
    local_count = max(block_atoms)
    if local_count > len(coordinates):
        msg = "GROMACS molecule topology contains more atoms than mapped coordinates"
        raise ValueError(msg)
    maximum = 0.0
    block_linkages = 0
    block_residues = {residue for residue, _ in block_atoms.values()}
    for first, second in block_bonds:
        if first not in block_atoms or second not in block_atoms:
            msg = "GROMACS [ bonds ] row references an undefined local atom"
            raise ValueError(msg)
        first_residue, first_name = block_atoms[first]
        second_residue, second_name = block_atoms[second]
        if first_residue == second_residue or {first_name, second_name} != {"O4", "C1"}:
            continue
        displacement = coordinates[second - 1] - coordinates[first - 1]
        fractional = displacement @ inverse
        fractional -= np.rint(fractional)
        distance = float(np.linalg.norm(fractional @ basis))
        if distance > _MAXIMUM_GLYCOSIDIC_BOND_ANGSTROM:
            msg = (
                "written GROMACS glycosidic bond is too long: "
                f"local atoms {first}--{second} ({first_name}--{second_name}) = "
                f"{distance:.3f} angstrom "
                f"(maximum {_MAXIMUM_GLYCOSIDIC_BOND_ANGSTROM:.1f})"
            )
            raise ValueError(msg)
        maximum = max(maximum, distance)
        block_linkages += 1
    if block_linkages != len(block_residues):
        msg = (
            "written periodic GROMACS molecule has an unexpected number of "
            f"inter-residue O4--C1 bonds: {block_linkages} for "
            f"{len(block_residues)} residues"
        )
        raise ValueError(msg)
    return block_linkages, maximum


def validate_periodic_gromacs_topology_bonds(
    system: GlycamMappedSystem,
    molecule_topology_path: str | Path,
    system_topology_path: str | Path | None = None,
) -> dict[str, float | int]:
    """Validate every instantiated O4--C1 bond in written GROMACS topology."""
    blocks = _read_gromacs_molecule_blocks(molecule_topology_path)
    definitions = {block[0]: block for block in blocks}
    inventory_path = system_topology_path or molecule_topology_path
    molecule_counts = _read_gromacs_molecule_counts(inventory_path)

    basis = np.asarray(system.basis_vectors, dtype=float)
    if basis.shape != (3, 3) or abs(np.linalg.det(basis)) < _BOX_DETERMINANT_TOLERANCE:
        msg = "topology bond validation requires an invertible 3x3 box"
        raise ValueError(msg)
    inverse = np.linalg.inv(basis)
    coordinates = np.asarray(system.cartesian_coordinates, dtype=float)
    maximum = 0.0
    linkage_count = 0
    validated_types: set[str] = set()
    molecule_count = 0
    coordinate_offset = 0
    for molecule_name, copies in molecule_counts:
        try:
            block = definitions[molecule_name]
        except KeyError as error:
            msg = f"[ molecules ] references undefined molecule type {molecule_name}"
            raise ValueError(msg) from error
        _, block_atoms, _ = block
        local_count = max(block_atoms)
        for _ in range(copies):
            molecule_coordinates = coordinates[coordinate_offset : coordinate_offset + local_count]
            count, block_maximum = _validate_gromacs_molecule_block(
                block, molecule_coordinates, basis, inverse
            )
            linkage_count += count
            maximum = max(maximum, block_maximum)
            coordinate_offset += local_count
            molecule_count += 1
        validated_types.add(molecule_name)
    if coordinate_offset != len(coordinates):
        msg = (
            "GROMACS [ molecules ] inventory covers "
            f"{coordinate_offset} atoms but mapped coordinates contain {len(coordinates)}"
        )
        raise ValueError(msg)
    return {
        "molecule_type_count": len(validated_types),
        "molecule_count": molecule_count,
        "bond_count": linkage_count,
        "maximum_distance_angstrom": maximum,
    }
