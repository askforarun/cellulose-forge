"""Map crystallographic cellulose coordinates onto a finite GLYCAM topology."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cellulose_forge.crystallography import FloatArray
from cellulose_forge.geometry import BulkCrystal

_GLYCAM_TO_CRYSTAL_NAME = {"H2O": "HO2", "H3O": "HO3", "H6O": "HO6"}
_REDUCING_ALIGNMENT_ATOMS = ("C1", "H1", "O5")
_NONREDUCING_ALIGNMENT_ATOMS = ("O4", "C4", "H4")
_RING_ALIGNMENT_ATOMS = ("C1", "C3", "C5")
_MISSING_CELLULOSE_II_HYDROGENS = frozenset({"HO2", "HO3", "HO6"})
_HYDROXYL_OXYGEN_BY_GLYCAM_HYDROGEN = {"H2O": "O2", "H3O": "O3", "H6O": "O6"}
_MAXIMUM_GLYCOSIDIC_BOND_ANGSTROM = 2.0


@dataclass(frozen=True)
class GlycamMappedSystem:
    """Coordinates and identifiers ordered exactly like repeated GLYCAM chains."""

    atom_names: tuple[str, ...]
    residue_names: tuple[str, ...]
    residue_numbers: tuple[int, ...]
    chain_numbers: tuple[int, ...]
    cartesian_coordinates: FloatArray
    basis_vectors: FloatArray

    @property
    def atom_count(self) -> int:
        """Return the mapped system atom count."""
        return len(self.atom_names)

    @property
    def chain_count(self) -> int:
        """Return the number of mapped cellulose chains."""
        return len(set(self.chain_numbers))


@dataclass(frozen=True)
class _ReferenceAtom:
    name: str
    residue_name: str
    residue_number: int
    coordinates: FloatArray


def _read_reference_pdb(path: str | Path) -> tuple[_ReferenceAtom, ...]:
    atoms: list[_ReferenceAtom] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        fields = line.split()
        atoms.append(
            _ReferenceAtom(
                name=fields[2],
                residue_name=fields[3],
                residue_number=int(fields[4]),
                coordinates=np.asarray([float(value) for value in fields[5:8]]),
            )
        )
    if not atoms:
        msg = f"no atoms were found in GLYCAM reference PDB: {path}"
        raise ValueError(msg)
    return tuple(atoms)


def _atoms_by_residue(
    atoms: tuple[_ReferenceAtom, ...],
) -> dict[int, tuple[_ReferenceAtom, ...]]:
    residue_numbers = dict.fromkeys(atom.residue_number for atom in atoms)
    return {
        residue_number: tuple(atom for atom in atoms if atom.residue_number == residue_number)
        for residue_number in residue_numbers
    }


def _coordinate_lookup(
    crystal: BulkCrystal,
    *,
    chain_number: int,
    residue_number: int,
) -> dict[str, FloatArray]:
    return {
        label: coordinates
        for label, coordinates, atom_chain, atom_residue in zip(
            crystal.atom_labels,
            crystal.cartesian_coordinates,
            crystal.chain_numbers,
            crystal.residue_numbers,
            strict=True,
        )
        if atom_chain == chain_number and atom_residue == residue_number
    }


def _require_crystallographic_atoms(
    crystal: BulkCrystal,
    reference: tuple[_ReferenceAtom, ...],
    *,
    constructed_names: frozenset[str] = frozenset(),
) -> None:
    """Reject force-field mapping when required experimental atoms are absent."""
    required = {
        _GLYCAM_TO_CRYSTAL_NAME.get(atom.name, atom.name)
        for atom in reference
        if atom.name not in constructed_names
    }
    missing = sorted(required.difference(crystal.atom_labels))
    if missing:
        msg = "missing crystallographic atom coordinates required by GLYCAM: " + ", ".join(missing)
        raise ValueError(msg)


def _rigid_transform(
    reference: dict[str, FloatArray],
    target: dict[str, FloatArray],
    anchor_names: tuple[str, str, str],
) -> tuple[FloatArray, FloatArray, FloatArray]:
    reference_points = np.vstack([reference[name] for name in anchor_names])
    target_points = np.vstack([target[name] for name in anchor_names])
    reference_center = np.mean(reference_points, axis=0)
    target_center = np.mean(target_points, axis=0)
    covariance = (reference_points - reference_center).T @ (target_points - target_center)
    left, _, right_transpose = np.linalg.svd(covariance)
    rotation = left @ right_transpose
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right_transpose
    return rotation, reference_center, target_center


def _transform_coordinate(
    coordinate: FloatArray,
    transform: tuple[FloatArray, FloatArray, FloatArray],
) -> FloatArray:
    rotation, reference_center, target_center = transform
    return (coordinate - reference_center) @ rotation + target_center


def _finite_chain_residue_order(
    crystal: BulkCrystal,
    *,
    chain_number: int,
    degree_polymerization: int,
) -> tuple[int, ...]:
    """Order residues along an open crystallographic O4-to-C1 path."""
    residues = tuple(range(1, degree_polymerization + 1))
    lookups = {
        residue: _coordinate_lookup(
            crystal,
            chain_number=chain_number,
            residue_number=residue,
        )
        for residue in residues
    }
    successors: dict[int, int] = {}
    for residue in residues:
        candidates = sorted(
            (
                float(np.linalg.norm(lookups[residue]["O4"] - lookups[other]["C1"])),
                other,
            )
            for other in residues
            if other != residue
        )
        if candidates and candidates[0][0] <= _MAXIMUM_GLYCOSIDIC_BOND_ANGSTROM:
            successors[residue] = candidates[0][1]
    predecessors = set(successors.values())
    starts = [residue for residue in residues if residue not in predecessors]
    if len(starts) != 1:
        msg = f"chain {chain_number} does not form one finite glycosidic path"
        raise ValueError(msg)
    order = [starts[0]]
    while order[-1] in successors:
        following = successors[order[-1]]
        if following in order:
            msg = f"chain {chain_number} finite residue connectivity contains a cycle"
            raise ValueError(msg)
        order.append(following)
    if len(order) != degree_polymerization or set(order) != set(residues):
        msg = f"chain {chain_number} finite residue connectivity is disconnected"
        raise ValueError(msg)
    return tuple(order)


def _periodic_chain_residue_order(
    crystal: BulkCrystal,
    *,
    chain_number: int,
    degree_polymerization: int,
) -> tuple[int, ...]:
    """Order residues around a periodic crystallographic O4-to-C1 cycle."""
    residues = tuple(range(1, degree_polymerization + 1))
    lookups = {
        residue: _coordinate_lookup(
            crystal,
            chain_number=chain_number,
            residue_number=residue,
        )
        for residue in residues
    }
    inverse_box = np.linalg.inv(crystal.basis_vectors)

    def distance(first: FloatArray, second: FloatArray) -> float:
        fractional = (second - first) @ inverse_box
        fractional -= np.round(fractional)
        return float(np.linalg.norm(fractional @ crystal.basis_vectors))

    successors: dict[int, int] = {}
    for residue in residues:
        candidate, bond_length = min(
            (
                (other, distance(lookups[residue]["O4"], lookups[other]["C1"]))
                for other in residues
                if other != residue
            ),
            key=lambda item: item[1],
        )
        if bond_length > _MAXIMUM_GLYCOSIDIC_BOND_ANGSTROM:
            msg = (
                f"chain {chain_number} has no periodic O4-C1 linkage from crystal "
                f"residue {residue} (nearest distance {bond_length:.3f} angstrom)"
            )
            raise ValueError(msg)
        successors[residue] = candidate
    order = [1]
    while len(order) < degree_polymerization:
        following = successors[order[-1]]
        if following in order:
            msg = f"chain {chain_number} residue connectivity does not form one cycle"
            raise ValueError(msg)
        order.append(following)
    if successors[order[-1]] != order[0] or set(order) != set(residues):
        msg = f"chain {chain_number} residue connectivity is not a closed periodic cycle"
        raise ValueError(msg)
    return tuple(order)


def map_crystal_to_finite_glycam(  # noqa: PLR0915
    crystal: BulkCrystal,
    reference_pdb: str | Path,
) -> GlycamMappedSystem:
    """Order crystal atoms like GLYCAM and construct the three finite-end atoms."""
    if not crystal.chain_numbers or not crystal.residue_numbers:
        msg = "GLYCAM mapping requires explicit chain and residue numbers"
        raise ValueError(msg)
    reference = _read_reference_pdb(reference_pdb)
    generated_names = _MISSING_CELLULOSE_II_HYDROGENS if crystal.allomorph == "II" else frozenset()
    _require_crystallographic_atoms(
        crystal,
        reference,
        constructed_names=frozenset(
            {"O1", "HO1", "H4O"}
            | {
                name
                for name, crystal_name in _GLYCAM_TO_CRYSTAL_NAME.items()
                if crystal_name in generated_names
            }
        ),
    )
    reference_residues = _atoms_by_residue(reference)
    chain_numbers = tuple(dict.fromkeys(crystal.chain_numbers))
    degrees = {
        chain: len(
            {
                residue
                for atom_chain, residue in zip(
                    crystal.chain_numbers, crystal.residue_numbers, strict=True
                )
                if atom_chain == chain
            }
        )
        for chain in chain_numbers
    }
    if len(set(degrees.values())) != 1:
        msg = "all cellulose chains must have the same degree of polymerization"
        raise ValueError(msg)
    degree_polymerization = next(iter(degrees.values()))
    if len(reference_residues) != degree_polymerization + 1:
        msg = "GLYCAM reference PDB degree of polymerization does not match the crystal"
        raise ValueError(msg)

    output_names: list[str] = []
    output_residue_names: list[str] = []
    output_residue_numbers: list[int] = []
    output_chain_numbers: list[int] = []
    output_coordinates: list[FloatArray] = []
    reducing_reference = {atom.name: atom.coordinates for atom in reference_residues[2]}
    cap_reference = {atom.name: atom.coordinates for atom in reference_residues[1]}
    nonreducing_reference_number = degree_polymerization + 1
    nonreducing_reference = {
        atom.name: atom.coordinates for atom in reference_residues[nonreducing_reference_number]
    }

    for chain_number in chain_numbers:
        crystal_residue_order = _finite_chain_residue_order(
            crystal,
            chain_number=chain_number,
            degree_polymerization=degree_polymerization,
        )
        reducing_target = _coordinate_lookup(
            crystal,
            chain_number=chain_number,
            residue_number=crystal_residue_order[0],
        )
        nonreducing_target = _coordinate_lookup(
            crystal,
            chain_number=chain_number,
            residue_number=crystal_residue_order[-1],
        )
        reducing_transform = _rigid_transform(
            reducing_reference,
            reducing_target,
            _REDUCING_ALIGNMENT_ATOMS,
        )
        nonreducing_transform = _rigid_transform(
            nonreducing_reference,
            nonreducing_target,
            _NONREDUCING_ALIGNMENT_ATOMS,
        )

        for reference_residue_number, reference_atoms in reference_residues.items():
            if reference_residue_number == 1:
                for atom in reference_atoms:
                    output_names.append(atom.name)
                    output_residue_names.append(atom.residue_name)
                    output_residue_numbers.append(1)
                    output_chain_numbers.append(chain_number)
                    output_coordinates.append(
                        _transform_coordinate(cap_reference[atom.name], reducing_transform)
                    )
                continue

            crystal_residue_number = crystal_residue_order[reference_residue_number - 2]
            target = _coordinate_lookup(
                crystal,
                chain_number=chain_number,
                residue_number=crystal_residue_number,
            )
            reference_lookup = {atom.name: atom.coordinates for atom in reference_atoms}
            hydrogen_transform = (
                _rigid_transform(reference_lookup, target, _RING_ALIGNMENT_ATOMS)
                if generated_names
                else None
            )
            for atom in reference_atoms:
                output_names.append(atom.name)
                output_residue_names.append(atom.residue_name)
                output_residue_numbers.append(reference_residue_number)
                output_chain_numbers.append(chain_number)
                if atom.name == "H4O":
                    coordinate = _transform_coordinate(atom.coordinates, nonreducing_transform)
                else:
                    crystal_name = _GLYCAM_TO_CRYSTAL_NAME.get(atom.name, atom.name)
                    if crystal_name in generated_names:
                        assert hydrogen_transform is not None
                        oxygen_name = _HYDROXYL_OXYGEN_BY_GLYCAM_HYDROGEN[atom.name]
                        aligned_hydrogen = _transform_coordinate(
                            atom.coordinates,
                            hydrogen_transform,
                        )
                        aligned_oxygen = _transform_coordinate(
                            reference_lookup[oxygen_name],
                            hydrogen_transform,
                        )
                        direction = aligned_hydrogen - aligned_oxygen
                        direction /= np.linalg.norm(direction)
                        bond_length = np.linalg.norm(
                            atom.coordinates - reference_lookup[oxygen_name]
                        )
                        coordinate = target[oxygen_name] + direction * bond_length
                    else:
                        coordinate = target[crystal_name]
                output_coordinates.append(coordinate)

    coordinates = np.vstack(output_coordinates)
    coordinates.setflags(write=False)
    return GlycamMappedSystem(
        atom_names=tuple(output_names),
        residue_names=tuple(output_residue_names),
        residue_numbers=tuple(output_residue_numbers),
        chain_numbers=tuple(output_chain_numbers),
        cartesian_coordinates=coordinates,
        basis_vectors=crystal.basis_vectors,
    )


def map_crystal_to_periodic_glycam(
    crystal: BulkCrystal,
    reference_pdb: str | Path,
) -> GlycamMappedSystem:
    """Order periodic crystallographic atoms exactly like a cyclic GLYCAM chain."""
    if not crystal.chain_numbers or not crystal.residue_numbers:
        msg = "periodic GLYCAM mapping requires explicit chain and residue numbers"
        raise ValueError(msg)
    reference = _read_reference_pdb(reference_pdb)
    generated_names = _MISSING_CELLULOSE_II_HYDROGENS if crystal.allomorph == "II" else frozenset()
    _require_crystallographic_atoms(
        crystal,
        reference,
        constructed_names=frozenset(
            name
            for name, crystal_name in _GLYCAM_TO_CRYSTAL_NAME.items()
            if crystal_name in generated_names
        ),
    )
    reference_residues = _atoms_by_residue(reference)
    chain_numbers = tuple(dict.fromkeys(crystal.chain_numbers))
    degree_by_chain = {
        chain: len(
            {
                residue
                for atom_chain, residue in zip(
                    crystal.chain_numbers, crystal.residue_numbers, strict=True
                )
                if atom_chain == chain
            }
        )
        for chain in chain_numbers
    }
    if len(set(degree_by_chain.values())) != 1:
        msg = "all periodic cellulose chains must have the same degree of polymerization"
        raise ValueError(msg)
    degree_polymerization = next(iter(degree_by_chain.values()))
    if len(reference_residues) != degree_polymerization:
        msg = "periodic GLYCAM reference degree of polymerization does not match the crystal"
        raise ValueError(msg)
    if any(atom.residue_name != "4GB" for atoms in reference_residues.values() for atom in atoms):
        msg = "periodic GLYCAM reference must contain only internal 4GB residues"
        raise ValueError(msg)

    output_names: list[str] = []
    output_residue_names: list[str] = []
    output_residue_numbers: list[int] = []
    output_chain_numbers: list[int] = []
    output_coordinates: list[FloatArray] = []
    for chain_number in chain_numbers:
        crystal_residue_order = (
            _periodic_chain_residue_order(
                crystal,
                chain_number=chain_number,
                degree_polymerization=degree_polymerization,
            )
            if crystal.allomorph == "II"
            else tuple(range(degree_polymerization, 0, -1))
        )
        for reference_residue_number, reference_atoms in reference_residues.items():
            crystal_residue_number = crystal_residue_order[reference_residue_number - 1]
            target = _coordinate_lookup(
                crystal,
                chain_number=chain_number,
                residue_number=crystal_residue_number,
            )
            reference_lookup = {atom.name: atom.coordinates for atom in reference_atoms}
            hydrogen_transform = (
                _rigid_transform(reference_lookup, target, _RING_ALIGNMENT_ATOMS)
                if generated_names
                else None
            )
            for atom in reference_atoms:
                crystal_name = _GLYCAM_TO_CRYSTAL_NAME.get(atom.name, atom.name)
                output_names.append(atom.name)
                output_residue_names.append("4GB")
                output_residue_numbers.append(reference_residue_number)
                output_chain_numbers.append(chain_number)
                if crystal_name in generated_names:
                    assert hydrogen_transform is not None
                    oxygen_name = _HYDROXYL_OXYGEN_BY_GLYCAM_HYDROGEN[atom.name]
                    aligned_hydrogen = _transform_coordinate(atom.coordinates, hydrogen_transform)
                    aligned_oxygen = _transform_coordinate(
                        reference_lookup[oxygen_name], hydrogen_transform
                    )
                    direction = aligned_hydrogen - aligned_oxygen
                    direction /= np.linalg.norm(direction)
                    bond_length = np.linalg.norm(atom.coordinates - reference_lookup[oxygen_name])
                    output_coordinates.append(target[oxygen_name] + direction * bond_length)
                else:
                    output_coordinates.append(target[crystal_name])

    coordinates = np.vstack(output_coordinates)
    coordinates.setflags(write=False)
    return GlycamMappedSystem(
        atom_names=tuple(output_names),
        residue_names=tuple(output_residue_names),
        residue_numbers=tuple(output_residue_numbers),
        chain_numbers=tuple(output_chain_numbers),
        cartesian_coordinates=coordinates,
        basis_vectors=crystal.basis_vectors,
    )
