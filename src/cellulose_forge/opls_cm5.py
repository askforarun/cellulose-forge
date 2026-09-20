"""Experimental adapter for an externally supplied pristine OPLS-CM5 topology."""

from __future__ import annotations

import shutil
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cellulose_forge.crystallography import FloatArray
from cellulose_forge.geometry import BulkCrystal
from cellulose_forge.glycam_mapping import GlycamMappedSystem
from cellulose_forge.gromacs import write_gro

_REFERENCE_DP = 20
_REFERENCE_ATOM_COUNT = 423
_CONSTRUCTED_TERMINAL_ATOMS = frozenset({"O1", "HO1", "HO4"})
_MINIMUM_TOPOLOGY_ATOM_FIELDS = 8
_MINIMUM_GRO_LINES = 3
_INTERACTION_ARITY = {"bonds": 2, "angles": 3, "dihedrals": 4, "pairs": 2}
_TEMPLATE_FIRST_RESIDUE = 5
_TEMPLATE_LAST_RESIDUE = 16
_INTERNAL_REFERENCE_RESIDUE = 10
_MINIMUM_PERIODIC_DP = 4
_RING_ALIGNMENT_ATOMS = ("C1", "C3", "C5")
_MISSING_CELLULOSE_II_HYDROGENS = frozenset({"HO2", "HO3", "HO6"})
_HYDROXYL_OXYGEN_BY_HYDROGEN = {"HO2": "O2", "HO3": "O3", "HO6": "O6"}


@dataclass(frozen=True)
class _TopologyAtom:
    atom_type: str
    atom_name: str
    charge: float
    mass: float


@dataclass(frozen=True)
class _CoordinateAtom:
    residue_number: int
    atom_name: str
    coordinates: FloatArray


@dataclass(frozen=True, order=True)
class _InteractionTemplate:
    section: str
    atoms: tuple[tuple[int, str], ...]
    parameters: tuple[str, ...]


@dataclass(frozen=True)
class PeriodicOplsTopology:
    """A cyclic cellulose-chain topology derived from internal DP20 terms."""

    degree_polymerization: int
    coordinate_atom_names: tuple[str, ...]
    reference_coordinates: FloatArray
    topology_atoms: tuple[_TopologyAtom, ...]
    interactions: tuple[_InteractionTemplate, ...]
    source_repeat_charge: float

    @property
    def atom_count(self) -> int:
        """Return atoms in one periodic chain."""
        return self.degree_polymerization * len(self.topology_atoms)

    @property
    def corrected_repeat_charge(self) -> float:
        """Return the repeat charge after symmetry-preserving correction."""
        correction = self.source_repeat_charge / len(self.topology_atoms)
        return sum(atom.charge - correction for atom in self.topology_atoms)

    @property
    def charge_correction_per_atom(self) -> float:
        """Return the uniform correction applied to each rounded source charge."""
        return self.source_repeat_charge / len(self.topology_atoms)


def _read_topology_atoms(path: str | Path) -> tuple[_TopologyAtom, ...]:
    atoms: list[_TopologyAtom] = []
    section = ""
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped[1:-1].strip()
            continue
        if section != "atoms" or not stripped or stripped.startswith(";"):
            continue
        fields = stripped.split()
        if len(fields) < _MINIMUM_TOPOLOGY_ATOM_FIELDS:
            msg = f"malformed [ atoms ] row in OPLS-CM5 topology: {line}"
            raise ValueError(msg)
        atoms.append(
            _TopologyAtom(
                atom_type=fields[1],
                atom_name=fields[4],
                charge=float(fields[6]),
                mass=float(fields[7]),
            )
        )
    if len(atoms) != _REFERENCE_ATOM_COUNT:
        msg = (
            "the experimental OPLS-CM5 adapter requires the supplied pristine "
            f"DP20 topology with {_REFERENCE_ATOM_COUNT} atoms; found {len(atoms)}"
        )
        raise ValueError(msg)
    total_charge = sum(atom.charge for atom in atoms)
    if not np.isclose(total_charge, 0.0, atol=1.0e-6):
        msg = f"OPLS-CM5 reference topology is not neutral (charge {total_charge:.8f})"
        raise ValueError(msg)
    return tuple(atoms)


def _read_topology_interactions(path: str | Path) -> dict[str, tuple[tuple[str, ...], ...]]:
    interactions: dict[str, list[tuple[str, ...]]] = defaultdict(list)
    section = ""
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped[1:-1].strip()
            continue
        if section not in _INTERACTION_ARITY or not stripped or stripped.startswith(";"):
            continue
        fields = tuple(stripped.split())
        arity = _INTERACTION_ARITY[section]
        if len(fields) <= arity:
            msg = f"malformed [{section}] row in OPLS-CM5 topology: {line}"
            raise ValueError(msg)
        interactions[section].append(fields)
    missing = sorted(set(_INTERACTION_ARITY).difference(interactions))
    if missing:
        msg = "OPLS-CM5 reference is missing interaction sections: " + ", ".join(missing)
        raise ValueError(msg)
    return {section: tuple(rows) for section, rows in interactions.items()}


def _read_first_reference_chain(path: str | Path) -> tuple[_CoordinateAtom, ...]:
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if len(lines) < _MINIMUM_GRO_LINES:
        msg = f"malformed OPLS-CM5 reference GRO file: {path}"
        raise ValueError(msg)
    atom_count = int(lines[1].strip())
    if atom_count < _REFERENCE_ATOM_COUNT:
        msg = "OPLS-CM5 reference GRO does not contain one complete DP20 chain"
        raise ValueError(msg)
    atoms = tuple(
        _CoordinateAtom(
            residue_number=int(line[:5]),
            atom_name=line[10:15].strip(),
            coordinates=np.asarray(
                [float(line[20:28]), float(line[28:36]), float(line[36:44])]
            )
            / 0.1,
        )
        for line in lines[2 : 2 + _REFERENCE_ATOM_COUNT]
    )
    residue_counts = {
        residue: sum(atom.residue_number == residue for atom in atoms)
        for residue in dict.fromkeys(atom.residue_number for atom in atoms)
    }
    expected = {1: 22, **dict.fromkeys(range(2, 20), 21), 20: 23}
    if residue_counts != expected:
        msg = f"unexpected OPLS-CM5 DP20 residue layout: {residue_counts}"
        raise ValueError(msg)
    return atoms


def derive_periodic_opls_cm5_topology(
    topology_path: str | Path,
    reference_gro: str | Path,
    *,
    degree_polymerization: int,
) -> PeriodicOplsTopology:
    """Derive a cyclic chain from interactions wholly inside the DP20 reference."""
    if degree_polymerization < _MINIMUM_PERIODIC_DP or degree_polymerization % 2:
        msg = "periodic OPLS-CM5 chains require an even degree of polymerization of at least 4"
        raise ValueError(msg)
    source_atoms = _read_topology_atoms(topology_path)
    coordinates = _read_first_reference_chain(reference_gro)
    if len(source_atoms) != len(coordinates):
        msg = "OPLS-CM5 topology and coordinate references have different atom counts"
        raise ValueError(msg)

    internal_indices = tuple(
        index
        for index, atom in enumerate(coordinates)
        if atom.residue_number == _INTERNAL_REFERENCE_RESIDUE
    )
    coordinate_atom_names = tuple(coordinates[index].atom_name for index in internal_indices)
    if len(set(coordinate_atom_names)) != len(coordinate_atom_names):
        msg = "internal OPLS-CM5 repeat contains duplicate coordinate atom names"
        raise ValueError(msg)
    topology_atoms = tuple(source_atoms[index] for index in internal_indices)
    reference_coordinates = np.vstack(
        [coordinates[index].coordinates for index in internal_indices]
    )
    reference_coordinates.setflags(write=False)
    source_repeat_charge = sum(atom.charge for atom in topology_atoms)

    templates: set[_InteractionTemplate] = set()
    for section, rows in _read_topology_interactions(topology_path).items():
        arity = _INTERACTION_ARITY[section]
        for row in rows:
            indices = tuple(int(value) - 1 for value in row[:arity])
            keyed_atoms = tuple(
                (coordinates[index].residue_number, coordinates[index].atom_name)
                for index in indices
            )
            residues = tuple(residue for residue, _ in keyed_atoms)
            if (
                min(residues) < _TEMPLATE_FIRST_RESIDUE
                or max(residues) > _TEMPLATE_LAST_RESIDUE
            ):
                continue
            anchor = min(residues)
            templates.add(
                _InteractionTemplate(
                    section=section,
                    atoms=tuple((residue - anchor, name) for residue, name in keyed_atoms),
                    parameters=tuple(row[arity:]),
                )
            )

    expected_template_counts = {"bonds": 22, "angles": 41, "dihedrals": 66, "pairs": 64}
    actual_template_counts = {
        section: sum(template.section == section for template in templates)
        for section in _INTERACTION_ARITY
    }
    if actual_template_counts != expected_template_counts:
        msg = (
            "supplied OPLS-CM5 internal interaction pattern changed: "
            f"expected {expected_template_counts}, found {actual_template_counts}"
        )
        raise ValueError(msg)
    return PeriodicOplsTopology(
        degree_polymerization=degree_polymerization,
        coordinate_atom_names=coordinate_atom_names,
        reference_coordinates=reference_coordinates,
        topology_atoms=topology_atoms,
        interactions=tuple(sorted(templates)),
        source_repeat_charge=source_repeat_charge,
    )


def _coordinate_lookup(
    crystal: BulkCrystal, *, chain_number: int, residue_number: int
) -> dict[str, FloatArray]:
    return {
        label: coordinates
        for label, coordinates, chain, residue in zip(
            crystal.atom_labels,
            crystal.cartesian_coordinates,
            crystal.chain_numbers,
            crystal.residue_numbers,
            strict=True,
        )
        if chain == chain_number and residue == residue_number
    }


def _rigid_transform(
    reference: dict[str, FloatArray],
    target: dict[str, FloatArray],
    anchors: tuple[str, str, str],
) -> tuple[FloatArray, FloatArray, FloatArray]:
    reference_points = np.vstack([reference[name] for name in anchors])
    target_points = np.vstack([target[name] for name in anchors])
    reference_center = np.mean(reference_points, axis=0)
    target_center = np.mean(target_points, axis=0)
    covariance = (reference_points - reference_center).T @ (target_points - target_center)
    left, _, right_transpose = np.linalg.svd(covariance)
    rotation = left @ right_transpose
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right_transpose
    return rotation, reference_center, target_center


def _transform(
    coordinate: FloatArray,
    transform: tuple[FloatArray, FloatArray, FloatArray],
) -> FloatArray:
    rotation, reference_center, target_center = transform
    return (coordinate - reference_center) @ rotation + target_center


def map_finite_dp20_opls_cm5(
    crystal: BulkCrystal,
    topology_path: str | Path,
    reference_gro: str | Path,
) -> GlycamMappedSystem:
    """Map a finite DP20 crystal onto the supplied OPLS-CM5 atom order."""
    topology_atoms = _read_topology_atoms(topology_path)
    reference_atoms = _read_first_reference_chain(reference_gro)
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
    if not degree_by_chain or set(degree_by_chain.values()) != {_REFERENCE_DP}:
        msg = "the supplied pristine OPLS-CM5 topology supports finite DP20 only"
        raise ValueError(msg)

    reference_by_residue = {
        residue: tuple(atom for atom in reference_atoms if atom.residue_number == residue)
        for residue in range(1, _REFERENCE_DP + 1)
    }
    output_names: list[str] = []
    output_residue_names: list[str] = []
    output_residue_numbers: list[int] = []
    output_chain_numbers: list[int] = []
    output_coordinates: list[FloatArray] = []

    for chain_number in chain_numbers:
        first_target = _coordinate_lookup(crystal, chain_number=chain_number, residue_number=1)
        last_target = _coordinate_lookup(
            crystal, chain_number=chain_number, residue_number=_REFERENCE_DP
        )
        first_reference = {
            atom.atom_name: atom.coordinates for atom in reference_by_residue[1]
        }
        last_reference = {
            atom.atom_name: atom.coordinates for atom in reference_by_residue[_REFERENCE_DP]
        }
        first_transform = _rigid_transform(
            first_reference, first_target, ("O4", "C4", "H4")
        )
        last_transform = _rigid_transform(
            last_reference, last_target, ("C1", "H1", "O5")
        )
        topology_offset = 0
        for residue_number in range(1, _REFERENCE_DP + 1):
            target = _coordinate_lookup(
                crystal,
                chain_number=chain_number,
                residue_number=residue_number,
            )
            for reference_atom in reference_by_residue[residue_number]:
                topology_atom = topology_atoms[topology_offset]
                topology_offset += 1
                if reference_atom.atom_name in _CONSTRUCTED_TERMINAL_ATOMS:
                    transform = first_transform if residue_number == 1 else last_transform
                    coordinate = _transform(reference_atom.coordinates, transform)
                else:
                    coordinate = target[reference_atom.atom_name]
                output_names.append(topology_atom.atom_name)
                output_residue_names.append("CNC")
                output_residue_numbers.append(
                    (chain_number - 1) * _REFERENCE_DP + residue_number
                )
                output_chain_numbers.append(chain_number)
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


def map_periodic_opls_cm5(
    crystal: BulkCrystal,
    topology: PeriodicOplsTopology,
) -> GlycamMappedSystem:
    """Order a periodic crystal exactly like the generated cyclic topology."""
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
    if not degree_by_chain or set(degree_by_chain.values()) != {
        topology.degree_polymerization
    }:
        msg = (
            "periodic OPLS-CM5 topology DP does not match the generated "
            f"crystal chains: topology DP {topology.degree_polymerization}, "
            f"crystal DPs {sorted(set(degree_by_chain.values()))}"
        )
        raise ValueError(msg)

    output_names: list[str] = []
    output_residue_names: list[str] = []
    output_residue_numbers: list[int] = []
    output_chain_numbers: list[int] = []
    output_coordinates: list[FloatArray] = []
    for chain_number in chain_numbers:
        for residue_number in range(1, topology.degree_polymerization + 1):
            crystal_residue_number = (
                topology.degree_polymerization - residue_number + 1
                if crystal.allomorph == "II" and chain_number % 2 == 0
                else residue_number
            )
            lookup = _coordinate_lookup(
                crystal,
                chain_number=chain_number,
                residue_number=crystal_residue_number,
            )
            allowed_missing = (
                _MISSING_CELLULOSE_II_HYDROGENS
                if crystal.allomorph == "II"
                else frozenset()
            )
            missing = set(topology.coordinate_atom_names).difference(lookup).difference(
                allowed_missing
            )
            if missing:
                msg = (
                    f"chain {chain_number} residue {residue_number} lacks OPLS-CM5 "
                    f"atoms: {', '.join(sorted(missing))}"
                )
                raise ValueError(msg)
            global_residue = (
                (chain_number - 1) * topology.degree_polymerization + residue_number
            )
            reference_lookup = dict(
                zip(
                    topology.coordinate_atom_names,
                    topology.reference_coordinates,
                    strict=True,
                )
            )
            hydrogen_transform = (
                _rigid_transform(reference_lookup, lookup, _RING_ALIGNMENT_ATOMS)
                if allowed_missing
                else None
            )
            for coordinate_name, reference_coordinate, topology_atom in zip(
                topology.coordinate_atom_names,
                topology.reference_coordinates,
                topology.topology_atoms,
                strict=True,
            ):
                output_names.append(topology_atom.atom_name)
                output_residue_names.append("CNC")
                output_residue_numbers.append(global_residue)
                output_chain_numbers.append(chain_number)
                if coordinate_name in allowed_missing:
                    assert hydrogen_transform is not None
                    oxygen_name = _HYDROXYL_OXYGEN_BY_HYDROGEN[coordinate_name]
                    aligned_hydrogen = _transform(
                        reference_coordinate, hydrogen_transform
                    )
                    aligned_oxygen = _transform(
                        reference_lookup[oxygen_name], hydrogen_transform
                    )
                    direction = aligned_hydrogen - aligned_oxygen
                    direction /= np.linalg.norm(direction)
                    bond_length = np.linalg.norm(
                        reference_coordinate - reference_lookup[oxygen_name]
                    )
                    output_coordinates.append(
                        lookup[oxygen_name] + direction * bond_length
                    )
                else:
                    output_coordinates.append(lookup[coordinate_name])

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


def _instantiate_periodic_interactions(
    topology: PeriodicOplsTopology,
) -> dict[str, tuple[tuple[str, ...], ...]]:
    atoms_per_residue = len(topology.topology_atoms)
    position_by_name = {
        name: position for position, name in enumerate(topology.coordinate_atom_names)
    }
    instantiated: dict[str, set[tuple[str, ...]]] = defaultdict(set)
    for anchor in range(topology.degree_polymerization):
        for template in topology.interactions:
            indices = tuple(
                str(
                    ((anchor + offset) % topology.degree_polymerization)
                    * atoms_per_residue
                    + position_by_name[name]
                    + 1
                )
                for offset, name in template.atoms
            )
            instantiated[template.section].add(indices + template.parameters)
    return {
        section: tuple(sorted(rows, key=lambda row: tuple(int(x) for x in row[:arity])))
        for section, arity in _INTERACTION_ARITY.items()
        if (rows := instantiated[section])
    }


def periodic_interaction_counts(
    topology: PeriodicOplsTopology,
) -> dict[str, int]:
    """Return generated bonded-interaction counts for one cyclic chain."""
    return {
        section: len(rows)
        for section, rows in _instantiate_periodic_interactions(topology).items()
    }


def render_periodic_opls_cm5_itp(topology: PeriodicOplsTopology) -> str:
    """Render a charge-neutral cyclic OPLS-CM5 chain as a GROMACS ITP."""
    lines = ["[ moleculetype ]", "; name  nrexcl", "CNC     3", "", "[ atoms ]"]
    lines.append("; nr type resnr residue atom cgnr charge mass")
    atom_number = 0
    for residue_number in range(1, topology.degree_polymerization + 1):
        for atom in topology.topology_atoms:
            atom_number += 1
            corrected_charge = atom.charge - topology.charge_correction_per_atom
            lines.append(
                f"{atom_number:6d} {atom.atom_type:<10s} {residue_number:6d} "
                f"CNC {atom.atom_name:<5s} {atom_number:6d} "
                f"{corrected_charge: .10f} {atom.mass:.6f}"
            )

    interactions = _instantiate_periodic_interactions(topology)
    for section in ("bonds", "angles", "dihedrals", "pairs"):
        lines.extend(("", f"[ {section} ]"))
        lines.extend(" ".join(row) for row in interactions[section])
    lines.append("")
    return "\n".join(lines)


def write_periodic_opls_cm5_system(
    system: GlycamMappedSystem,
    topology: PeriodicOplsTopology,
    output_directory: str | Path,
) -> tuple[Path, Path, Path]:
    """Write a generated periodic OPLS-CM5 GROMACS system."""
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=True)
    gro_path = directory / "structure.gro"
    itp_path = directory / "cellulose_opls_cm5.itp"
    top_path = directory / "system.top"
    write_gro(system, gro_path)
    itp_path.write_text(render_periodic_opls_cm5_itp(topology), encoding="utf-8")
    top_path.write_text(
        '#include "oplsaa.ff/forcefield.itp"\n'
        '#include "cellulose_opls_cm5.itp"\n\n'
        "[ system ]\nPeriodic cellulose OPLS-CM5\n\n"
        "[ molecules ]\nCNC  " + str(system.chain_count) + "\n",
        encoding="utf-8",
    )
    return gro_path, itp_path, top_path


def write_opls_cm5_system(
    system: GlycamMappedSystem,
    topology_template: str | Path,
    output_directory: str | Path,
) -> tuple[Path, Path, Path]:
    """Write GRO, ITP, and TOP files using an externally supplied reference."""
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=True)
    gro_path = directory / "structure.gro"
    itp_path = directory / "cellulose_opls_cm5.itp"
    top_path = directory / "system.top"
    write_gro(system, gro_path)
    shutil.copyfile(topology_template, itp_path)
    top_path.write_text(
        '#include "oplsaa.ff/forcefield.itp"\n'
        '#include "cellulose_opls_cm5.itp"\n\n'
        "[ system ]\nCellulose I-beta OPLS-CM5\n\n"
        "[ molecules ]\nCNC  " + str(system.chain_count) + "\n",
        encoding="utf-8",
    )
    return gro_path, itp_path, top_path
