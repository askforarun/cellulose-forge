"""Experimental adapter for the published pristine cellulose OPLS-CM5 topology."""

from __future__ import annotations

import shutil
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
            "the experimental OPLS-CM5 adapter requires the published pristine "
            f"DP20 topology with {_REFERENCE_ATOM_COUNT} atoms; found {len(atoms)}"
        )
        raise ValueError(msg)
    total_charge = sum(atom.charge for atom in atoms)
    if not np.isclose(total_charge, 0.0, atol=1.0e-6):
        msg = f"OPLS-CM5 reference topology is not neutral (charge {total_charge:.8f})"
        raise ValueError(msg)
    return tuple(atoms)


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
    """Map a finite DP20 crystal onto the exact published OPLS-CM5 atom order."""
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
        msg = "the published pristine OPLS-CM5 reference currently supports finite DP20 only"
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
