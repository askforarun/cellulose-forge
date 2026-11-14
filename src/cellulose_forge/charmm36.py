"""Native CHARMM36 topology construction and crystallographic mapping."""

from __future__ import annotations

import shutil
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cellulose_forge.cellulose_ii import cellulose_ii_preferred_acceptor
from cellulose_forge.clashes import validate_finite_periodic_image_contact
from cellulose_forge.crystallography import FloatArray, load_allomorph
from cellulose_forge.geometry import BulkCrystal, build_bulk_crystal
from cellulose_forge.glycam_mapping import GlycamMappedSystem
from cellulose_forge.gromacs import write_gro
from cellulose_forge.hydroxyl import HydroxylHydrogenPlacer
from cellulose_forge.restraints import position_restraints_section

_MINIMUM_ATOM_FIELDS = 8
_INTERACTION_ARITY = {"bonds": 2, "pairs": 2, "angles": 3, "dihedrals": 4}
_MINIMUM_PERIODIC_DP = 4
_MINIMUM_REFERENCE_DP = 6
_RING_ALIGNMENT_ATOMS = ("C1", "C3", "C5")
_MISSING_CELLULOSE_II_HYDROGENS = frozenset({"HO2", "HO3", "HO6"})
_HYDROXYL_OXYGEN_BY_HYDROGEN = {"HO2": "O2", "HO3": "O3", "HO6": "O6"}
_MAXIMUM_GLYCOSIDIC_BOND_ANGSTROM = 2.0
# CHARMM36 BGLCNA names for the N-acetyl group and the crystal labels used by
# the chitin definitions. Glucose (BGLC) has none of these names.
_CRYSTAL_NAME_BY_CHARMM_NAME = {
    "N": "N2",
    "HN": "HN2",
    "C": "C7",
    "O": "O7",
    "CT": "C8",
    "HT1": "H81",
    "HT2": "H82",
    "HT3": "H83",
}
_MINIMUM_ATOMTYPE_FIELDS = 3
_NATIVE_REFERENCE_DP = 8
_PATCH_DELETIONS = frozenset({"O1", "HO1", "HO4"})
_PATCH_ATOMS = {
    "C4": ("CC3161", 0.09),
    "O4": ("OC301", -0.36),
    "C1": ("CC3162", 0.29),
}
# Terminal geometry from the finite cellulose example distributed with doGlycans
# (CHARMM-GUI input POL_charmm-gui.pdb). Coordinates are used only as a local
# rigid-body geometry reference; force-field parameters come from the installed
# CHARMM36 release.
_REDUCING_END_REFERENCE = {
    "C1": np.array([7.209, 1.834, 9.925]),
    "C3": np.array([9.697, 1.806, 10.456]),
    "C5": np.array([7.961, 2.248, 12.274]),
    "O1": np.array([6.250, 2.279, 8.977]),
    "HO1": np.array([6.039, 1.562, 8.374]),
}
_NONREDUCING_END_REFERENCE = {
    "C1": np.array([66.699, 29.681, 33.127]),
    "C3": np.array([68.997, 30.551, 33.788]),
    "C5": np.array([66.899, 31.257, 35.058]),
    "HO4": np.array([70.018, 32.142, 35.442]),
}


@dataclass(frozen=True)
class Charmm36TemplateAtom:
    """One atom in an externally supplied CHARMM36 chain template."""

    residue_number: int
    residue_name: str
    atom_name: str
    charge: float
    atom_type: str
    mass: float


@dataclass(frozen=True)
class _IndexedInteraction:
    section: str
    atom_indices: tuple[int, ...]
    parameters: tuple[str, ...]


@dataclass(frozen=True, order=True)
class _InteractionTemplate:
    section: str
    atoms: tuple[tuple[int, str], ...]
    parameters: tuple[str, ...]


@dataclass(frozen=True)
class Charmm36ChainTemplate:
    """Validated metadata for one CHARMM36 cellulose-chain ITP."""

    molecule_name: str
    atoms: tuple[Charmm36TemplateAtom, ...]
    interactions: tuple[_IndexedInteraction, ...]

    @property
    def degree_polymerization(self) -> int:
        """Return the number of residues in the chain template."""
        return len({atom.residue_number for atom in self.atoms})

    @property
    def total_charge(self) -> float:
        """Return the template charge in elementary-charge units."""
        return sum(atom.charge for atom in self.atoms)


@dataclass(frozen=True)
class PeriodicCharmm36Topology:
    """An arbitrary-DP cyclic topology derived from a periodic reference ITP."""

    degree_polymerization: int
    molecule_name: str
    topology_atoms: tuple[Charmm36TemplateAtom, ...]
    interactions: tuple[_InteractionTemplate, ...]

    @property
    def atom_count(self) -> int:
        """Return atoms in one periodic chain."""
        return self.degree_polymerization * len(self.topology_atoms)

    @property
    def total_charge(self) -> float:
        """Return the generated chain charge."""
        return self.degree_polymerization * sum(atom.charge for atom in self.topology_atoms)


@dataclass(frozen=True)
class FiniteCharmm36Topology:
    """One explicitly terminated CHARMM36 cellulose chain."""

    degree_polymerization: int
    molecule_name: str
    topology_atoms: tuple[Charmm36TemplateAtom, ...]
    interactions: tuple[_IndexedInteraction, ...]

    @property
    def atom_count(self) -> int:
        """Return atoms in one finite chain."""
        return len(self.topology_atoms)

    @property
    def total_charge(self) -> float:
        """Return the generated chain charge."""
        return sum(atom.charge for atom in self.topology_atoms)


def read_charmm36_chain_template(  # noqa: PLR0912
    path: str | Path,
    *,
    require_uniform_residues: bool = True,
) -> Charmm36ChainTemplate:
    """Read and validate a GROMACS CHARMM36 cellulose-chain ITP."""
    molecule_name: str | None = None
    atoms: list[Charmm36TemplateAtom] = []
    interactions: list[_IndexedInteraction] = []
    section = ""
    sections: set[str] = set()
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        stripped = line.split(";", maxsplit=1)[0].strip()
        if stripped.startswith("#"):
            continue
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped[1:-1].strip()
            sections.add(section)
            continue
        if not stripped:
            continue
        fields = stripped.split()
        if section == "moleculetype" and molecule_name is None:
            molecule_name = fields[0]
        elif section == "atoms":
            if len(fields) < _MINIMUM_ATOM_FIELDS:
                msg = f"malformed CHARMM36 [ atoms ] row: {line}"
                raise ValueError(msg)
            atoms.append(
                Charmm36TemplateAtom(
                    residue_number=int(fields[2]),
                    residue_name=fields[3],
                    atom_name=fields[4],
                    charge=float(fields[6]),
                    atom_type=fields[1],
                    mass=float(fields[7]),
                )
            )
        elif section in _INTERACTION_ARITY:
            arity = _INTERACTION_ARITY[section]
            if len(fields) <= arity:
                msg = f"malformed CHARMM36 [{section}] row: {line}"
                raise ValueError(msg)
            interactions.append(
                _IndexedInteraction(
                    section=section,
                    atom_indices=tuple(int(value) - 1 for value in fields[:arity]),
                    parameters=tuple(fields[arity:]),
                )
            )
    if molecule_name is None or not atoms:
        msg = "CHARMM36 template must contain [ moleculetype ] and [ atoms ]"
        raise ValueError(msg)
    required_sections = {"bonds", "pairs", "angles", "dihedrals"}
    missing_sections = sorted(required_sections.difference(sections))
    if missing_sections:
        msg = "CHARMM36 template is missing sections: " + ", ".join(missing_sections)
        raise ValueError(msg)
    residues = tuple(dict.fromkeys(atom.residue_number for atom in atoms))
    if residues != tuple(range(1, len(residues) + 1)):
        msg = f"CHARMM36 template residues must be consecutively numbered: {residues}"
        raise ValueError(msg)
    names_by_residue = {
        residue: tuple(atom.atom_name for atom in atoms if atom.residue_number == residue)
        for residue in residues
    }
    if require_uniform_residues and len(set(names_by_residue.values())) != 1:
        msg = "periodic CHARMM36 template must have the same atoms in every residue"
        raise ValueError(msg)
    template = Charmm36ChainTemplate(
        molecule_name=molecule_name,
        atoms=tuple(atoms),
        interactions=tuple(interactions),
    )
    if not np.isclose(template.total_charge, 0.0, atol=1.0e-6):
        msg = f"CHARMM36 chain template is not neutral (charge {template.total_charge:.8f})"
        raise ValueError(msg)
    return template


def derive_periodic_charmm36_topology(
    path: str | Path,
    *,
    degree_polymerization: int,
) -> PeriodicCharmm36Topology:
    """Derive an arbitrary even-DP cyclic topology from a periodic reference."""
    if degree_polymerization < _MINIMUM_PERIODIC_DP or degree_polymerization % 2:
        msg = "periodic CHARMM36 chains require an even DP of at least 4"
        raise ValueError(msg)
    source = read_charmm36_chain_template(path)
    return _derive_periodic_topology(source, degree_polymerization)


def _derive_periodic_topology(
    source: Charmm36ChainTemplate,
    degree_polymerization: int,
) -> PeriodicCharmm36Topology:
    """Derive one cyclic repeat pattern from a validated periodic source."""
    if source.degree_polymerization < _MINIMUM_REFERENCE_DP:
        msg = "CHARMM36 reference must contain at least six periodic residues"
        raise ValueError(msg)
    reference_residue = source.degree_polymerization // 2
    topology_atoms = tuple(
        atom for atom in source.atoms if atom.residue_number == reference_residue
    )
    repeat_charge = sum(atom.charge for atom in topology_atoms)
    if not np.isclose(repeat_charge, 0.0, atol=1.0e-6):
        msg = f"CHARMM36 internal repeat is not neutral ({repeat_charge:.8f})"
        raise ValueError(msg)

    lower = 3
    upper = source.degree_polymerization - 2
    templates: set[_InteractionTemplate] = set()
    for interaction in source.interactions:
        keyed_atoms = tuple(
            (
                source.atoms[index].residue_number,
                source.atoms[index].atom_name,
            )
            for index in interaction.atom_indices
        )
        residues = tuple(residue for residue, _ in keyed_atoms)
        if min(residues) < lower or max(residues) > upper:
            continue
        anchor = min(residues)
        templates.add(
            _InteractionTemplate(
                section=interaction.section,
                atoms=tuple((residue - anchor, name) for residue, name in keyed_atoms),
                parameters=interaction.parameters,
            )
        )
    source_counts = {
        section: sum(item.section == section for item in source.interactions)
        for section in _INTERACTION_ARITY
    }
    expected_counts = {
        section: count // source.degree_polymerization for section, count in source_counts.items()
    }
    if any(count % source.degree_polymerization for count in source_counts.values()):
        msg = f"CHARMM36 reference interactions are not periodic: {source_counts}"
        raise ValueError(msg)
    actual_counts = {
        section: sum(item.section == section for item in templates)
        for section in _INTERACTION_ARITY
    }
    if actual_counts != expected_counts:
        msg = (
            "CHARMM36 repeating interaction pattern is incomplete: "
            f"expected {expected_counts}, found {actual_counts}"
        )
        raise ValueError(msg)
    return PeriodicCharmm36Topology(
        degree_polymerization=degree_polymerization,
        molecule_name=f"CARB_{degree_polymerization}",
        topology_atoms=topology_atoms,
        interactions=tuple(sorted(templates)),
    )


def _read_rtp_residue(
    force_field_directory: Path,
    residue: str,
) -> tuple[
    list[tuple[str, str, float]],
    list[tuple[str, str]],
    list[tuple[str, str, str, str]],
]:
    """Read one residue's atoms, bonds, and impropers from the installed RTP."""
    block = ""
    section = ""
    atoms: list[tuple[str, str, float]] = []
    bonds: list[tuple[str, str]] = []
    impropers: list[tuple[str, str, str, str]] = []
    for raw_line in (
        (force_field_directory / "carb.rtp")
        .read_text(encoding="utf-8", errors="replace")
        .splitlines()
    ):
        line = raw_line.split(";", maxsplit=1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            name = line[1:-1].strip()
            if name == residue:
                block = residue
                section = ""
            elif block == residue and name in {"atoms", "bonds", "impropers"}:
                section = name
            elif block == residue:
                break
            continue
        fields = line.split()
        if block != residue:
            continue
        if section == "atoms":
            atoms.append((fields[0], fields[1], float(fields[2])))
        elif section == "bonds":
            bonds.append((fields[0], fields[1]))
        elif section == "impropers":
            impropers.append((fields[0], fields[1], fields[2], fields[3]))
    if not atoms or not bonds:
        msg = f"could not read {residue} atoms and bonds from {force_field_directory / 'carb.rtp'}"
        raise ValueError(msg)
    return atoms, bonds, impropers


def _read_residue(
    force_field_directory: Path,
    residue: str,
) -> tuple[
    list[tuple[str, str, float]],
    list[tuple[str, str]],
    list[tuple[str, str, str, str]],
]:
    """Read a residue; BGLC (no impropers) keeps its dedicated reader."""
    if residue == "BGLC":
        atoms, bonds = _read_bglc_rtp(force_field_directory)
        return atoms, bonds, []
    return _read_rtp_residue(force_field_directory, residue)


def _read_bglc_rtp(
    force_field_directory: Path,
) -> tuple[list[tuple[str, str, float]], list[tuple[str, str]]]:
    """Read the BGLC atoms and bonds from the installed GROMACS RTP."""
    atoms, bonds, _ = _read_rtp_residue(force_field_directory, "BGLC")
    return atoms, bonds


def _read_atom_type_masses(force_field_directory: Path) -> dict[str, float]:
    """Read atom-type masses from the installed nonbonded parameter table."""
    masses: dict[str, float] = {}
    section = ""
    for raw_line in (
        (force_field_directory / "ffnonbonded.itp")
        .read_text(encoding="utf-8", errors="replace")
        .splitlines()
    ):
        line = raw_line.split(";", maxsplit=1)[0].strip()
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1].strip()
            continue
        fields = line.split()
        if section == "atomtypes" and len(fields) >= _MINIMUM_ATOMTYPE_FIELDS:
            try:
                masses.setdefault(fields[0], float(fields[2]))
            except ValueError:
                continue
    return masses


def _enumerate_graph_interactions(
    bonds: set[tuple[int, int]],
) -> dict[str, set[tuple[int, ...]]]:
    """Generate CHARMM angles, proper torsions, and graph-distance-three pairs."""
    adjacency: dict[int, set[int]] = defaultdict(set)
    for first, second in bonds:
        adjacency[first].add(second)
        adjacency[second].add(first)

    angles: set[tuple[int, ...]] = set()
    for center, neighbors in adjacency.items():
        ordered = sorted(neighbors)
        for position, first in enumerate(ordered):
            for second in ordered[position + 1 :]:
                angle_row = (first, center, second)
                angles.add(min(angle_row, angle_row[::-1]))

    dihedrals: set[tuple[int, ...]] = set()
    for second, second_neighbors in adjacency.items():
        for third in second_neighbors:
            for first in second_neighbors - {third}:
                for fourth in adjacency[third] - {second}:
                    if first != fourth:
                        dihedral_row = (first, second, third, fourth)
                        dihedrals.add(min(dihedral_row, dihedral_row[::-1]))

    pairs: set[tuple[int, ...]] = set()
    for first in adjacency:
        frontier = {first}
        visited = {first}
        for _ in range(3):
            frontier = {
                neighbor for current in frontier for neighbor in adjacency[current]
            }.difference(visited)
            visited.update(frontier)
        pairs.update(tuple(sorted((first, fourth))) for fourth in frontier)
    return {
        "bonds": {tuple(sorted(bond)) for bond in bonds},
        "pairs": pairs,
        "angles": angles,
        "dihedrals": dihedrals,
    }


def native_periodic_charmm36_topology(
    force_field_directory: str | Path,
    *,
    degree_polymerization: int,
    residue: str = "BGLC",
) -> PeriodicCharmm36Topology:
    """Construct a cyclic chain topology from an RTP residue and the 14bb patch.

    ``residue`` is ``BGLC`` for cellulose or ``BGLCNA`` for chitin. Residue
    impropers are added per residue as function-2 dihedrals.
    """
    if degree_polymerization < _MINIMUM_PERIODIC_DP or degree_polymerization % 2:
        msg = "periodic CHARMM36 chains require an even DP of at least 4"
        raise ValueError(msg)
    directory = Path(force_field_directory)
    rtp_atoms, rtp_bonds, rtp_impropers = _read_residue(directory, residue)
    masses = _read_atom_type_masses(directory)
    patched_atoms: list[Charmm36TemplateAtom] = []
    for name, atom_type, charge in rtp_atoms:
        if name in _PATCH_DELETIONS:
            continue
        patched_type, patched_charge = _PATCH_ATOMS.get(name, (atom_type, charge))
        try:
            mass = masses[patched_type]
        except KeyError as error:
            msg = f"CHARMM36 atom type {patched_type} has no mass in ffnonbonded.itp"
            raise ValueError(msg) from error
        patched_atoms.append(
            Charmm36TemplateAtom(1, residue, name, patched_charge, patched_type, mass)
        )
    names = {atom.atom_name for atom in patched_atoms}
    patched_bonds = {
        (first, second) for first, second in rtp_bonds if first in names and second in names
    }
    atoms_per_residue = len(patched_atoms)
    index_by_residue_name = {
        (source_residue, atom.atom_name): source_residue * atoms_per_residue + position
        for source_residue in range(_NATIVE_REFERENCE_DP)
        for position, atom in enumerate(patched_atoms)
    }
    indexed_bonds: set[tuple[int, int]] = set()
    for source_residue in range(_NATIVE_REFERENCE_DP):
        for first, second in patched_bonds:
            indexed_bonds.add(
                (
                    index_by_residue_name[(source_residue, first)],
                    index_by_residue_name[(source_residue, second)],
                )
            )
        following = (source_residue + 1) % _NATIVE_REFERENCE_DP
        indexed_bonds.add(
            (
                index_by_residue_name[(source_residue, "O4")],
                index_by_residue_name[(following, "C1")],
            )
        )
    graph = _enumerate_graph_interactions(indexed_bonds)
    parameters = {"bonds": "1", "pairs": "1", "angles": "5", "dihedrals": "9"}
    source_atoms = tuple(
        Charmm36TemplateAtom(
            residue + 1,
            atom.residue_name,
            atom.atom_name,
            atom.charge,
            atom.atom_type,
            atom.mass,
        )
        for residue in range(_NATIVE_REFERENCE_DP)
        for atom in patched_atoms
    )
    source_interactions = tuple(
        _IndexedInteraction(section, row, (parameters[section],))
        for section in ("bonds", "pairs", "angles", "dihedrals")
        for row in sorted(graph[section])
    ) + tuple(
        _IndexedInteraction(
            "dihedrals",
            tuple(index_by_residue_name[(position, name)] for name in improper),
            ("2",),
        )
        for position in range(_NATIVE_REFERENCE_DP)
        for improper in rtp_impropers
    )
    source = Charmm36ChainTemplate("CARB_8", source_atoms, source_interactions)
    return _derive_periodic_topology(source, degree_polymerization)


def native_finite_charmm36_topology(  # noqa: PLR0912
    force_field_directory: str | Path,
    *,
    degree_polymerization: int,
    residue: str = "BGLC",
) -> FiniteCharmm36Topology:
    """Construct a neutral, explicitly terminated chain of one RTP residue.

    ``residue`` is ``BGLC`` for cellulose or ``BGLCNA`` for chitin. Residue
    impropers are added per residue as function-2 dihedrals.
    """
    if degree_polymerization < 1:
        msg = "finite CHARMM36 chains require a positive DP"
        raise ValueError(msg)
    directory = Path(force_field_directory)
    rtp_atoms, rtp_bonds, rtp_impropers = _read_residue(directory, residue)
    masses = _read_atom_type_masses(directory)
    topology_atoms: list[Charmm36TemplateAtom] = []
    index_by_residue_name: dict[tuple[int, str], int] = {}
    for position in range(1, degree_polymerization + 1):
        linked_from_previous = position > 1
        linked_to_next = position < degree_polymerization
        for name, base_atom_type, base_charge in rtp_atoms:
            if linked_from_previous and name in {"O1", "HO1"}:
                continue
            if linked_to_next and name == "HO4":
                continue
            atom_type, charge = base_atom_type, base_charge
            if linked_from_previous and name == "C1":
                atom_type, charge = _PATCH_ATOMS[name]
            if linked_to_next and name in {"C4", "O4"}:
                atom_type, charge = _PATCH_ATOMS[name]
            try:
                mass = masses[atom_type]
            except KeyError as error:
                msg = f"CHARMM36 atom type {atom_type} has no mass in ffnonbonded.itp"
                raise ValueError(msg) from error
            index_by_residue_name[(position, name)] = len(topology_atoms)
            topology_atoms.append(
                Charmm36TemplateAtom(position, residue, name, charge, atom_type, mass)
            )

    bonds: set[tuple[int, int]] = set()
    for position in range(1, degree_polymerization + 1):
        for first, second in rtp_bonds:
            first_key = (position, first)
            second_key = (position, second)
            if first_key in index_by_residue_name and second_key in index_by_residue_name:
                bonds.add((index_by_residue_name[first_key], index_by_residue_name[second_key]))
        if position < degree_polymerization:
            bonds.add(
                (
                    index_by_residue_name[(position, "O4")],
                    index_by_residue_name[(position + 1, "C1")],
                )
            )
    graph = _enumerate_graph_interactions(bonds)
    functions = {"bonds": "1", "pairs": "1", "angles": "5", "dihedrals": "9"}
    interactions = tuple(
        _IndexedInteraction(section, row, (functions[section],))
        for section in ("bonds", "pairs", "angles", "dihedrals")
        for row in sorted(graph[section])
    ) + tuple(
        _IndexedInteraction(
            "dihedrals",
            tuple(index_by_residue_name[(position, name)] for name in improper),
            ("2",),
        )
        for position in range(1, degree_polymerization + 1)
        for improper in rtp_impropers
    )
    topology = FiniteCharmm36Topology(
        degree_polymerization,
        f"CARB_{degree_polymerization}_FINITE",
        tuple(topology_atoms),
        interactions,
    )
    if not np.isclose(topology.total_charge, 0.0, atol=1.0e-6):
        msg = f"finite CHARMM36 chain is not neutral ({topology.total_charge:.8f})"
        raise ValueError(msg)
    return topology


def _crystal_name(charmm_name: str) -> str:
    """Return the crystal-definition label for a CHARMM36 atom name."""
    return _CRYSTAL_NAME_BY_CHARMM_NAME.get(charmm_name, charmm_name)


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
) -> tuple[FloatArray, FloatArray, FloatArray]:
    reference_points = np.vstack([reference[name] for name in _RING_ALIGNMENT_ATOMS])
    target_points = np.vstack([target[name] for name in _RING_ALIGNMENT_ATOMS])
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


def _hydrogen_reference() -> dict[str, FloatArray]:
    reference = build_bulk_crystal(load_allomorph("i_beta"), (1, 1, 1))
    return _coordinate_lookup(reference, chain_number=1, residue_number=1)


def _chain_residue_order(
    crystal: BulkCrystal,
    *,
    chain_number: int,
    degree_polymerization: int,
) -> tuple[int, ...]:
    """Order crystal residues by the periodic O4-to-C1 glycosidic path."""
    residues = tuple(range(1, degree_polymerization + 1))
    lookups = {
        residue: _coordinate_lookup(crystal, chain_number=chain_number, residue_number=residue)
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


def _finite_chain_residue_order(
    crystal: BulkCrystal,
    *,
    chain_number: int,
    degree_polymerization: int,
) -> tuple[int, ...]:
    """Order residues along the open O4-to-C1 glycosidic path."""
    residues = tuple(range(1, degree_polymerization + 1))
    lookups = {
        residue: _coordinate_lookup(crystal, chain_number=chain_number, residue_number=residue)
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


def map_periodic_charmm36(
    crystal: BulkCrystal,
    topology: PeriodicCharmm36Topology,
) -> GlycamMappedSystem:
    """Order a periodic crystal exactly like a CHARMM36 chain template."""
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
    if not degree_by_chain or set(degree_by_chain.values()) != {topology.degree_polymerization}:
        msg = (
            "CHARMM36 topology DP does not match the generated crystal chains: "
            f"topology DP {topology.degree_polymerization}, "
            f"crystal DPs {sorted(set(degree_by_chain.values()))}"
        )
        raise ValueError(msg)

    output_names: list[str] = []
    output_residue_names: list[str] = []
    output_residue_numbers: list[int] = []
    output_chain_numbers: list[int] = []
    output_coordinates: list[FloatArray] = []
    generated_names = _MISSING_CELLULOSE_II_HYDROGENS if crystal.allomorph == "II" else frozenset()
    hydrogen_reference = _hydrogen_reference() if generated_names else {}
    placer = (
        HydroxylHydrogenPlacer(
            crystal.atom_labels, crystal.cartesian_coordinates, crystal.basis_vectors
        )
        if generated_names
        else None
    )
    for chain_number in chain_numbers:
        crystal_residue_order = _chain_residue_order(
            crystal,
            chain_number=chain_number,
            degree_polymerization=topology.degree_polymerization,
        )
        for residue_number in range(1, topology.degree_polymerization + 1):
            crystal_residue_number = crystal_residue_order[residue_number - 1]
            lookup = _coordinate_lookup(
                crystal,
                chain_number=chain_number,
                residue_number=crystal_residue_number,
            )
            expected_names = {_crystal_name(atom.atom_name) for atom in topology.topology_atoms}
            missing = expected_names.difference(lookup).difference(generated_names)
            if missing:
                msg = (
                    f"chain {chain_number} residue {residue_number} lacks CHARMM36 "
                    f"atoms: {', '.join(sorted(missing))}"
                )
                raise ValueError(msg)
            hydrogen_transform = (
                _rigid_transform(hydrogen_reference, lookup) if generated_names else None
            )
            global_residue = (chain_number - 1) * topology.degree_polymerization + residue_number
            for atom in topology.topology_atoms:
                output_names.append(atom.atom_name)
                output_residue_names.append(atom.residue_name)
                output_residue_numbers.append(global_residue)
                output_chain_numbers.append(chain_number)
                if atom.atom_name in generated_names:
                    assert hydrogen_transform is not None
                    oxygen_name = _HYDROXYL_OXYGEN_BY_HYDROGEN[atom.atom_name]
                    aligned_hydrogen = _transform(
                        hydrogen_reference[atom.atom_name], hydrogen_transform
                    )
                    assert placer is not None
                    placer.request(
                        len(output_coordinates),
                        lookup["C" + oxygen_name[1:]],
                        lookup[oxygen_name],
                        aligned_hydrogen,
                        donor=oxygen_name,
                        preferred_acceptor="O5" if oxygen_name == "O3" else None,
                        preferred_acceptor_coordinate=cellulose_ii_preferred_acceptor(
                            crystal,
                            chain_number=chain_number,
                            residue_number=crystal_residue_number,
                            donor=oxygen_name,
                        ),
                    )
                    output_coordinates.append(np.zeros(3))
                else:
                    output_coordinates.append(lookup[_crystal_name(atom.atom_name)])

    if placer is not None:
        for index, coordinate in placer.resolve().items():
            output_coordinates[index] = coordinate
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


def map_finite_charmm36(
    crystal: BulkCrystal,
    topology: FiniteCharmm36Topology,
) -> GlycamMappedSystem:
    """Add terminal atoms and order a finite crystal for native CHARMM36."""
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
    if not degree_by_chain or set(degree_by_chain.values()) != {topology.degree_polymerization}:
        msg = (
            "CHARMM36 topology DP does not match the generated crystal chains: "
            f"topology DP {topology.degree_polymerization}, "
            f"crystal DPs {sorted(set(degree_by_chain.values()))}"
        )
        raise ValueError(msg)

    atoms_by_residue: dict[int, tuple[Charmm36TemplateAtom, ...]] = {
        residue: tuple(atom for atom in topology.topology_atoms if atom.residue_number == residue)
        for residue in range(1, topology.degree_polymerization + 1)
    }
    output_names: list[str] = []
    output_residue_names: list[str] = []
    output_residue_numbers: list[int] = []
    output_chain_numbers: list[int] = []
    output_coordinates: list[FloatArray] = []
    generated_hydrogens = (
        _MISSING_CELLULOSE_II_HYDROGENS if crystal.allomorph == "II" else frozenset()
    )
    hydrogen_reference = _hydrogen_reference() if generated_hydrogens else {}
    placer = (
        HydroxylHydrogenPlacer(
            crystal.atom_labels, crystal.cartesian_coordinates, crystal.basis_vectors
        )
        if generated_hydrogens
        else None
    )
    for chain_number in chain_numbers:
        crystal_residue_order = _finite_chain_residue_order(
            crystal,
            chain_number=chain_number,
            degree_polymerization=topology.degree_polymerization,
        )
        for residue_number in range(1, topology.degree_polymerization + 1):
            crystal_residue_number = crystal_residue_order[residue_number - 1]
            lookup = _coordinate_lookup(
                crystal,
                chain_number=chain_number,
                residue_number=crystal_residue_number,
            )
            reducing_transform = (
                _rigid_transform(_REDUCING_END_REFERENCE, lookup) if residue_number == 1 else None
            )
            nonreducing_transform = (
                _rigid_transform(_NONREDUCING_END_REFERENCE, lookup)
                if residue_number == topology.degree_polymerization
                else None
            )
            hydrogen_transform = (
                _rigid_transform(hydrogen_reference, lookup) if generated_hydrogens else None
            )
            global_residue = (chain_number - 1) * topology.degree_polymerization + residue_number
            for atom in atoms_by_residue[residue_number]:
                output_names.append(atom.atom_name)
                output_residue_names.append(atom.residue_name)
                output_residue_numbers.append(global_residue)
                output_chain_numbers.append(chain_number)
                if atom.atom_name in {"O1", "HO1"}:
                    assert reducing_transform is not None
                    coordinate = _transform(
                        _REDUCING_END_REFERENCE[atom.atom_name], reducing_transform
                    )
                elif atom.atom_name == "HO4":
                    assert nonreducing_transform is not None
                    coordinate = _transform(
                        _NONREDUCING_END_REFERENCE[atom.atom_name], nonreducing_transform
                    )
                elif atom.atom_name in generated_hydrogens:
                    assert hydrogen_transform is not None
                    oxygen_name = _HYDROXYL_OXYGEN_BY_HYDROGEN[atom.atom_name]
                    aligned_hydrogen = _transform(
                        hydrogen_reference[atom.atom_name], hydrogen_transform
                    )
                    assert placer is not None
                    placer.request(
                        len(output_coordinates),
                        lookup["C" + oxygen_name[1:]],
                        lookup[oxygen_name],
                        aligned_hydrogen,
                        donor=oxygen_name,
                        preferred_acceptor="O5" if oxygen_name == "O3" else None,
                        preferred_acceptor_coordinate=cellulose_ii_preferred_acceptor(
                            crystal,
                            chain_number=chain_number,
                            residue_number=crystal_residue_number,
                            donor=oxygen_name,
                        ),
                    )
                    coordinate = np.zeros(3)
                else:
                    coordinate = lookup[_crystal_name(atom.atom_name)]
                output_coordinates.append(coordinate)
    if placer is not None:
        for index, coordinate in placer.resolve().items():
            output_coordinates[index] = coordinate
    coordinates = np.vstack(output_coordinates)
    coordinates.setflags(write=False)
    system = GlycamMappedSystem(
        atom_names=tuple(output_names),
        residue_names=tuple(output_residue_names),
        residue_numbers=tuple(output_residue_numbers),
        chain_numbers=tuple(output_chain_numbers),
        cartesian_coordinates=coordinates,
        basis_vectors=crystal.basis_vectors,
    )
    validate_finite_periodic_image_contact(system)
    return system


def _instantiate_periodic_interactions(
    topology: PeriodicCharmm36Topology,
) -> dict[str, tuple[tuple[str, ...], ...]]:
    atoms_per_residue = len(topology.topology_atoms)
    position_by_name = {
        atom.atom_name: position for position, atom in enumerate(topology.topology_atoms)
    }
    instantiated: dict[str, set[tuple[str, ...]]] = defaultdict(set)
    for anchor in range(topology.degree_polymerization):
        for template in topology.interactions:
            indices = tuple(
                str(
                    ((anchor + offset) % topology.degree_polymerization) * atoms_per_residue
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


def periodic_charmm36_interaction_counts(
    topology: PeriodicCharmm36Topology,
) -> dict[str, int]:
    """Return interaction counts for one generated cyclic chain."""
    return {
        section: len(rows) for section, rows in _instantiate_periodic_interactions(topology).items()
    }


def render_periodic_charmm36_itp(topology: PeriodicCharmm36Topology) -> str:
    """Render a generated cyclic CHARMM36 cellulose chain."""
    lines = [
        "[ moleculetype ]",
        "; name nrexcl",
        f"{topology.molecule_name} 3",
        "",
        "[ atoms ]",
        "; nr type resnr residue atom cgnr charge mass",
    ]
    atom_number = 0
    for residue_number in range(1, topology.degree_polymerization + 1):
        for atom in topology.topology_atoms:
            atom_number += 1
            lines.append(
                f"{atom_number:6d} {atom.atom_type:<10s} {residue_number:6d} "
                f"{atom.residue_name:<5s} {atom.atom_name:<5s} {atom_number:6d} "
                f"{atom.charge: .8f} {atom.mass:.6f}"
            )
    interactions = _instantiate_periodic_interactions(topology)
    for section in ("bonds", "pairs", "angles", "dihedrals"):
        lines.extend(("", f"[ {section} ]"))
        lines.extend(" ".join(row) for row in interactions[section])
    lines.append("")
    return "\n".join(lines)


def finite_charmm36_interaction_counts(
    topology: FiniteCharmm36Topology,
) -> dict[str, int]:
    """Return interaction counts for one generated finite chain."""
    return {
        section: sum(item.section == section for item in topology.interactions)
        for section in _INTERACTION_ARITY
    }


def render_finite_charmm36_itp(topology: FiniteCharmm36Topology) -> str:
    """Render a generated, explicitly terminated CHARMM36 chain."""
    lines = [
        "[ moleculetype ]",
        "; name nrexcl",
        f"{topology.molecule_name} 3",
        "",
        "[ atoms ]",
        "; nr type resnr residue atom cgnr charge mass",
    ]
    for atom_number, atom in enumerate(topology.topology_atoms, start=1):
        lines.append(
            f"{atom_number:6d} {atom.atom_type:<10s} {atom.residue_number:6d} "
            f"{atom.residue_name:<5s} {atom.atom_name:<5s} {atom_number:6d} "
            f"{atom.charge: .8f} {atom.mass:.6f}"
        )
    for section in ("bonds", "pairs", "angles", "dihedrals"):
        lines.extend(("", f"[ {section} ]"))
        for interaction in topology.interactions:
            if interaction.section == section:
                indices = (str(index + 1) for index in interaction.atom_indices)
                lines.append(" ".join((*indices, *interaction.parameters)))
    lines.append("")
    return "\n".join(lines)


def write_charmm36_system(  # noqa: PLR0913
    system: GlycamMappedSystem,
    topology: PeriodicCharmm36Topology | FiniteCharmm36Topology,
    force_field_directory: str | Path,
    output_directory: str | Path,
    *,
    restrained_chain_count: int = 0,
    restraint_force: tuple[float, float, float] | None = None,
    polymer: str = "cellulose",
) -> tuple[Path, Path, Path, Path]:
    """Write a CHARMM36 GROMACS system and local FF copy.

    Cellulose keeps the historical ``cellulose_charmm36.itp`` name; chitin
    writes ``chitin_charmm36.itp``.
    """
    source_force_field = Path(force_field_directory)
    if not (source_force_field / "forcefield.itp").is_file():
        msg = "CHARMM36 force-field directory must contain forcefield.itp"
        raise ValueError(msg)
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=True)
    gro_path = directory / "structure.gro"
    itp_path = directory / f"{polymer}_charmm36.itp"
    top_path = directory / "system.top"
    force_field_path = directory / "charmm36.ff"
    write_gro(system, gro_path)
    free_chain_count = system.chain_count - restrained_chain_count
    if restrained_chain_count < 0 or free_chain_count < 0:
        msg = "restrained chain count must be between zero and the total chain count"
        raise ValueError(msg)
    base_itp = (
        render_periodic_charmm36_itp(topology)
        if isinstance(topology, PeriodicCharmm36Topology)
        else render_finite_charmm36_itp(topology)
    )
    itp_text = base_itp
    molecule_rows: list[str] = []
    if free_chain_count:
        molecule_rows.append(f"{topology.molecule_name}  {free_chain_count}")
    if restrained_chain_count:
        if restraint_force is None:
            msg = "restraint force constants are required for restrained chains"
            raise ValueError(msg)
        restrained_name = f"{topology.molecule_name}_RESTRAINED"
        restrained_itp = base_itp.replace(f"{topology.molecule_name} 3", f"{restrained_name} 3", 1)
        restraint_atom_names = (
            tuple(
                atom.atom_name
                for _ in range(topology.degree_polymerization)
                for atom in topology.topology_atoms
            )
            if isinstance(topology, PeriodicCharmm36Topology)
            else tuple(atom.atom_name for atom in topology.topology_atoms)
        )
        restrained_itp += position_restraints_section(
            restraint_atom_names, restraint_force, define=f"POSRES_{polymer.upper()}"
        )
        itp_text += "\n" + restrained_itp
        molecule_rows.append(f"{restrained_name}  {restrained_chain_count}")
    itp_path.write_text(itp_text, encoding="utf-8")
    shutil.copytree(source_force_field, force_field_path, dirs_exist_ok=True)
    top_path.write_text(
        '#include "charmm36.ff/forcefield.itp"\n'
        f'#include "{itp_path.name}"\n\n'
        + (
            "[ system ]\nPeriodic cellulose CHARMM36\n\n"
            if polymer == "cellulose"
            else f"[ system ]\n{polymer} CHARMM36\n\n"
        )
        + "[ molecules ]\n"
        + "\n".join(molecule_rows)
        + "\n",
        encoding="utf-8",
    )
    return gro_path, itp_path, top_path, force_field_path
