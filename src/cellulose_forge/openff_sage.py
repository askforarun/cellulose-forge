"""Experimental OpenFF Sage construction for crystalline polysaccharides.

OpenFF is imported lazily so the core CelluloseForge installation retains its
NumPy-only coordinate path.  Molecular graphs are constructed directly from
CelluloseForge atom names and bonds; no AmberTools or GLYCAM template is used.
"""

from __future__ import annotations

import functools
import hashlib
import importlib
import importlib.metadata
import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from cellulose_forge.crystallography import FloatArray
from cellulose_forge.geometry import BulkCrystal
from cellulose_forge.glycam_mapping import GlycamMappedSystem
from cellulose_forge.openmm import (
    OPENMM_RUNNER_FILE,
    OPENMM_STATE_FILE,
    chain_axis,
    render_openmm_runner,
    validation_platform_for_box,
)
from cellulose_forge.polymer_mapping import (
    chain_template,
    map_crystal_to_chain_template,
    minimum_image,
    reduced_periodic_box,
)
from cellulose_forge.polymer_profiles import polymer_profile
from cellulose_forge.run_inputs import dispersion_correction as _dispersion_policy
from cellulose_forge.run_inputs import write_gromacs_run_inputs

SAGE_FORCE_FIELD = "openff-2.3.0.offxml"
NAGL_CHARGE_MODEL = "openff-gnn-am1bcc-1.0.0.pt"
VALIDATION_SCOPE = (
    "model construction validated; force-field accuracy for the selected polymer not assessed"
)

_INTERNAL_ATOMS = (
    ("C1", 6),
    ("H1", 1),
    ("O5", 8),
    ("C5", 6),
    ("H5", 1),
    ("C6", 6),
    ("H62", 1),
    ("H61", 1),
    ("O6", 8),
    ("HO6", 1),
    ("C4", 6),
    ("H4", 1),
    ("C3", 6),
    ("H3", 1),
    ("O3", 8),
    ("HO3", 1),
    ("C2", 6),
    ("H2", 1),
    ("O2", 8),
    ("HO2", 1),
    ("O4", 8),
)
_INTERNAL_BONDS = (
    ("C1", "H1"),
    ("C1", "O5"),
    ("C1", "C2"),
    ("O5", "C5"),
    ("C5", "H5"),
    ("C5", "C6"),
    ("C5", "C4"),
    ("C6", "H62"),
    ("C6", "H61"),
    ("C6", "O6"),
    ("O6", "HO6"),
    ("C4", "H4"),
    ("C4", "C3"),
    ("C4", "O4"),
    ("C3", "H3"),
    ("C3", "O3"),
    ("C3", "C2"),
    ("O3", "HO3"),
    ("C2", "H2"),
    ("C2", "O2"),
    ("O2", "HO2"),
)
_CHITIN_INTERNAL_ATOMS = (
    *(atom for atom in _INTERNAL_ATOMS if atom[0] not in {"O2", "HO2"}),
    ("N2", 7),
    ("HN2", 1),
    ("C7", 6),
    ("O7", 8),
    ("C8", 6),
    ("H81", 1),
    ("H82", 1),
    ("H83", 1),
)
_CHITIN_INTERNAL_BONDS = (
    *(bond for bond in _INTERNAL_BONDS if not ({"O2", "HO2"} & set(bond))),
    ("C2", "N2"),
    ("N2", "HN2"),
    ("N2", "C7"),
    ("C7", "O7"),
    ("C7", "C8"),
    ("C8", "H81"),
    ("C8", "H82"),
    ("C8", "H83"),
)
_STEREOCHEMISTRY = {"C1": "S", "C2": "R", "C3": "R", "C4": "S", "C5": "R"}
_REFERENCE_PARITY = {"C1": 1, "C2": 1, "C3": -1, "C4": 1, "C5": 1}
_STEREOCENTRES = frozenset(_STEREOCHEMISTRY)
_MINIMUM_DP = 2
_CHARGE_TOLERANCE_E = 1.0e-6
_BOX_DEGENERACY_TOLERANCE = 1.0e-12


@dataclass(frozen=True)
class SagePolymerProfile:
    """Residue graph and validation metadata for one Sage polymer."""

    name: str
    internal_atoms: tuple[tuple[str, int], ...]
    internal_bonds: tuple[tuple[str, str], ...]
    bond_orders: Mapping[frozenset[str], int]
    stereochemistry: dict[str, str]
    terminal_stereochemistry: Mapping[tuple[str, str], str]
    reference_parity: dict[str, int]

    def __post_init__(self) -> None:
        """Reject incomplete chemistry profiles at definition time."""
        try:
            polymer_profile(self.name)
        except ValueError as error:
            raise ValueError(f"Sage profile names unknown polymer {self.name!r}") from error
        expected_bonds = {frozenset(bond) for bond in self.internal_bonds}
        recorded_bonds = set(self.bond_orders)
        missing = expected_bonds - recorded_bonds
        extra = recorded_bonds - expected_bonds
        if missing or extra:
            raise ValueError(
                "Sage profile bond orders must cover every internal bond exactly; "
                f"missing={sorted(map(sorted, missing))}, extra={sorted(map(sorted, extra))}"
            )
        required_terminal = {("reducing", "C1"), ("nonreducing", "C3")}
        if set(self.terminal_stereochemistry) != required_terminal:
            raise ValueError(
                "Sage profile terminal stereochemistry must define reducing C1 and nonreducing C3"
            )

    @property
    def residue_name(self) -> str:
        """Return the canonical residue name from the shared polymer profile."""
        return polymer_profile(self.name).residue_name


_POLYMER_PROFILES = {
    "cellulose": SagePolymerProfile(
        name="cellulose",
        internal_atoms=_INTERNAL_ATOMS,
        internal_bonds=_INTERNAL_BONDS,
        bond_orders={frozenset(bond): 1 for bond in _INTERNAL_BONDS},
        stereochemistry=_STEREOCHEMISTRY,
        terminal_stereochemistry={("reducing", "C1"): "R", ("nonreducing", "C3"): "S"},
        reference_parity=_REFERENCE_PARITY,
    ),
    "chitin": SagePolymerProfile(
        name="chitin",
        internal_atoms=_CHITIN_INTERNAL_ATOMS,
        internal_bonds=_CHITIN_INTERNAL_BONDS,
        bond_orders={
            frozenset(bond): 2 if frozenset(bond) == frozenset({"C7", "O7"}) else 1
            for bond in _CHITIN_INTERNAL_BONDS
        },
        stereochemistry=_STEREOCHEMISTRY,
        terminal_stereochemistry={("reducing", "C1"): "R", ("nonreducing", "C3"): "R"},
        reference_parity=_REFERENCE_PARITY,
    ),
}


@dataclass(frozen=True)
class SageBuild:
    """A charged OpenFF molecule plus mapped multi-chain coordinates."""

    molecule: Any
    mapped: GlycamMappedSystem
    degree_polymerization: int
    chain_boundary: str
    charge_provenance: dict[str, object]
    stereochemistry: dict[str, object]


def dispersion_correction_for_shape(shape: str) -> bool:
    """Return the shared dispersion policy for an OpenFF Sage build."""
    return _dispersion_policy("openff", shape)


def periodic_box_heights(box: FloatArray) -> tuple[float, float, float]:
    """Return perpendicular repeat distances for the three lattice faces."""
    vectors = np.asarray(box, dtype=np.float64)
    volume = abs(float(np.linalg.det(vectors)))
    heights = tuple(
        volume / float(np.linalg.norm(np.cross(vectors[(axis + 1) % 3], vectors[(axis + 2) % 3])))
        for axis in range(3)
    )
    if not all(np.isfinite(heights)) or min(heights) <= 0:
        raise ValueError("OpenFF output requires a non-degenerate periodic box")
    return (float(heights[0]), float(heights[1]), float(heights[2]))


def shortest_lattice_vector(box: FloatArray) -> float:
    """Return the shortest non-zero translation of a three-dimensional lattice."""
    vectors = np.asarray(box, dtype=np.float64)
    if (
        not np.isfinite(vectors).all()
        or abs(float(np.linalg.det(vectors))) <= _BOX_DEGENERACY_TOLERANCE
    ):
        raise ValueError("OpenFF output requires a non-degenerate periodic box")
    translations = (
        np.asarray((i, j, k), dtype=np.float64) @ vectors
        for i in range(-2, 3)
        for j in range(-2, 3)
        for k in range(-2, 3)
        if (i, j, k) != (0, 0, 0)
    )
    return min(float(np.linalg.norm(translation)) for translation in translations)


def require_openff() -> tuple[Any, Any, Any]:
    """Import optional OpenFF components or raise an actionable error."""
    try:
        toolkit = importlib.import_module("openff.toolkit")
    except ImportError as error:
        raise RuntimeError(
            "OpenFF Sage support requires the optional OpenFF environment; "
            "install it with 'conda env create -f environment-openff.yml'"
        ) from error
    return toolkit.ForceField, toolkit.Molecule, toolkit.Topology


def _polymer_profile(polymer: str) -> SagePolymerProfile:
    try:
        return _POLYMER_PROFILES[polymer]
    except KeyError as error:
        raise ValueError(f"OpenFF Sage does not support polymer {polymer!r}") from error


def build_sage_molecule(  # noqa: PLR0912
    degree_polymerization: int,
    chain_boundary: str,
    polymer: str = "cellulose",
) -> Any:
    """Construct a finite or cyclic beta-(1->4) polysaccharide graph directly."""
    profile = _polymer_profile(polymer)
    if degree_polymerization < _MINIMUM_DP:
        raise ValueError(f"OpenFF {polymer} chains require DP >= 2")
    if chain_boundary not in {"finite", "periodic"}:
        raise ValueError("OpenFF chain boundary must be finite or periodic")
    _, molecule_type, _ = require_openff()
    molecule = molecule_type()
    molecule.name = f"{polymer}_dp{degree_polymerization}_{chain_boundary}"
    indices: dict[tuple[int, str], int] = {}
    for residue in range(1, degree_polymerization + 1):
        for name, atomic_number in profile.internal_atoms:
            stereo = profile.stereochemistry.get(name)
            if chain_boundary == "finite" and residue == 1 and name == "C1":
                stereo = profile.terminal_stereochemistry["reducing", "C1"]
            if chain_boundary == "finite" and residue == degree_polymerization and name == "C3":
                stereo = profile.terminal_stereochemistry["nonreducing", "C3"]
            indices[residue, name] = molecule.add_atom(
                atomic_number,
                0,
                False,
                stereo,
                name,
                {
                    "residue_name": profile.residue_name,
                    "residue_number": residue,
                    "chain_id": "A",
                },
            )
        if chain_boundary == "finite" and residue == 1:
            for name, atomic_number in (("O1", 8), ("HO1", 1)):
                indices[residue, name] = molecule.add_atom(
                    atomic_number,
                    0,
                    False,
                    None,
                    name,
                    {
                        "residue_name": profile.residue_name,
                        "residue_number": residue,
                        "chain_id": "A",
                    },
                )
        if chain_boundary == "finite" and residue == degree_polymerization:
            indices[residue, "HO4"] = molecule.add_atom(
                1,
                0,
                False,
                None,
                "HO4",
                {
                    "residue_name": profile.residue_name,
                    "residue_number": residue,
                    "chain_id": "A",
                },
            )
        for left, right in profile.internal_bonds:
            try:
                order = profile.bond_orders[frozenset({left, right})]
            except KeyError as error:
                raise ValueError(
                    f"Sage {profile.name} profile lacks a bond order for {left}-{right}"
                ) from error
            molecule.add_bond(indices[residue, left], indices[residue, right], order, False)
    for residue in range(1, degree_polymerization):
        molecule.add_bond(indices[residue, "O4"], indices[residue + 1, "C1"], 1, False)
    if chain_boundary == "periodic":
        molecule.add_bond(indices[degree_polymerization, "O4"], indices[1, "C1"], 1, False)
    else:
        molecule.add_bond(indices[1, "C1"], indices[1, "O1"], 1, False)
        molecule.add_bond(indices[1, "O1"], indices[1, "HO1"], 1, False)
        molecule.add_bond(
            indices[degree_polymerization, "O4"],
            indices[degree_polymerization, "HO4"],
            1,
            False,
        )
    return molecule


def _centre_parity(positions: FloatArray, atom: Any, box: FloatArray) -> int:
    neighbours = sorted(atom.bonded_atoms, key=lambda item: (item.symbol, item.name))[:3]
    centre = positions[atom.molecule_atom_index]
    vectors = [
        minimum_image(positions[item.molecule_atom_index] - centre, box) for item in neighbours
    ]
    value = float(np.dot(vectors[0], np.cross(vectors[1], vectors[2])))
    if not np.isfinite(value) or value == 0.0:
        raise ValueError("undefined cellulose stereocentre parity")
    return int(np.sign(value))


def validate_sage_stereochemistry(
    molecule: Any,
    mapped: GlycamMappedSystem,
    polymer: str = "cellulose",
) -> dict[str, object]:
    """Check mapped C1--C5 centres against the polymer's coordinate reference."""
    profile = _polymer_profile(polymer)
    failures: list[tuple[int, int, str]] = []
    for chain_index, chain in enumerate(dict.fromkeys(mapped.chain_numbers)):
        start = chain_index * molecule.n_atoms
        positions = mapped.cartesian_coordinates[start : start + molecule.n_atoms]
        for atom in molecule.atoms:
            if (
                atom.name in profile.reference_parity
                and _centre_parity(positions, atom, mapped.basis_vectors)
                != profile.reference_parity[atom.name]
            ):
                failures.append(  # noqa: PERF401
                    (chain, int(atom.metadata["residue_number"]), atom.name)
                )
    result: dict[str, object] = {
        "reference": polymer_profile(polymer).validation_reference,
        "centres_checked": mapped.chain_count
        * len({int(atom.metadata["residue_number"]) for atom in molecule.atoms})
        * len(profile.reference_parity),
        "failures": failures,
        "passed": not failures,
    }
    if failures:
        raise ValueError(f"OpenFF {polymer} stereochemistry check failed: {failures[:5]}")
    return result


def validate_declared_stereochemistry(
    molecule: Any,
    mapped: GlycamMappedSystem,
    chain_boundary: str,
    polymer: str = "cellulose",
) -> dict[str, object]:
    """Confirm every declared R/S label against the mapped coordinates.

    RDKit's global stereo perception treats a covalently periodic chain as one
    macrocycle and its cost grows explosively with chain length (seconds per
    chain at DP20), so it is not run on the built chains. Instead the CIP
    priority order of each centre's four neighbours, which depends only on the
    chemistry within a few bonds and is therefore the same for every residue
    of a kind (reducing, internal, non-reducing) at any chain length, is
    computed once by RDKit on a short chain of the same polymer and boundary.
    Each centre's handedness is then read from the minimum-image vectors to its
    bonded neighbours, so wrapped periodic coordinates need no unwrapping. The
    cost is linear in the number of centres.
    """
    if chain_boundary not in _REFERENCE_DEGREE:
        raise ValueError(
            "OpenFF stereochemistry validation requires chain_boundary to be 'finite' or 'periodic'"
        )
    boundary = chain_boundary
    degree = len({int(atom.metadata["residue_number"]) for atom in molecule.atoms})
    priorities = _neighbour_priorities(polymer, boundary)
    index = {
        (int(atom.metadata["residue_number"]), atom.name): atom.molecule_atom_index
        for atom in molecule.atoms
    }
    centres: list[int] = []
    neighbours: list[list[int]] = []
    declared: list[str] = []
    keys: list[tuple[int, str]] = []
    for atom in molecule.atoms:
        if atom.stereochemistry is None:
            continue
        residue = int(atom.metadata["residue_number"])
        order = priorities.get((_residue_class(residue, degree, boundary), atom.name))
        if order is None:
            raise ValueError(
                f"no CIP priority reference for {polymer} {atom.name} in residue {residue}"
            )
        try:
            ranked = [
                index[_neighbour_residue(residue, offset, degree, boundary), name]
                for offset, name in order
            ]
        except KeyError as error:
            raise ValueError(
                f"OpenFF {polymer} molecular graph is incompatible with chain_boundary={boundary!r}"
            ) from error
        bonded = {neighbour.molecule_atom_index for neighbour in atom.bonded_atoms}
        if set(ranked) != bonded:
            raise ValueError(
                f"{polymer} {atom.name} in residue {residue}: CIP reference neighbours "
                "do not match the molecular graph"
            )
        centres.append(atom.molecule_atom_index)
        neighbours.append(ranked)
        declared.append(str(atom.stereochemistry))
        keys.append((residue, atom.name))
    chain_numbers = tuple(dict.fromkeys(mapped.chain_numbers))
    coordinates = np.asarray(mapped.cartesian_coordinates, dtype=np.float64).reshape(
        len(chain_numbers), molecule.n_atoms, 3
    )
    centre_xyz = coordinates[:, centres, None, :]
    neighbour_xyz = coordinates[:, np.asarray(neighbours), :]
    box = np.asarray(mapped.basis_vectors, dtype=np.float64)
    fractional = (neighbour_xyz - centre_xyz) @ np.linalg.inv(box)
    vectors = (fractional - np.round(fractional)) @ box
    handedness = np.einsum(
        "...i,...i->...", vectors[..., 0, :], np.cross(vectors[..., 1, :], vectors[..., 2, :])
    )
    if not np.all(np.isfinite(handedness)) or np.any(handedness == 0.0):
        raise ValueError(f"OpenFF {polymer} stereocentre geometry is degenerate")
    # Priorities 1-3 traversed clockwise with priority 4 pointing away is R,
    # which is a negative triple product of the three bond vectors.
    perceived = np.where(handedness < 0.0, "R", "S")
    failures: list[tuple[int, int, str, str, str]] = []
    signatures: set[tuple[str, ...]] = set()
    for chain_index, chain in enumerate(chain_numbers):
        labels = tuple(str(label) for label in perceived[chain_index])
        signatures.add(labels)
        for (residue, name), expected, actual in zip(keys, declared, labels, strict=True):
            if actual != expected:
                failures.append((chain, residue, name, expected, actual))
    if failures:
        raise ValueError(
            f"OpenFF {polymer} declared stereochemistry disagrees with mapped "
            f"coordinates: {failures[:5]}"
        )
    return {
        "method": (
            "RDKit CIP neighbour priorities from a short reference chain; handedness "
            "from minimum-image bond vectors of every centre"
        ),
        "chains_examined": len(chain_numbers),
        "distinct_chain_types_checked": len(signatures),
        "labelled_centres_per_chain": len(centres),
        "failures": failures,
        "passed": True,
    }


_REFERENCE_DEGREE = {"finite": 3, "periodic": 5}


def _residue_class(residue: int, degree: int, boundary: str) -> str:
    if boundary == "periodic":
        return "internal"
    if residue == 1:
        return "reducing"
    return "nonreducing" if residue == degree else "internal"


def _neighbour_residue(residue: int, offset: int, degree: int, boundary: str) -> int:
    target = residue + offset
    return (target - 1) % degree + 1 if boundary == "periodic" else target


@functools.cache
def _neighbour_priorities(
    polymer: str, boundary: str
) -> dict[tuple[str, str], tuple[tuple[int, str], ...]]:
    """Return each centre's neighbours in descending CIP priority, by residue kind.

    Neighbours are keyed by residue offset (-1, 0 or +1 along the chain) and atom
    name. RDKit ranks a short reference chain built from the same polymer graph;
    a ranking that cannot order all four neighbours is refused.
    """
    chem = importlib.import_module("rdkit.Chem")
    degree = _REFERENCE_DEGREE[boundary]
    reference = build_sage_molecule(degree, boundary, polymer)
    editable = chem.RWMol()
    for atom in reference.atoms:
        rdkit_atom = chem.Atom(int(atom.atomic_number))
        rdkit_atom.SetNoImplicit(True)
        editable.AddAtom(rdkit_atom)
    orders = {1: chem.BondType.SINGLE, 2: chem.BondType.DOUBLE}
    for bond in reference.bonds:
        editable.AddBond(
            int(bond.atom1.molecule_atom_index),
            int(bond.atom2.molecule_atom_index),
            orders[int(bond.bond_order)],
        )
    ranked = editable.GetMol()
    chem.SanitizeMol(ranked)
    chem.AssignStereochemistry(ranked, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    table: dict[tuple[str, str], tuple[tuple[int, str], ...]] = {}
    for atom in reference.atoms:
        if atom.stereochemistry is None:
            continue
        residue = int(atom.metadata["residue_number"])
        bonded = list(atom.bonded_atoms)
        rank = {
            neighbour.molecule_atom_index: int(
                ranked.GetAtomWithIdx(neighbour.molecule_atom_index).GetProp("_CIPRank")
            )
            for neighbour in bonded
        }
        if len(set(rank.values())) != len(bonded):
            raise ValueError(f"{polymer} {atom.name} neighbours cannot be CIP-ranked")
        order = []
        for neighbour in sorted(bonded, key=lambda item: -rank[item.molecule_atom_index]):
            offset = int(neighbour.metadata["residue_number"]) - residue
            if boundary == "periodic" and abs(offset) > 1:
                offset = -1 if offset > 0 else 1
            order.append((offset, neighbour.name))
        key = (_residue_class(residue, degree, boundary), atom.name)
        if table.setdefault(key, tuple(order)) != tuple(order):
            raise ValueError(f"{polymer} {atom.name} CIP order differs between equivalent residues")
    return table


def build_sage_system(
    crystal: BulkCrystal,
    chain_boundary: str,
    polymer: str | None = None,
) -> SageBuild:
    """Build, map, charge, and validate one reusable Sage polymer chain."""
    selected_polymer = polymer or ("chitin" if crystal.residue_name == "NAG" else "cellulose")
    profile = _polymer_profile(selected_polymer)
    if not crystal.chain_numbers or not crystal.residue_numbers:
        raise ValueError("OpenFF mapping requires explicit chain and residue numbers")
    first_chain = next(iter(crystal.chain_numbers))
    degree = len(
        {
            residue
            for chain, residue in zip(crystal.chain_numbers, crystal.residue_numbers, strict=True)
            if chain == first_chain
        }
    )
    molecule = build_sage_molecule(degree, chain_boundary, selected_polymer)
    mapped = map_crystal_to_chain_template(
        crystal,
        chain_template(
            (atom.name, int(atom.metadata["residue_number"])) for atom in molecule.atoms
        ),
        chain_boundary,
        profile.residue_name,
    )
    declared_stereo = validate_declared_stereochemistry(
        molecule, mapped, chain_boundary, selected_polymer
    )
    unit = importlib.import_module("openff.units").unit

    molecule.add_conformer(mapped.cartesian_coordinates[: molecule.n_atoms] * unit.angstrom)
    molecule.assign_partial_charges(NAGL_CHARGE_MODEL)
    charges = np.asarray(molecule.partial_charges.m, dtype=np.float64)
    if not np.isclose(charges.sum(), 0.0, atol=_CHARGE_TOLERANCE_E):
        raise ValueError(
            f"OpenFF NAGL {selected_polymer} chain charge is not zero: {charges.sum():.8g} e"
        )
    periodic_charge_check = True
    if chain_boundary == "periodic":
        by_name: dict[str, list[float]] = {}
        for atom, charge in zip(molecule.atoms, charges, strict=True):
            by_name.setdefault(atom.name, []).append(float(charge))
        periodic_charge_check = all(
            max(values) - min(values) <= _CHARGE_TOLERANCE_E for values in by_name.values()
        )
        if not periodic_charge_check:
            raise ValueError("periodic OpenFF charges differ between equivalent residues")
    stereo = validate_sage_stereochemistry(molecule, mapped, selected_polymer)
    stereo["declared_labels"] = declared_stereo
    charge_hash = hashlib.sha256(charges.tobytes()).hexdigest()
    provenance = {
        "method": NAGL_CHARGE_MODEL,
        "openff_toolkit": importlib.metadata.version("openff-toolkit"),
        "openff_nagl": importlib.metadata.version("openff-nagl"),
        "openff_nagl_models": importlib.metadata.version("openff-nagl-models"),
        "net_charge_e": float(charges.sum()),
        "periodic_equivalent_atom_charges": periodic_charge_check,
        "charge_sha256": charge_hash,
        "polymer": selected_polymer,
    }
    return SageBuild(molecule, mapped, degree, chain_boundary, provenance, stereo)


def write_sage_outputs(
    build: SageBuild,
    output_directory: Path,
    engines: tuple[str, ...],
    *,
    dispersion_correction: bool,
) -> dict[str, object]:
    """Write GROMACS and native serialized OpenMM outputs."""
    force_field_type, _, topology_type = require_openff()
    openmm = importlib.import_module("openmm")
    unit = importlib.import_module("openff.units").unit
    pdbx_file = importlib.import_module("openmm.app").PDBxFile

    output_directory.mkdir(parents=True, exist_ok=True)
    topology = topology_type.from_molecules([build.molecule] * build.mapped.chain_count)
    topology.box_vectors = reduced_periodic_box(build.mapped.basis_vectors) * unit.angstrom
    interchange = force_field_type(SAGE_FORCE_FIELD).create_interchange(
        topology, charge_from_molecules=[build.molecule]
    )
    interchange.positions = build.mapped.cartesian_coordinates * unit.angstrom
    write_gromacs = "gromacs" in engines
    if write_gromacs:
        prefix = str(output_directory / "structure")
        interchange.to_gromacs(prefix, decimal=8)
        (output_directory / "structure.top").replace(output_directory / "system.top")
        # Interchange writes this helper for its internal point-energy driver.
        # It is not part of the supported CelluloseForge output and its
        # defaults may differ from the geometry-aware inputs written below.
        (output_directory / "structure_pointenergy.mdp").unlink(missing_ok=True)

    periodic = build.chain_boundary == "periodic"
    system = interchange.to_openmm_system(
        combine_nonbonded_forces=False,
        add_constrained_forces=True,
        ewald_tolerance=5.0e-4,
    )
    periodic_forces: list[str] = []
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.setUseDispersionCorrection(dispersion_correction)
            if periodic:
                force.setExceptionsUsePeriodicBoundaryConditions(True)
                periodic_forces.append("NonbondedForce exceptions")
        elif isinstance(force, openmm.CustomNonbondedForce):
            force.setUseLongRangeCorrection(dispersion_correction)
        elif periodic and hasattr(force, "setUsesPeriodicBoundaryConditions"):
            force.setUsesPeriodicBoundaryConditions(True)
            periodic_forces.append(force.__class__.__name__)
    if "openmm" in engines:
        (output_directory / "system.xml").write_text(
            openmm.XmlSerializer.serialize(system), encoding="utf-8"
        )
        omm_topology = interchange.to_openmm_topology()
        with (output_directory / "topology.cif").open("w", encoding="utf-8") as handle:
            pdbx_file.writeFile(
                omm_topology,
                interchange.positions.to_openmm(),
                handle,
                keepIds=True,
            )
        # The starting state carries positions and box vectors; the runner is
        # the same configurable OpenMM runner used for GLYCAM and CHARMM36.
        validation_platform = validation_platform_for_box(build.mapped.basis_vectors)
        context = openmm.Context(
            system,
            openmm.VerletIntegrator(0.001),
            openmm.Platform.getPlatformByName(validation_platform),
        )
        context.setPeriodicBoxVectors(*interchange.box.to_openmm())
        context.setPositions(interchange.positions.to_openmm())
        (output_directory / OPENMM_STATE_FILE).write_text(
            openmm.XmlSerializer.serialize(context.getState(getPositions=True)),
            encoding="utf-8",
        )
        del context
        runner = output_directory / OPENMM_RUNNER_FILE
        runner.write_text(
            render_openmm_runner(
                chain_axis=chain_axis(build.mapped),
                periodic_chains=periodic,
                has_restraints=False,
            ),
            encoding="utf-8",
        )
        runner.chmod(0o755)
    outputs: dict[str, str] = {}
    if write_gromacs:
        outputs.update(
            {
                "gromacs_coordinates": "structure.gro",
                "gromacs_topology": "system.top",
                **write_gromacs_run_inputs(
                    output_directory,
                    force_field="openff",
                    periodic_chains=periodic,
                    dispersion=dispersion_correction,
                ),
            }
        )
    if "openmm" in engines:
        outputs.update(
            {
                "openmm_system": "system.xml",
                "openmm_state": OPENMM_STATE_FILE,
                "openmm_topology": "topology.cif",
                "openmm_runner": OPENMM_RUNNER_FILE,
            }
        )
    return {
        **outputs,
        "openmm_written": "openmm" in engines,
        "periodic_forces": sorted(set(periodic_forces)),
        "dispersion_correction": dispersion_correction,
        "validation_platform": validation_platform_for_box(build.mapped.basis_vectors),
    }


def write_sage_report(path: Path, report: dict[str, object]) -> None:
    """Write stable JSON for an experimental Sage build."""
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
