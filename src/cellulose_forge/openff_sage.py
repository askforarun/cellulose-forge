"""Experimental OpenFF Sage construction for crystalline cellulose.

OpenFF is imported lazily so the core CelluloseForge installation retains its
NumPy-only coordinate path.  Molecular graphs are constructed directly from
CelluloseForge atom names and bonds; no AmberTools or GLYCAM template is used.
"""

from __future__ import annotations

import hashlib
import importlib
import importlib.metadata
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from cellulose_forge.crystallography import FloatArray
from cellulose_forge.geometry import BulkCrystal
from cellulose_forge.glycam_mapping import GlycamMappedSystem
from cellulose_forge.hydroxyl import HydroxylHydrogenPlacer
from cellulose_forge.openmm import (
    OPENMM_RUNNER_FILE,
    OPENMM_STATE_FILE,
    chain_axis,
    render_openmm_runner,
)
from cellulose_forge.run_inputs import dispersion_correction as _dispersion_policy
from cellulose_forge.run_inputs import write_gromacs_run_inputs

SAGE_FORCE_FIELD = "openff-2.3.0.offxml"
NAGL_CHARGE_MODEL = "openff-gnn-am1bcc-1.0.0.pt"
VALIDATION_SCOPE = "model construction validated; force-field accuracy for cellulose not assessed"

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
_STEREOCHEMISTRY = {"C1": "S", "C2": "R", "C3": "R", "C4": "S", "C5": "R"}
_REFERENCE_PARITY = {"C1": 1, "C2": 1, "C3": -1, "C4": 1, "C5": 1}
_STEREOCENTRES = frozenset(_STEREOCHEMISTRY)
_MAX_LINK_ANGSTROM = 2.0
_MINIMUM_DP = 2
_CHARGE_TOLERANCE_E = 1.0e-6
_BOX_DEGENERACY_TOLERANCE = 1.0e-12
_MISSING_II_HYDROGENS = frozenset({"HO2", "HO3", "HO6"})

_INTERNAL_REFERENCE = {
    "C1": (2.720, 11.610, 12.966),
    "H1": (2.053, 11.447, 12.268),
    "O5": (3.942, 12.030, 12.413),
    "C5": (3.939, 13.328, 11.848),
    "H5": (3.382, 13.347, 11.042),
    "C6": (5.429, 13.601, 11.483),
    "H62": (5.837, 12.780, 11.164),
    "H61": (5.907, 13.890, 12.273),
    "O6": (5.535, 14.593, 10.478),
    "HO6": (6.378, 14.431, 10.000),
    "C4": (3.454, 14.380, 12.860),
    "H4": (4.094, 14.469, 13.594),
    "C3": (2.078, 13.975, 13.392),
    "H3": (1.435, 13.946, 12.653),
    "O3": (1.602, 14.853, 14.397),
    "HO3": (1.714, 15.767, 14.053),
    "C2": (2.194, 12.595, 14.004),
    "H2": (2.841, 12.639, 14.743),
    "O2": (0.944, 12.153, 14.530),
    "HO2": (1.042, 11.190, 14.695),
    "O4": (3.368, 15.606, 12.110),
}
_REDUCING_REFERENCE = {
    "HO1": (3.163, 4.521, 12.573),
    "O1": (3.479, 5.230, 12.007),
    "C1": (3.658, 6.420, 12.839),
    "H1": (4.325, 6.257, 13.537),
    "O5": (2.436, 6.840, 13.392),
}
_NONREDUCING_REFERENCE = {
    "O4": (3.368, 36.366, 12.110),
    "C4": (3.454, 35.140, 12.860),
    "H4": (4.094, 35.229, 13.594),
    "HO4": (2.949, 37.032, 12.679),
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


def reduced_periodic_box(box: FloatArray) -> FloatArray:
    """Return an equivalent triclinic lattice in OpenMM reduced form."""
    a, b, c = np.asarray(box, dtype=np.float64).copy()
    if (
        abs(a[0]) < _BOX_DEGENERACY_TOLERANCE
        or abs(b[1]) < _BOX_DEGENERACY_TOLERANCE
    ):
        raise ValueError("OpenFF output requires a non-degenerate triclinic box")
    c -= b * round(float(c[1] / b[1]))
    c -= a * round(float(c[0] / a[0]))
    b -= a * round(float(b[0] / a[0]))
    reduced = np.vstack((a, b, c))
    reduced.setflags(write=False)
    return reduced


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


def _lookup(crystal: BulkCrystal, chain: int, residue: int) -> dict[str, FloatArray]:
    return {
        name: xyz
        for name, xyz, atom_chain, atom_residue in zip(
            crystal.atom_labels,
            crystal.cartesian_coordinates,
            crystal.chain_numbers,
            crystal.residue_numbers,
            strict=True,
        )
        if atom_chain == chain and atom_residue == residue
    }


def _minimum_image(delta: FloatArray, box: FloatArray) -> FloatArray:
    fractional = delta @ np.linalg.inv(box)
    fractional -= np.round(fractional)
    return np.asarray(fractional @ box, dtype=np.float64)


def _residue_order(crystal: BulkCrystal, chain: int, periodic: bool) -> tuple[int, ...]:
    residues = tuple(
        dict.fromkeys(
            residue
            for atom_chain, residue in zip(
                crystal.chain_numbers, crystal.residue_numbers, strict=True
            )
            if atom_chain == chain
        )
    )
    lookups = {residue: _lookup(crystal, chain, residue) for residue in residues}
    successors: dict[int, int] = {}
    for residue in residues:
        candidates = []
        for other in residues:
            if other == residue:
                continue
            delta = lookups[other]["C1"] - lookups[residue]["O4"]
            if periodic:
                delta = _minimum_image(delta, crystal.basis_vectors)
            candidates.append((float(np.linalg.norm(delta)), other))
        distance, other = min(candidates)
        if distance <= _MAX_LINK_ANGSTROM:
            successors[residue] = other
    if periodic:
        start = residues[0]
    else:
        predecessors = set(successors.values())
        starts = [residue for residue in residues if residue not in predecessors]
        if len(starts) != 1:
            raise ValueError(f"chain {chain} does not form one finite cellulose path")
        start = starts[0]
    order = [start]
    while len(order) < len(residues) and order[-1] in successors:
        following = successors[order[-1]]
        if following in order:
            break
        order.append(following)
    if len(order) != len(residues):
        raise ValueError(f"chain {chain} cellulose connectivity is incomplete")
    if periodic and successors.get(order[-1]) != order[0]:
        raise ValueError(f"chain {chain} cellulose connectivity is not cyclic")
    return tuple(order)


def _rigid_transform(
    reference: dict[str, tuple[float, float, float]],
    target: dict[str, FloatArray],
    anchors: tuple[str, str, str],
) -> tuple[FloatArray, FloatArray, FloatArray]:
    source = np.vstack([reference[name] for name in anchors])
    destination = np.vstack([target[name] for name in anchors])
    source_center = source.mean(axis=0)
    destination_center = destination.mean(axis=0)
    left, _, right = np.linalg.svd((source - source_center).T @ (destination - destination_center))
    rotation = left @ right
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right
    return rotation, source_center, destination_center


def _transform(
    coordinate: tuple[float, float, float],
    transform: tuple[FloatArray, FloatArray, FloatArray],
) -> FloatArray:
    rotation, source, destination = transform
    return np.asarray(
        (np.asarray(coordinate) - source) @ rotation + destination,
        dtype=np.float64,
    )


def build_sage_molecule(  # noqa: PLR0912
    degree_polymerization: int, chain_boundary: str
) -> Any:
    """Construct a finite or cyclic beta-(1->4)-D-glucan graph directly."""
    if degree_polymerization < _MINIMUM_DP:
        raise ValueError("OpenFF cellulose chains require DP >= 2")
    if chain_boundary not in {"finite", "periodic"}:
        raise ValueError("OpenFF chain boundary must be finite or periodic")
    _, molecule_type, _ = require_openff()
    molecule = molecule_type()
    molecule.name = f"cellulose_dp{degree_polymerization}_{chain_boundary}"
    indices: dict[tuple[int, str], int] = {}
    for residue in range(1, degree_polymerization + 1):
        for name, atomic_number in _INTERNAL_ATOMS:
            stereo = _STEREOCHEMISTRY.get(name)
            if chain_boundary == "finite" and residue == 1 and name == "C1":
                stereo = "R"
            if chain_boundary == "finite" and residue == degree_polymerization and name == "C3":
                stereo = "S"
            indices[residue, name] = molecule.add_atom(
                atomic_number,
                0,
                False,
                stereo,
                name,
                {"residue_name": "BGC", "residue_number": residue, "chain_id": "A"},
            )
        if chain_boundary == "finite" and residue == 1:
            for name, atomic_number in (("O1", 8), ("HO1", 1)):
                indices[residue, name] = molecule.add_atom(
                    atomic_number,
                    0,
                    False,
                    None,
                    name,
                    {"residue_name": "BGC", "residue_number": residue, "chain_id": "A"},
                )
        if chain_boundary == "finite" and residue == degree_polymerization:
            indices[residue, "HO4"] = molecule.add_atom(
                1,
                0,
                False,
                None,
                "HO4",
                {"residue_name": "BGC", "residue_number": residue, "chain_id": "A"},
            )
        for left, right in _INTERNAL_BONDS:
            molecule.add_bond(indices[residue, left], indices[residue, right], 1, False)
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


def _mapped_coordinates(crystal: BulkCrystal, molecule: Any, boundary: str) -> GlycamMappedSystem:
    periodic = boundary == "periodic"
    chain_numbers = tuple(dict.fromkeys(crystal.chain_numbers))
    output: list[FloatArray] = []
    names: list[str] = []
    residues: list[int] = []
    chains: list[int] = []
    placer = (
        HydroxylHydrogenPlacer(
            crystal.atom_labels, crystal.cartesian_coordinates, crystal.basis_vectors
        )
        if crystal.allomorph == "II"
        else None
    )
    for chain in chain_numbers:
        order = _residue_order(crystal, chain, periodic)
        targets = {
            index + 1: _lookup(crystal, chain, residue) for index, residue in enumerate(order)
        }
        reducing = (
            _rigid_transform(_REDUCING_REFERENCE, targets[1], ("C1", "H1", "O5"))
            if not periodic
            else None
        )
        nonreducing = (
            _rigid_transform(
                _NONREDUCING_REFERENCE,
                targets[len(order)],
                ("O4", "C4", "H4"),
            )
            if not periodic
            else None
        )
        hydrogen_transforms = (
            {
                residue: _rigid_transform(_INTERNAL_REFERENCE, target, ("C1", "C3", "C5"))
                for residue, target in targets.items()
            }
            if crystal.allomorph == "II"
            else {}
        )
        for atom in molecule.atoms:
            residue = int(atom.metadata["residue_number"])
            name = atom.name
            target = targets[residue]
            if name in {"O1", "HO1"}:
                assert reducing is not None
                coordinate = _transform(_REDUCING_REFERENCE[name], reducing)
            elif name == "HO4":
                assert nonreducing is not None
                coordinate = _transform(_NONREDUCING_REFERENCE[name], nonreducing)
            elif crystal.allomorph == "II" and name in _MISSING_II_HYDROGENS:
                oxygen = "O" + name[2]
                carbon = "C" + name[2]
                transform = hydrogen_transforms[residue]
                aligned_h = _transform(_INTERNAL_REFERENCE[name], transform)
                assert placer is not None
                placer.request(
                    len(output),
                    target[carbon],
                    target[oxygen],
                    aligned_h,
                    donor=oxygen,
                    preferred_acceptor="O5" if oxygen == "O3" else None,
                )
                coordinate = np.zeros(3)
            else:
                coordinate = target[name]
            output.append(coordinate)
            names.append(name)
            residues.append(residue)
            chains.append(chain)
    if placer is not None:
        for index, coordinate in placer.resolve().items():
            output[index] = coordinate
    coordinates = np.vstack(output)
    coordinates.setflags(write=False)
    return GlycamMappedSystem(
        atom_names=tuple(names),
        residue_names=tuple("BGC" for _ in names),
        residue_numbers=tuple(residues),
        chain_numbers=tuple(chains),
        cartesian_coordinates=coordinates,
        basis_vectors=crystal.basis_vectors,
    )


def _centre_parity(positions: FloatArray, atom: Any, box: FloatArray) -> int:
    neighbours = sorted(atom.bonded_atoms, key=lambda item: (item.symbol, item.name))[:3]
    centre = positions[atom.molecule_atom_index]
    vectors = [
        _minimum_image(positions[item.molecule_atom_index] - centre, box) for item in neighbours
    ]
    value = float(np.dot(vectors[0], np.cross(vectors[1], vectors[2])))
    if not np.isfinite(value) or value == 0.0:
        raise ValueError("undefined cellulose stereocentre parity")
    return int(np.sign(value))


def validate_sage_stereochemistry(molecule: Any, mapped: GlycamMappedSystem) -> dict[str, object]:
    """Check every mapped C1--C5 centre against an independent reference."""
    failures: list[tuple[int, int, str]] = []
    for chain_index, chain in enumerate(dict.fromkeys(mapped.chain_numbers)):
        start = chain_index * molecule.n_atoms
        positions = mapped.cartesian_coordinates[start : start + molecule.n_atoms]
        for atom in molecule.atoms:
            if (
                atom.name in _STEREOCENTRES
                and _centre_parity(positions, atom, mapped.basis_vectors)
                != _REFERENCE_PARITY[atom.name]
            ):
                failures.append(  # noqa: PERF401
                    (chain, int(atom.metadata["residue_number"]), atom.name)
                )
    result = {
        "reference": "independent PubChem finite DP6 internal-residue signed parity",
        "centres_checked": mapped.chain_count
        * len({int(atom.metadata["residue_number"]) for atom in molecule.atoms})
        * len(_STEREOCENTRES),
        "failures": failures,
        "passed": not failures,
    }
    if failures:
        raise ValueError(f"OpenFF cellulose stereochemistry check failed: {failures[:5]}")
    return result


def build_sage_system(crystal: BulkCrystal, chain_boundary: str) -> SageBuild:
    """Build, map, charge, and validate one reusable Sage cellulose chain."""
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
    molecule = build_sage_molecule(degree, chain_boundary)
    mapped = _mapped_coordinates(crystal, molecule, chain_boundary)
    unit = importlib.import_module("openff.units").unit

    molecule.add_conformer(mapped.cartesian_coordinates[: molecule.n_atoms] * unit.angstrom)
    molecule.assign_partial_charges(NAGL_CHARGE_MODEL)
    charges = np.asarray(molecule.partial_charges.m, dtype=np.float64)
    if not np.isclose(charges.sum(), 0.0, atol=_CHARGE_TOLERANCE_E):
        raise ValueError(f"OpenFF NAGL chain charge is not zero: {charges.sum():.8g} e")
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
    stereo = validate_sage_stereochemistry(molecule, mapped)
    charge_hash = hashlib.sha256(charges.tobytes()).hexdigest()
    provenance = {
        "method": NAGL_CHARGE_MODEL,
        "openff_toolkit": importlib.metadata.version("openff-toolkit"),
        "openff_nagl": importlib.metadata.version("openff-nagl"),
        "openff_nagl_models": importlib.metadata.version("openff-nagl-models"),
        "net_charge_e": float(charges.sum()),
        "periodic_equivalent_atom_charges": periodic_charge_check,
        "charge_sha256": charge_hash,
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
        context = openmm.Context(
            system,
            openmm.VerletIntegrator(0.001),
            openmm.Platform.getPlatformByName("Reference"),
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
    }


def write_sage_report(path: Path, report: dict[str, object]) -> None:
    """Write stable JSON for an experimental Sage build."""
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
