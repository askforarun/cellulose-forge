"""OpenMM and LAMMPS outputs for experimental OPLS-AA and OPLS-AA-SEI cellulose.

Both writers start from the same :class:`~cellulose_forge.opls.OplsTopology`
as the GROMACS writer, so all three engines share one parameter source.

OpenMM
    OpenMM's GROMACS reader rejects the function-2 pairs used for OPLS-AA-SEI,
    so the ``System`` is built directly.  Coulomb interactions use a PME
    ``NonbondedForce``; Lennard-Jones uses a ``CustomNonbondedForce`` with the
    OPLS geometric combination rule.  Every 1-4 and SEI pair is a
    ``NonbondedForce`` exception with its own charge product, sigma and
    epsilon, and is excluded from the custom force.

LAMMPS
    The data file is written directly (no OpenMM needed).  OPLS-AA uses
    ``pair_modify mix geometric``, ``special_bonds lj/coul 0 0 0.5`` and
    ``dihedral_style opls``, whose K1-K3 are Kony's V1-V3.  LAMMPS cannot
    scale individual 1-5 or 1-6 pairs, and ``bond_style special`` rejects
    long-range Coulomb, so each OPLS-AA-SEI pair is written as a tabulated
    pseudo-bond carrying its full Lennard-Jones plus scaled Coulomb energy.
    The pseudo-bond excludes the pair from the pair style, and the
    ``angle yes dihedral yes`` keywords of ``special_bonds`` keep the
    pseudo-bonds from creating spurious 1-3 or 1-4 exclusions.
"""

from __future__ import annotations

import json
import math
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from cellulose_forge.openmm import (
    OPENMM_RUNNER_FILE,
    OPENMM_STATE_FILE,
    OPENMM_SYSTEM_FILE,
    OPENMM_TOPOLOGY_FILE,
    _resolve_python,
    _validated_metadata,
    chain_axis,
    render_openmm_runner,
)
from cellulose_forge.opls import (
    _TORSIONS,
    KCAL_TO_KJ,
    OplsBuild,
    angle_parameters,
    atom_mass,
    bond_parameters,
    lj_pair,
    lj_parameters,
    nonbonded_exceptions,
    rb_coefficients,
)

_ANGSTROM_PER_NM = 10.0
# Coulomb constant of LAMMPS real units (force->qqr2e), kcal mol-1 angstrom e-2.
_LAMMPS_QQR2E = 332.06371
NONBONDED_CUTOFF_NM = 1.0
_SERIALIZED_EWALD_ERROR_TOLERANCE = 5.0e-4
_VALIDATION_EWALD_ERROR_TOLERANCE = 1.0e-5
# Tabulated SEI pseudo-bonds: range and number of points (angstrom).  The
# largest intra-residue 1-5/1-6 distance is about 6 angstrom; LAMMPS stops with
# an error if a pair ever leaves the table.  0.001 angstrom spacing.
SEI_TABLE_RANGE_ANGSTROM = (1.0, 8.0)
SEI_TABLE_POINTS = 7001
SEI_TABLE_EQUILIBRIUM_ANGSTROM = 4.0
SEI_TABLE_FILE = "sei_pairs.table"
LAMMPS_DATA_FILE = "system.data"
LAMMPS_INPUT_FILE = "in.cellulose"
_BOX_TILT_TOLERANCE = 1.0e-6


def validation_platform(build: OplsBuild) -> str:
    """OpenMM platform for the build-time energy check.

    The double-precision Reference platform is used unless the third box
    vector is tilted (cellulose I-alpha bulk crystals).  For such boxes the
    OpenMM 8.6.1 Reference Lennard-Jones energy disagreed with an exact
    nearest-image sum by 0.14--0.40 % in local measurements, while the CPU
    platform agreed within 1e-6 (validation/measure_openmm_tilted_box.py), so
    the CPU platform is used instead.
    """
    box = reduced_box_angstrom(build)
    tilted = abs(box[2][0]) > _BOX_TILT_TOLERANCE or abs(box[2][1]) > _BOX_TILT_TOLERANCE
    return "CPU" if tilted else "Reference"


def reduced_box_angstrom(build: OplsBuild) -> NDArray[np.float64]:
    """Return the box in the reduced form required by OpenMM and LAMMPS."""
    from cellulose_forge.openff_sage import reduced_periodic_box  # noqa: PLC0415

    return np.asarray(reduced_periodic_box(build.mapped.basis_vectors), dtype=np.float64)


def chain_template(build: OplsBuild) -> dict[str, object]:
    """Per-chain parameters in OpenMM units (nm, kJ/mol, radians)."""
    topology = build.topology
    torsions = _TORSIONS[topology.variant]
    hydrogen = {index for index, atom in enumerate(topology.atoms) if atom.element == "H"}
    exceptions = []
    for first, second, coulomb, lennard_jones, kind in nonbonded_exceptions(topology):
        sigma, epsilon = lj_pair(topology, first, second)
        exceptions.append(
            [
                first,
                second,
                coulomb * topology.charges[first] * topology.charges[second],
                sigma if lennard_jones else 1.0,
                lennard_jones * epsilon,
                kind,
            ]
        )
    return {
        "names": [atom.name for atom in topology.atoms],
        "elements": [atom.element for atom in topology.atoms],
        "residues": [atom.residue for atom in topology.atoms],
        "masses": [atom_mass(atom_type) for atom_type in topology.atom_types],
        "charges": list(topology.charges),
        "lj": [list(lj_parameters(atom_type)) for atom_type in topology.atom_types],
        "bonds": [
            [a, b, *bond_parameters(topology, a, b), a in hydrogen or b in hydrogen]
            for a, b in topology.bonds
        ],
        "angles": [
            [a, b, c, math.radians(theta), constant]
            for a, b, c in topology.angles
            for theta, constant in [angle_parameters(topology, a, b, c)]
        ],
        "torsions": [
            [a, b, c, d, *rb_coefficients(torsions[name])]
            for a, b, c, d, name in topology.dihedrals
        ],
        "exceptions": exceptions,
    }


_OPENMM_BUILDER = r"""
import json
import sys

import numpy
import openmm
from openmm import app, unit

request = json.loads(open(sys.argv[1], encoding="utf-8").read())
template = request["template"]
chains = request["chain_count"]
per_chain = len(template["masses"])
box = [openmm.Vec3(*row) for row in request["box_nm"]]
positions = [openmm.Vec3(*row) for row in request["positions_nm"]] * unit.nanometer
periodic = request["periodic_chains"]
cutoff = request["nonbonded_cutoff_nm"]
LJ = (
    "4*epsilon*((sigma/r)^12-(sigma/r)^6); "
    "sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)"
)


def make_system(constrain_hydrogen_bonds, ewald_error_tolerance):
    system = openmm.System()
    system.setDefaultPeriodicBoxVectors(*box)
    bonds = openmm.HarmonicBondForce()
    angles = openmm.HarmonicAngleForce()
    torsions = openmm.RBTorsionForce()
    coulomb = openmm.NonbondedForce()
    coulomb.setNonbondedMethod(openmm.NonbondedForce.PME)
    coulomb.setCutoffDistance(cutoff)
    coulomb.setEwaldErrorTolerance(ewald_error_tolerance)
    coulomb.setUseDispersionCorrection(False)
    lennard_jones = openmm.CustomNonbondedForce(LJ)
    lennard_jones.addPerParticleParameter("sigma")
    lennard_jones.addPerParticleParameter("epsilon")
    lennard_jones.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    lennard_jones.setCutoffDistance(cutoff)
    lennard_jones.setUseLongRangeCorrection(request["dispersion_correction"])
    for chain in range(chains):
        offset = chain * per_chain
        for mass, charge, (sigma, epsilon) in zip(
            template["masses"], template["charges"], template["lj"]
        ):
            system.addParticle(mass)
            coulomb.addParticle(charge, 1.0, 0.0)
            lennard_jones.addParticle([sigma, epsilon])
        for a, b, length, constant, hydrogen in template["bonds"]:
            if constrain_hydrogen_bonds and hydrogen:
                system.addConstraint(offset + a, offset + b, length)
            else:
                bonds.addBond(offset + a, offset + b, length, constant)
        for a, b, c, theta, constant in template["angles"]:
            angles.addAngle(offset + a, offset + b, offset + c, theta, constant)
        for a, b, c, d, *coefficients in template["torsions"]:
            torsions.addTorsion(offset + a, offset + b, offset + c, offset + d, *coefficients)
        for a, b, charge_product, sigma, epsilon, _ in template["exceptions"]:
            coulomb.addException(offset + a, offset + b, charge_product, sigma, epsilon)
            lennard_jones.addExclusion(offset + a, offset + b)
    enabled = []
    for name, force in (
        ("HarmonicBondForce", bonds),
        ("HarmonicAngleForce", angles),
        ("RBTorsionForce", torsions),
        ("NonbondedForce", coulomb),
        ("OPLSLennardJones", lennard_jones),
    ):
        force.setName(name)
        system.addForce(force)
    if periodic:
        coulomb.setExceptionsUsePeriodicBoundaryConditions(True)
        enabled.append("NonbondedForce exceptions")
        for force in (bonds, angles, torsions):
            force.setUsesPeriodicBoundaryConditions(True)
            enabled.append(force.getName())
    system.addForce(openmm.CMMotionRemover())
    return system, enabled


system, periodic_forces = make_system(False, request["validation_ewald_error_tolerance"])
names = []
for group, force in enumerate(system.getForces()):
    force.setForceGroup(group)
    names.append(force.getName() or force.__class__.__name__)
platform = openmm.Platform.getPlatformByName(request["validation_platform"])
context = openmm.Context(system, openmm.VerletIntegrator(0.001), platform)
context.setPeriodicBoxVectors(*box)
context.setPositions(positions)
initial_state = context.getState(getPositions=True, getEnergy=True)


def group_energy(group):
    state = context.getState(getEnergy=True, groups={group})
    return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


energies = {name: group_energy(group) for group, name in enumerate(names)}
coulomb_group = names.index("NonbondedForce")
coulomb = system.getForce(coulomb_group)
saved = [coulomb.getParticleParameters(index) for index in range(coulomb.getNumParticles())]
for index in range(coulomb.getNumParticles()):
    coulomb.setParticleParameters(index, 0.0, 1.0, 0.0)
coulomb.updateParametersInContext(context)
energies["NonbondedForce exceptions (1-4)"] = group_energy(coulomb_group)
for index, parameters in enumerate(saved):
    coulomb.setParticleParameters(index, *parameters)
coulomb.updateParametersInContext(context)

coordinates = numpy.asarray(request["positions_nm"])
box_matrix = numpy.asarray(request["box_nm"])
maximum_deviation = 0.0
bond_count = 0
for force in system.getForces():
    if isinstance(force, openmm.HarmonicBondForce):
        for index in range(force.getNumBonds()):
            first, second, length, _ = force.getBondParameters(index)
            delta = coordinates[second] - coordinates[first]
            for axis in (2, 1, 0):
                delta = delta - box_matrix[axis] * round(delta[axis] / box_matrix[axis][axis])
            deviation = abs(float(numpy.linalg.norm(delta)) - length.value_in_unit(unit.nanometer))
            maximum_deviation = max(maximum_deviation, deviation)
            bond_count += 1

constrained, _ = make_system(True, request["ewald_error_tolerance"])
constrained_context = openmm.Context(constrained, openmm.VerletIntegrator(0.001), platform)
constrained_context.setPeriodicBoxVectors(*box)
constrained_context.setPositions(positions)
openmm.LocalEnergyMinimizer.minimize(constrained_context, 10.0, request["minimization_iterations"])
minimized = constrained_context.getState(getEnergy=True).getPotentialEnergy()

with open(request["system_path"], "w", encoding="utf-8") as handle:
    handle.write(openmm.XmlSerializer.serialize(constrained))
with open(request["state_path"], "w", encoding="utf-8") as handle:
    handle.write(openmm.XmlSerializer.serialize(initial_state))
topology = app.Topology()
atoms = []
for chain_index in range(chains):
    chain = topology.addChain(id=str(chain_index + 1))
    residues = {}
    for name, element, residue_number in zip(
        template["names"], template["elements"], template["residues"]
    ):
        if residue_number not in residues:
            residues[residue_number] = topology.addResidue("BGC", chain, id=str(residue_number))
        atoms.append(topology.addAtom(
            name, app.Element.getBySymbol(element), residues[residue_number]
        ))
    offset = chain_index * per_chain
    for a, b, *_ in template["bonds"]:
        topology.addBond(atoms[offset + a], atoms[offset + b])
topology.setPeriodicBoxVectors(box)
with open(request["topology_path"], "w", encoding="utf-8") as handle:
    app.PDBxFile.writeFile(topology, positions, handle, keepIds=True)

print(json.dumps({
    "openmm_version": openmm.__version__,
    "validation_platform": request["validation_platform"],
    "particle_count": system.getNumParticles(),
    "constraint_count": constrained.getNumConstraints(),
    "validation_constraints": "none",
    "serialized_constraints": "HBonds",
    "nonbonded_cutoff_nm": cutoff,
    "lennard_jones": "CustomNonbondedForce, geometric combination rule, plain cutoff",
    "dispersion_correction": request["dispersion_correction"],
    "validation_ewald_error_tolerance": request["validation_ewald_error_tolerance"],
    "serialized_ewald_error_tolerance": request["ewald_error_tolerance"],
    "forces": names,
    "periodic_forces": periodic_forces,
    "initial_energies_kj_mol": energies,
    "initial_potential_energy_kj_mol": initial_state.getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    ),
    "minimized_potential_energy_kj_mol": minimized.value_in_unit(unit.kilojoule_per_mole),
    "minimization_iterations_requested": request["minimization_iterations"],
    "harmonic_bond_count": bond_count,
    "maximum_bond_length_deviation_nm": maximum_deviation,
}))
"""


@dataclass(frozen=True)
class OplsOpenMMFiles:
    system: Path
    state: Path
    topology: Path
    runner: Path


def write_opls_openmm(  # noqa: PLR0913
    build: OplsBuild,
    output_directory: Path,
    *,
    dispersion_correction: bool,
    python_executable: str | Path | None = None,
    minimization_iterations: int = 10,
    validation_ewald_error_tolerance: float = _VALIDATION_EWALD_ERROR_TOLERANCE,
) -> tuple[OplsOpenMMFiles, dict[str, object]]:
    """Build, check and serialize the OpenMM system for an OPLS build."""
    executable = _resolve_python(python_executable)
    output_directory.mkdir(parents=True, exist_ok=True)
    files = OplsOpenMMFiles(
        output_directory / OPENMM_SYSTEM_FILE,
        output_directory / OPENMM_STATE_FILE,
        output_directory / OPENMM_TOPOLOGY_FILE,
        output_directory / OPENMM_RUNNER_FILE,
    )
    periodic = build.topology.chain_boundary == "periodic"
    request = {
        "template": chain_template(build),
        "chain_count": build.mapped.chain_count,
        "box_nm": (reduced_box_angstrom(build) / _ANGSTROM_PER_NM).tolist(),
        "positions_nm": (
            np.asarray(build.mapped.cartesian_coordinates, dtype=np.float64) / _ANGSTROM_PER_NM
        ).tolist(),
        "periodic_chains": periodic,
        "nonbonded_cutoff_nm": NONBONDED_CUTOFF_NM,
        "dispersion_correction": dispersion_correction,
        "ewald_error_tolerance": _SERIALIZED_EWALD_ERROR_TOLERANCE,
        "validation_ewald_error_tolerance": validation_ewald_error_tolerance,
        "validation_platform": validation_platform(build),
        "minimization_iterations": minimization_iterations,
        "system_path": str(files.system.resolve()),
        "state_path": str(files.state.resolve()),
        "topology_path": str(files.topology.resolve()),
    }
    with tempfile.TemporaryDirectory(prefix="cellulose_forge_opls_openmm_") as temporary:
        request_path = Path(temporary) / "request.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        result = subprocess.run(
            [executable, "-c", _OPENMM_BUILDER, str(request_path)],
            check=False,
            capture_output=True,
            text=True,
        )
    if result.returncode != 0:
        raise RuntimeError(
            "OpenMM failed to build or validate the OPLS system "
            f"(exit code {result.returncode}):\n{result.stdout}\n{result.stderr}"
        )
    metadata = _validated_metadata(result.stdout, periodic_chains=periodic)
    axis = chain_axis(build.mapped)
    files.runner.write_text(
        render_openmm_runner(chain_axis=axis, periodic_chains=periodic, has_restraints=False),
        encoding="utf-8",
    )
    files.runner.chmod(0o755)
    metadata.update(
        {
            "python_executable": executable,
            "chain_axis": axis,
            "periodic_chains": periodic,
            "status": "passed",
        }
    )
    return files, metadata


@dataclass(frozen=True)
class OplsLammpsFiles:
    data: Path
    input: Path
    table: Path | None


def _wrapped_with_images(
    coordinates: NDArray[np.float64], box: NDArray[np.float64], bonds: list[tuple[int, int]]
) -> tuple[NDArray[np.float64], NDArray[np.int_]]:
    """Wrap atoms into the box and derive image flags from the bond graph."""
    fractional = coordinates @ np.linalg.inv(box)
    wrapped = fractional - np.floor(fractional)
    count = len(coordinates)
    neighbours: list[list[int]] = [[] for _ in range(count)]
    for first, second in bonds:
        neighbours[first].append(second)
        neighbours[second].append(first)
    images = np.zeros((count, 3), dtype=int)
    visited = np.zeros(count, dtype=bool)
    for root in range(count):
        if visited[root]:
            continue
        visited[root] = True
        stack = [root]
        while stack:
            first = stack.pop()
            for second in neighbours[first]:
                if not visited[second]:
                    # A periodic chain closes a cycle with a net image shift;
                    # LAMMPS evaluates each bond by minimum image regardless.
                    images[second] = images[first] - np.rint(
                        wrapped[second] - wrapped[first]
                    ).astype(int)
                    visited[second] = True
                    stack.append(second)
    return wrapped @ box, images


def _sei_table(
    sigma_nm: float, epsilon_kj: float, charge_product: float, scale: float
) -> list[str]:
    """Rows of one LAMMPS bond table: full LJ plus scaled Coulomb (real units)."""
    sigma = sigma_nm * _ANGSTROM_PER_NM
    epsilon = epsilon_kj / KCAL_TO_KJ
    coulomb = scale * _LAMMPS_QQR2E * charge_product
    low, high = SEI_TABLE_RANGE_ANGSTROM
    rows = []
    for index, r in enumerate(np.linspace(low, high, SEI_TABLE_POINTS), start=1):
        s6 = (sigma / r) ** 6
        energy = 4.0 * epsilon * (s6 * s6 - s6) + coulomb / r
        force = (24.0 * epsilon * (2.0 * s6 * s6 - s6) + coulomb / r) / r
        rows.append(f"{index} {r:.6f} {energy:.12e} {force:.12e}")
    return rows


_FloatKey = tuple[float, ...]


def _type_id(table: dict[_FloatKey, int], key: _FloatKey, offset: int = 0) -> int:
    return table.setdefault(key, offset + len(table) + 1)


@dataclass(frozen=True)
class _ChainRecords:
    """LAMMPS type tables and per-chain records (zero-based atom indices)."""

    atom_types: dict[str, int]
    harmonic_bonds: dict[_FloatKey, int]
    sei_tables: dict[_FloatKey, int]
    angle_types: dict[_FloatKey, int]
    dihedral_types: dict[_FloatKey, int]
    atom_type_of: list[int]
    bonds: list[tuple[int, ...]]
    angles: list[tuple[int, ...]]
    dihedrals: list[tuple[int, ...]]


def _chain_records(build: OplsBuild) -> _ChainRecords:
    topology = build.topology
    torsions = _TORSIONS[topology.variant]
    atom_types: dict[str, int] = {}
    atom_type_of = [
        atom_types.setdefault(name, len(atom_types) + 1) for name in topology.atom_types
    ]
    harmonic: dict[_FloatKey, int] = {}
    bonds: list[tuple[int, ...]] = [
        (_type_id(harmonic, bond_parameters(topology, a, b)), a, b) for a, b in topology.bonds
    ]
    # SEI pseudo-bond types follow the harmonic types.
    tables: dict[_FloatKey, int] = {}
    for a, b, scale, _ in topology.scaled_pairs:
        sigma, epsilon = lj_pair(topology, a, b)
        key = (
            round(sigma, 9), round(epsilon, 9),
            round(topology.charges[a] * topology.charges[b], 9), round(scale, 12),
        )
        bonds.append((_type_id(tables, key, len(harmonic)), a, b))
    angle_types: dict[_FloatKey, int] = {}
    angles: list[tuple[int, ...]] = [
        (_type_id(angle_types, angle_parameters(topology, a, b, c)), a, b, c)
        for a, b, c in topology.angles
    ]
    dihedral_types: dict[_FloatKey, int] = {}
    dihedrals: list[tuple[int, ...]] = [
        (_type_id(dihedral_types, torsions[name]), a, b, c, d)
        for a, b, c, d, name in topology.dihedrals
    ]
    return _ChainRecords(
        atom_types, harmonic, tables, angle_types, dihedral_types,
        atom_type_of, bonds, angles, dihedrals,
    )


def _coefficient_lines(records: _ChainRecords) -> list[str]:
    hybrid = bool(records.sei_tables)
    lines = ["", "Masses", ""]
    lines += [f"{i} {atom_mass(name):.6f} # {name}" for name, i in records.atom_types.items()]
    lines += ["", "Pair Coeffs # lj/cut/coul/long", ""]
    for name, index in records.atom_types.items():
        sigma, epsilon = lj_parameters(name)
        lines.append(f"{index} {epsilon / KCAL_TO_KJ:.12g} {sigma * _ANGSTROM_PER_NM:.12g}")
    lines += ["", "Bond Coeffs" + (" # hybrid" if hybrid else " # harmonic"), ""]
    style = "harmonic " if hybrid else ""
    for (length, constant), index in records.harmonic_bonds.items():
        # GROMACS E = kb/2 (b - b0)^2 in kJ mol-1 nm-2; LAMMPS E = K (r - r0)^2.
        lines.append(
            f"{index} {style}{constant / 2.0 / KCAL_TO_KJ / 100.0:.12g} {length * 10.0:.12g}"
        )
    for number, index in enumerate(records.sei_tables.values(), start=1):
        lines.append(f"{index} table {SEI_TABLE_FILE} SEI_{number}")
    lines += ["", "Angle Coeffs # harmonic", ""]
    for (theta, constant), index in records.angle_types.items():
        lines.append(f"{index} {constant / 2.0 / KCAL_TO_KJ:.12g} {theta:.12g}")
    # dihedral_style opls K1-K3 are Kony's V1-V3 in kcal/mol.
    lines += ["", "Dihedral Coeffs # opls", ""]
    for (v1, v2, v3), index in records.dihedral_types.items():
        lines.append(f"{index} {v1:.12g} {v2:.12g} {v3:.12g} 0.0")
    return lines


def _write_sei_table(path: Path, tables: dict[_FloatKey, int]) -> None:
    lines = [
        "# CelluloseForge OPLS-AA-SEI 1-5/1-6 pairs: full LJ plus scaled Coulomb",
        "# r (angstrom), E (kcal/mol), F (kcal/mol/angstrom)",
    ]
    for number, (sigma, epsilon, charge_product, scale) in enumerate(tables, start=1):
        # EQ replaces LAMMPS's default length estimate (the table energy
        # minimum, at the far end for repulsive pairs) used only for its
        # communication and box-size checks.
        header = f"N {SEI_TABLE_POINTS} EQ {SEI_TABLE_EQUILIBRIUM_ANGSTROM}"
        lines += ["", f"SEI_{number}", header, ""]
        lines += _sei_table(sigma, epsilon, charge_product, scale)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _lammps_input(variant: str, *, hybrid: bool, dispersion_correction: bool) -> str:
    bond_style = f"hybrid harmonic table linear {SEI_TABLE_POINTS}" if hybrid else "harmonic"
    tail = "yes" if dispersion_correction else "no"
    sei_note = (
        "# OPLS-AA-SEI 1-5/1-6 pairs are tabulated pseudo-bonds (sei_pairs.table):\n"
        "# their full Lennard-Jones and scaled Coulomb energies appear in ebond.\n"
        "# Expected, harmless LAMMPS warnings: table forces flagged at inflection\n"
        "# points, and 'Bond/angle/dihedral extent > half of periodic box length'\n"
        "# (a heuristic that treats a ~5 A pseudo-bond as part of a dihedral).\n"
        if hybrid
        else ""
    )
    return f"""# CelluloseForge experimental {variant} cellulose (single-point check)
units real
atom_style full
boundary p p p
bond_style {bond_style}
angle_style harmonic
dihedral_style opls
pair_style lj/cut/coul/long 10.0
pair_modify mix geometric tail {tail}
# 1-4 pairs at 0.5; 'angle yes dihedral yes' restricts the 1-3 and 1-4 lists to
# real angles and dihedrals.
special_bonds lj/coul 0.0 0.0 0.5 angle yes dihedral yes
{sei_note}read_data system.data
kspace_style pppm 1.0e-5
thermo 100
thermo_style custom step temp pe ebond eangle edihed evdwl ecoul elong
timestep 1.0
# This file is an initial-energy validation input, not a dynamics protocol.
# Covalently periodic chains form closed image cycles, so LAMMPS may warn
# about inconsistent image flags; bonds are evaluated by minimum image.
run 0
# Explicit termination avoids an exit-time fault in some packaged LAMMPS builds.
quit 0
"""


def write_opls_lammps(
    build: OplsBuild,
    output_directory: Path,
    *,
    dispersion_correction: bool,
) -> tuple[OplsLammpsFiles, dict[str, object]]:
    """Write ``system.data`` and ``in.cellulose`` (plus the SEI pair table)."""
    topology = build.topology
    per_chain = len(topology.atoms)
    chains = build.mapped.chain_count
    box = reduced_box_angstrom(build)
    if any(abs(value) > _BOX_TILT_TOLERANCE for value in (box[0][1], box[0][2], box[1][2])):
        raise ValueError("LAMMPS output requires a restricted triclinic box")
    records = _chain_records(build)
    covalent = [
        (chain * per_chain + a, chain * per_chain + b)
        for chain in range(chains)
        for a, b in topology.bonds
    ]
    positions, images = _wrapped_with_images(
        np.asarray(build.mapped.cartesian_coordinates, dtype=np.float64), box, covalent
    )
    bond_type_count = len(records.harmonic_bonds) + len(records.sei_tables)
    lines = [
        f"CelluloseForge experimental {topology.variant} cellulose", "",
        f"{chains * per_chain} atoms", f"{chains * len(records.bonds)} bonds",
        f"{chains * len(records.angles)} angles",
        f"{chains * len(records.dihedrals)} dihedrals", "",
        f"{len(records.atom_types)} atom types", f"{bond_type_count} bond types",
        f"{len(records.angle_types)} angle types",
        f"{len(records.dihedral_types)} dihedral types", "",
        f"0.0 {box[0][0]:.12g} xlo xhi", f"0.0 {box[1][1]:.12g} ylo yhi",
        f"0.0 {box[2][2]:.12g} zlo zhi",
        f"{box[1][0]:.12g} {box[2][0]:.12g} {box[2][1]:.12g} xy xz yz",
        *_coefficient_lines(records),
        "", "Atoms # full", "",
    ]
    for index, ((x, y, z), (ix, iy, iz)) in enumerate(zip(positions, images, strict=True)):
        atom = index % per_chain
        lines.append(
            f"{index + 1} {index // per_chain + 1} {records.atom_type_of[atom]} "
            f"{topology.charges[atom]:.6f} {x:.12g} {y:.12g} {z:.12g} {ix} {iy} {iz}"
        )
    for title, chain_records in (
        ("Bonds", records.bonds), ("Angles", records.angles), ("Dihedrals", records.dihedrals)
    ):
        lines += ["", title, ""]
        lines += [
            f"{chain * len(chain_records) + number} {kind} "
            + " ".join(str(member + chain * per_chain + 1) for member in members)
            for chain in range(chains)
            for number, (kind, *members) in enumerate(chain_records, start=1)
        ]

    hybrid = bool(records.sei_tables)
    output_directory.mkdir(parents=True, exist_ok=True)
    files = OplsLammpsFiles(
        output_directory / LAMMPS_DATA_FILE,
        output_directory / LAMMPS_INPUT_FILE,
        output_directory / SEI_TABLE_FILE if hybrid else None,
    )
    files.data.write_text("\n".join(lines) + "\n", encoding="utf-8")
    if files.table is not None:
        _write_sei_table(files.table, records.sei_tables)
    files.input.write_text(
        _lammps_input(
            topology.variant, hybrid=hybrid, dispersion_correction=dispersion_correction
        ),
        encoding="utf-8",
    )
    metadata: dict[str, object] = {
        "data": files.data.name,
        "input": files.input.name,
        "atom_count": chains * per_chain,
        "bond_count": chains * len(topology.bonds),
        "sei_pseudo_bond_count": chains * len(topology.scaled_pairs),
        "angle_count": chains * len(records.angles),
        "dihedral_count": chains * len(records.dihedrals),
        "atom_type_count": len(records.atom_types),
        "special_bonds": "lj/coul 0.0 0.0 0.5 angle yes dihedral yes",
        "dispersion_correction": dispersion_correction,
        "input_scope": "single-point run 0 validation; not a dynamics protocol",
    }
    if files.table is not None:
        metadata["sei_table"] = files.table.name
        metadata["sei_table_types"] = len(records.sei_tables)
    return files, metadata
