"""Serialized OpenMM systems for GLYCAM and CHARMM36 cellulose builds.

OpenMM is driven through a separate Python interpreter so that the core package
keeps NumPy as its only runtime dependency. The builder reads the GROMACS-format
topology that CelluloseForge has already written, enables periodic treatment of
every bonded term and of nonbonded exceptions for covalently periodic chains,
checks the initial interaction energies, and serializes the complete
``System`` and starting ``State``. The generated runner needs only OpenMM.
"""

from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import numpy as np

from cellulose_forge.glycam_mapping import GlycamMappedSystem

OPENMM_SYSTEM_FILE = "system.xml"
OPENMM_STATE_FILE = "state.xml"
OPENMM_TOPOLOGY_FILE = "topology.cif"
OPENMM_RUNNER_FILE = "run_openmm.py"
POSITION_RESTRAINT_SCALE_PARAMETER = "cellulose_restraint_scale"

_ANGSTROM_TO_NANOMETER = 0.1
_AXIS_NAMES = ("x", "y", "z")
# Matches the initial bonded-energy criterion of the GROMACS integration tests.
_MAXIMUM_BONDED_ENERGY_PER_ATOM_KJ_MOL = 50.0
# A box-length bond across a periodic boundary deviates by nanometres.
_MAXIMUM_BOND_DEVIATION_NM = 0.05
# OpenMM's default PME accuracy for simulation, and GROMACS's default ewald-rtol
# for the energy evaluation used to compare engines.
_SERIALIZED_EWALD_ERROR_TOLERANCE = 5.0e-4
_VALIDATION_EWALD_ERROR_TOLERANCE = 1.0e-5


@dataclass(frozen=True)
class OpenMMPositionRestraints:
    """Heavy-atom harmonic restraints expressed in OpenMM units."""

    atom_indices: tuple[int, ...]
    force_constants_kj_mol_nm2: tuple[float, float, float]


@dataclass(frozen=True)
class OpenMMSystemFiles:
    """Paths of a serialized OpenMM system written by CelluloseForge."""

    system: Path
    state: Path
    topology: Path
    runner: Path


_BUILDER = r"""
import json
import sys

import numpy
import openmm
from openmm import unit
from openmm.app import HBonds, PME, GromacsGroFile, GromacsTopFile, PDBxFile

request = json.loads(open(sys.argv[1], encoding="utf-8").read())
gro = GromacsGroFile(request["gromacs_coordinates"])
box_vectors = gro.getPeriodicBoxVectors()
top = GromacsTopFile(
    request["gromacs_topology"],
    periodicBoxVectors=box_vectors,
    includeDir=request["include_directory"],
)
periodic = request["periodic_chains"]
restraints = request["position_restraints"]


def make_system(constraints, ewald_error_tolerance):
    created = top.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=request["nonbonded_cutoff_nm"] * unit.nanometer,
        constraints=constraints,
        ewaldErrorTolerance=ewald_error_tolerance,
    )
    if created.getNumParticles() != request["expected_atom_count"]:
        raise SystemExit(
            f"OpenMM read {created.getNumParticles()} particles; "
            f"expected {request['expected_atom_count']}"
        )
    enabled = []
    for force in created.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.setUseDispersionCorrection(request["dispersion_correction"])
            if periodic:
                force.setExceptionsUsePeriodicBoundaryConditions(True)
                enabled.append("NonbondedForce exceptions")
        elif periodic and hasattr(force, "setUsesPeriodicBoundaryConditions"):
            force.setUsesPeriodicBoundaryConditions(True)
            enabled.append(force.__class__.__name__)
    if restraints is not None:
        restraint = openmm.CustomExternalForce(
            "0.5*cellulose_restraint_scale*(kx*(x-x0)^2+ky*(y-y0)^2+kz*(z-z0)^2)"
        )
        restraint.addGlobalParameter("cellulose_restraint_scale", 1.0)
        for name in ("x0", "y0", "z0", "kx", "ky", "kz"):
            restraint.addPerParticleParameter(name)
        kx, ky, kz = restraints["force_constants_kj_mol_nm2"]
        for index in restraints["atom_indices"]:
            x0, y0, z0 = gro.positions[index].value_in_unit(unit.nanometer)
            restraint.addParticle(index, [x0, y0, z0, kx, ky, kz])
        restraint.setName("CellulosePositionRestraints")
        created.addForce(restraint)
    return created, enabled


# Interaction energies are checked without constraints so that every bonded term,
# including X-H bonds, is evaluated, and with a tight PME tolerance so that the
# nonbonded total can be compared with other engines. The serialized system
# constrains X-H bonds and uses the production PME tolerance.
system, periodic_forces = make_system(None, request["validation_ewald_error_tolerance"])
names = []
for group, force in enumerate(system.getForces()):
    force.setForceGroup(group)
    names.append(force.getName() or force.__class__.__name__)

platform = openmm.Platform.getPlatformByName("Reference")
context = openmm.Context(system, openmm.VerletIntegrator(0.001), platform)
context.setPeriodicBoxVectors(*box_vectors)
context.setPositions(gro.positions)
initial_state = context.getState(getPositions=True, getEnergy=True)

def group_energy(group):
    state = context.getState(getEnergy=True, groups={group})
    return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

energies = {}
for group, name in enumerate(names):
    energies[name] = energies.get(name, 0.0) + group_energy(group)

nonbonded_group = next(
    group for group, force in enumerate(system.getForces())
    if isinstance(force, openmm.NonbondedForce)
)
nonbonded = system.getForce(nonbonded_group)
saved_parameters = [
    nonbonded.getParticleParameters(index) for index in range(nonbonded.getNumParticles())
]
for index in range(nonbonded.getNumParticles()):
    nonbonded.setParticleParameters(index, 0.0, 1.0, 0.0)
nonbonded.updateParametersInContext(context)
energies["NonbondedForce exceptions (1-4)"] = group_energy(nonbonded_group)
for index, (charge, sigma, epsilon) in enumerate(saved_parameters):
    nonbonded.setParticleParameters(index, charge, sigma, epsilon)
nonbonded.updateParametersInContext(context)

positions = numpy.asarray(
    initial_state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
)
box = numpy.asarray(
    [[component for component in vector.value_in_unit(unit.nanometer)]
     for vector in box_vectors]
)

def minimum_image(delta):
    for axis in (2, 1, 0):
        delta = delta - box[axis] * round(delta[axis] / box[axis][axis])
    return delta

maximum_deviation = 0.0
bond_count = 0
for force in system.getForces():
    if isinstance(force, openmm.HarmonicBondForce):
        for index in range(force.getNumBonds()):
            first, second, length, _ = force.getBondParameters(index)
            delta = minimum_image(positions[second] - positions[first])
            distance = float(numpy.linalg.norm(delta))
            deviation = abs(distance - length.value_in_unit(unit.nanometer))
            maximum_deviation = max(maximum_deviation, deviation)
            bond_count += 1

initial = initial_state.getPotentialEnergy()
constrained, _ = make_system(HBonds, request["ewald_error_tolerance"])
constrained_context = openmm.Context(
    constrained, openmm.VerletIntegrator(0.001), platform
)
constrained_context.setPeriodicBoxVectors(*box_vectors)
constrained_context.setPositions(gro.positions)
openmm.LocalEnergyMinimizer.minimize(
    constrained_context, 10.0, request["minimization_iterations"]
)
minimized = constrained_context.getState(getEnergy=True).getPotentialEnergy()

with open(request["system_path"], "w", encoding="utf-8") as handle:
    handle.write(openmm.XmlSerializer.serialize(constrained))
with open(request["state_path"], "w", encoding="utf-8") as handle:
    handle.write(openmm.XmlSerializer.serialize(initial_state))
top.topology.setPeriodicBoxVectors(box_vectors)
with open(request["topology_path"], "w", encoding="utf-8") as handle:
    PDBxFile.writeFile(top.topology, gro.positions, handle, keepIds=True)

print(json.dumps({
    "openmm_version": openmm.__version__,
    "validation_platform": "Reference",
    "particle_count": system.getNumParticles(),
    "constraint_count": constrained.getNumConstraints(),
    "validation_constraints": "none",
    "serialized_constraints": "HBonds",
    "nonbonded_cutoff_nm": request["nonbonded_cutoff_nm"],
    "dispersion_correction": request["dispersion_correction"],
    "validation_ewald_error_tolerance": request["validation_ewald_error_tolerance"],
    "serialized_ewald_error_tolerance": request["ewald_error_tolerance"],
    "forces": names,
    "periodic_forces": periodic_forces,
    "initial_energies_kj_mol": energies,
    "initial_potential_energy_kj_mol": initial.value_in_unit(unit.kilojoule_per_mole),
    "minimized_potential_energy_kj_mol": minimized.value_in_unit(
        unit.kilojoule_per_mole
    ),
    "minimization_iterations_requested": request["minimization_iterations"],
    "harmonic_bond_count": bond_count,
    "maximum_bond_length_deviation_nm": maximum_deviation,
}))
"""


_RUNNER = '''#!/usr/bin/env python3
"""Run a CelluloseForge cellulose system with OpenMM.

The serialized system already contains every force-field term, the periodic
treatment required by covalently periodic chains, and any position restraints.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import openmm
from openmm import unit
from openmm.app import DCDReporter, PDBxFile, Simulation, StateDataReporter

CHAIN_AXIS = __CHAIN_AXIS__
PERIODIC_CHAINS = __PERIODIC_CHAINS__
HAS_RESTRAINTS = __HAS_RESTRAINTS__
AXES = "xyz"
BOX_TILT_TOLERANCE = 1e-9


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", help="OpenMM platform name, e.g. CUDA or CPU")
    parser.add_argument("--temperature", type=float, default=300.0, help="kelvin")
    parser.add_argument("--friction", type=float, default=1.0, help="1/ps")
    parser.add_argument("--timestep-fs", type=float, default=2.0)
    parser.add_argument("--minimize-iterations", type=int, default=1000)
    parser.add_argument("--steps", type=int, default=0)
    parser.add_argument("--report-interval", type=int, default=1000)
    parser.add_argument(
        "--no-restraints",
        action="store_true",
        help="switch off generated cellulose position restraints",
    )
    parser.add_argument(
        "--pressure-bar",
        type=float,
        help="enable an anisotropic Monte Carlo barostat at this pressure",
    )
    parser.add_argument(
        "--pressure-axes",
        default=None,
        help="Cartesian box axes scaled by the barostat (default: all but the chain axis)",
    )
    args = parser.parse_args()

    with Path("system.xml").open(encoding="utf-8") as handle:
        system = openmm.XmlSerializer.deserialize(handle.read())
    with Path("state.xml").open(encoding="utf-8") as handle:
        state = openmm.XmlSerializer.deserialize(handle.read())
    topology = PDBxFile("topology.cif").topology

    if args.pressure_bar is not None:
        axes = args.pressure_axes or "".join(
            axis for axis in AXES if axis != CHAIN_AXIS
        )
        if HAS_RESTRAINTS and not args.no_restraints:
            parser.error(
                "position restraints use fixed reference coordinates that a barostat "
                "does not rescale; add --no-restraints for pressure coupling"
            )
        if PERIODIC_CHAINS and CHAIN_AXIS in axes:
            parser.error(
                "covalently periodic chains cross the box along the chain axis; "
                "OpenMM barostats move molecules rigidly, so that axis must stay fixed"
            )
        box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
        if any(abs(float(box[row][column])) > BOX_TILT_TOLERANCE for row in range(3)
               for column in range(3) if row != column):
            parser.error("the anisotropic barostat requires a rectangular box")
        system.addForce(
            openmm.MonteCarloAnisotropicBarostat(
                openmm.Vec3(args.pressure_bar, args.pressure_bar, args.pressure_bar)
                * unit.bar,
                args.temperature * unit.kelvin,
                "x" in axes,
                "y" in axes,
                "z" in axes,
            )
        )

    integrator = openmm.LangevinMiddleIntegrator(
        args.temperature * unit.kelvin,
        args.friction / unit.picosecond,
        args.timestep_fs * unit.femtosecond,
    )
    platform = (
        openmm.Platform.getPlatformByName(args.platform) if args.platform else None
    )
    simulation = Simulation(topology, system, integrator, platform)
    simulation.context.setState(state)
    if HAS_RESTRAINTS and args.no_restraints:
        simulation.context.setParameter("__RESTRAINT_PARAMETER__", 0.0)
    if args.minimize_iterations:
        simulation.minimizeEnergy(maxIterations=args.minimize_iterations)
    if args.steps:
        simulation.context.setVelocitiesToTemperature(args.temperature * unit.kelvin)
        simulation.reporters.append(DCDReporter("trajectory.dcd", args.report_interval))
        for target in (sys.stdout, "state_data.csv"):
            simulation.reporters.append(
                StateDataReporter(
                    target,
                    args.report_interval,
                    step=True,
                    time=True,
                    potentialEnergy=True,
                    temperature=True,
                    volume=True,
                    speed=True,
                )
            )
        simulation.step(args.steps)
    final = simulation.context.getState(getPositions=True, getEnergy=True)
    with Path("final_state.xml").open("w", encoding="utf-8") as handle:
        handle.write(openmm.XmlSerializer.serialize(final))
    with Path("final.cif").open("w", encoding="utf-8") as handle:
        PDBxFile.writeFile(topology, final.getPositions(), handle, keepIds=True)
    energy = final.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    print(f"Final potential energy: {energy:.6f} kJ/mol")


if __name__ == "__main__":
    main()
'''


def _resolve_python(executable: str | Path | None) -> str:
    """Return the Python interpreter used to load OpenMM."""
    requested = str(executable or sys.executable)
    resolved = shutil.which(requested)
    if resolved is not None:
        return resolved
    path = Path(requested)
    if path.is_file():
        return str(path)
    raise FileNotFoundError(f"Python executable for OpenMM was not found: {requested}")


def chain_axis(system: GlycamMappedSystem) -> str:
    """Return the Cartesian axis most closely aligned with the cellulose chains."""
    first_chain = system.chain_numbers[0]
    anomeric = sorted(
        (residue, coordinates)
        for name, residue, chain, coordinates in zip(
            system.atom_names,
            system.residue_numbers,
            system.chain_numbers,
            np.asarray(system.cartesian_coordinates),
            strict=True,
        )
        if chain == first_chain and name == "C1"
    )
    if len(anomeric) < 2:  # noqa: PLR2004
        raise ValueError("cannot determine the chain axis from fewer than two C1 atoms")
    direction = anomeric[-1][1] - anomeric[0][1]
    return _AXIS_NAMES[int(np.argmax(np.abs(direction)))]


def heavy_atom_restraints(
    system: GlycamMappedSystem,
    chain_ids: tuple[int, ...],
    force_constants_kj_mol_nm2: tuple[float, float, float],
) -> OpenMMPositionRestraints:
    """Select zero-based heavy-atom indices of the restrained chains."""
    selected = set(chain_ids)
    indices = tuple(
        index
        for index, (name, chain) in enumerate(
            zip(system.atom_names, system.chain_numbers, strict=True)
        )
        if chain in selected and not name.upper().startswith("H")
    )
    return OpenMMPositionRestraints(indices, force_constants_kj_mol_nm2)


def render_openmm_runner(*, chain_axis: str, periodic_chains: bool, has_restraints: bool) -> str:
    """Return the standalone runner text for one generated system."""
    return (
        _RUNNER.replace("__CHAIN_AXIS__", repr(chain_axis))
        .replace("__PERIODIC_CHAINS__", repr(periodic_chains))
        .replace("__HAS_RESTRAINTS__", repr(has_restraints))
        .replace("__RESTRAINT_PARAMETER__", POSITION_RESTRAINT_SCALE_PARAMETER)
    )


def _check_box_supports_cutoff(system: GlycamMappedSystem, cutoff_nm: float) -> None:
    """Reject boxes whose perpendicular width is below twice the PME cutoff."""
    vectors = np.asarray(system.basis_vectors, dtype=np.float64) * _ANGSTROM_TO_NANOMETER
    volume = abs(float(np.linalg.det(vectors)))
    widths = [
        volume / float(np.linalg.norm(np.cross(vectors[(axis + 1) % 3], vectors[(axis + 2) % 3])))
        for axis in range(3)
    ]
    narrowest = min(widths)
    if narrowest < 2.0 * cutoff_nm:
        raise ValueError(
            f"the periodic box is only {narrowest:.2f} nm wide, less than twice the "
            f"{cutoff_nm:g} nm nonbonded cutoff required by OpenMM; build a larger "
            "system or add vacuum"
        )


def _validated_metadata(stdout: str, *, periodic_chains: bool) -> dict[str, object]:
    """Parse and check the builder's JSON report."""
    try:
        metadata = cast(dict[str, object], json.loads(stdout.strip().splitlines()[-1]))
    except (IndexError, json.JSONDecodeError) as error:
        raise RuntimeError(f"could not parse OpenMM builder output: {stdout}") from error
    for key in (
        "initial_potential_energy_kj_mol",
        "minimized_potential_energy_kj_mol",
        "maximum_bond_length_deviation_nm",
    ):
        value = metadata.get(key)
        if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
            raise ValueError(f"OpenMM returned a non-finite {key}: {value}")
    energies = cast(dict[str, float], metadata.get("initial_energies_kj_mol", {}))
    for name, value in energies.items():
        if not math.isfinite(float(value)):
            raise ValueError(f"OpenMM returned a non-finite {name} energy: {value}")
    deviation = float(cast(float, metadata["maximum_bond_length_deviation_nm"]))
    if deviation > _MAXIMUM_BOND_DEVIATION_NM:
        raise ValueError(
            f"a bonded distance deviates by {deviation:.3f} nm from its equilibrium "
            "length in OpenMM; a covalent term is not using the correct periodic image"
        )
    bonded = sum(
        float(value)
        for name, value in energies.items()
        if name in {"HarmonicBondForce", "HarmonicAngleForce", "PeriodicTorsionForce"}
    )
    particles = int(cast(int, metadata["particle_count"]))
    if bonded / particles > _MAXIMUM_BONDED_ENERGY_PER_ATOM_KJ_MOL:
        raise ValueError(
            f"initial OpenMM bonded energy is {bonded / particles:.1f} kJ/mol per atom, "
            f"above the {_MAXIMUM_BONDED_ENERGY_PER_ATOM_KJ_MOL:g} kJ/mol limit"
        )
    if periodic_chains:
        periodic_forces = cast(list[str], metadata.get("periodic_forces", []))
        if "NonbondedForce exceptions" not in periodic_forces:
            raise ValueError("periodic nonbonded exceptions were not enabled in OpenMM")
    return metadata


def write_openmm_system(  # noqa: PLR0913
    system: GlycamMappedSystem,
    gromacs_topology: str | Path,
    gromacs_coordinates: str | Path,
    output_directory: str | Path,
    *,
    periodic_chains: bool,
    dispersion_correction: bool,
    restraints: OpenMMPositionRestraints | None = None,
    include_directory: str | Path | None = None,
    python_executable: str | Path | None = None,
    minimization_iterations: int = 10,
    nonbonded_cutoff_nm: float = 1.0,
) -> tuple[OpenMMSystemFiles, dict[str, object]]:
    """Serialize a validated OpenMM system, state, topology, and runner.

    ``dispersion_correction`` should be enabled only for homogeneous periodic
    systems; slabs and fibrils contain vacuum, which the isotropic long-range
    correction does not describe.
    """
    _check_box_supports_cutoff(system, nonbonded_cutoff_nm)
    executable = _resolve_python(python_executable)
    topology_path = Path(gromacs_topology).resolve()
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=True)
    files = OpenMMSystemFiles(
        system=directory / OPENMM_SYSTEM_FILE,
        state=directory / OPENMM_STATE_FILE,
        topology=directory / OPENMM_TOPOLOGY_FILE,
        runner=directory / OPENMM_RUNNER_FILE,
    )
    request = {
        "gromacs_topology": str(topology_path),
        "gromacs_coordinates": str(Path(gromacs_coordinates).resolve()),
        "include_directory": str(
            Path(include_directory).resolve()
            if include_directory is not None
            else topology_path.parent
        ),
        "expected_atom_count": system.atom_count,
        "periodic_chains": periodic_chains,
        "position_restraints": (
            None
            if restraints is None
            else {
                "atom_indices": list(restraints.atom_indices),
                "force_constants_kj_mol_nm2": list(restraints.force_constants_kj_mol_nm2),
            }
        ),
        "nonbonded_cutoff_nm": nonbonded_cutoff_nm,
        "dispersion_correction": dispersion_correction,
        "ewald_error_tolerance": _SERIALIZED_EWALD_ERROR_TOLERANCE,
        "validation_ewald_error_tolerance": _VALIDATION_EWALD_ERROR_TOLERANCE,
        "minimization_iterations": minimization_iterations,
        "system_path": str(files.system.resolve()),
        "state_path": str(files.state.resolve()),
        "topology_path": str(files.topology.resolve()),
    }
    with tempfile.TemporaryDirectory(prefix="cellulose_forge_openmm_") as temporary:
        request_path = Path(temporary) / "request.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        result = subprocess.run(
            [executable, "-c", _BUILDER, str(request_path)],
            check=False,
            capture_output=True,
            text=True,
        )
    if result.returncode != 0:
        raise RuntimeError(
            "OpenMM failed to build or validate the generated system "
            f"(exit code {result.returncode}):\n{result.stdout}\n{result.stderr}"
        )
    metadata = _validated_metadata(result.stdout, periodic_chains=periodic_chains)
    axis = chain_axis(system)
    files.runner.write_text(
        render_openmm_runner(
            chain_axis=axis,
            periodic_chains=periodic_chains,
            has_restraints=restraints is not None,
        ),
        encoding="utf-8",
    )
    files.runner.chmod(0o755)
    metadata.update(
        {
            "python_executable": executable,
            "chain_axis": axis,
            "periodic_chains": periodic_chains,
            "position_restraint_atom_count": (
                0 if restraints is None else len(restraints.atom_indices)
            ),
            "position_restraint_parameter": (
                None if restraints is None else POSITION_RESTRAINT_SCALE_PARAMETER
            ),
            "status": "passed",
        }
    )
    return files, metadata
