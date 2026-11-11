"""OpenMM output and engine-level validation for GLYCAM systems."""

from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
from pathlib import Path
from typing import cast

_RUNNER = '''#!/usr/bin/env python3
"""Run a generated CelluloseForge GLYCAM system with OpenMM."""

from __future__ import annotations

import argparse
import sys

from openmm import LangevinMiddleIntegrator
from openmm import unit
from openmm.app import (
    AmberInpcrdFile,
    AmberPrmtopFile,
    DCDReporter,
    HBonds,
    PME,
    PDBFile,
    Simulation,
    StateDataReporter,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--timestep-fs", type=float, default=2.0)
    parser.add_argument("--minimize-iterations", type=int, default=1000)
    parser.add_argument("--steps", type=int, default=0)
    parser.add_argument("--report-interval", type=int, default=1000)
    args = parser.parse_args()

    topology = AmberPrmtopFile("structure.prmtop")
    coordinates = AmberInpcrdFile("structure.inpcrd")
    system = topology.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
    )
    integrator = LangevinMiddleIntegrator(
        args.temperature * unit.kelvin,
        1.0 / unit.picosecond,
        args.timestep_fs * unit.femtosecond,
    )
    simulation = Simulation(topology.topology, system, integrator)
    simulation.context.setPositions(coordinates.positions)
    if coordinates.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*coordinates.boxVectors)
    simulation.minimizeEnergy(maxIterations=args.minimize_iterations)
    if args.steps:
        simulation.context.setVelocitiesToTemperature(args.temperature * unit.kelvin)
        simulation.reporters.append(DCDReporter("trajectory.dcd", args.report_interval))
        simulation.reporters.append(
            StateDataReporter(
                sys.stdout,
                args.report_interval,
                step=True,
                potentialEnergy=True,
                temperature=True,
                speed=True,
            )
        )
        simulation.step(args.steps)
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    with open("openmm_final.pdb", "w", encoding="utf-8") as handle:
        PDBFile.writeFile(topology.topology, state.getPositions(), handle)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    print(f"Final potential energy: {energy:.6f} kJ/mol")


if __name__ == "__main__":
    main()
'''


_VALIDATION_CODE = '''
import json
import sys

import openmm
from openmm import LangevinMiddleIntegrator, unit
from openmm.app import AmberInpcrdFile, AmberPrmtopFile, HBonds, PME, Simulation

topology = AmberPrmtopFile(sys.argv[1])
coordinates = AmberInpcrdFile(sys.argv[2])
system = topology.createSystem(
    nonbondedMethod=PME,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=HBonds,
)
integrator = LangevinMiddleIntegrator(
    300.0 * unit.kelvin,
    1.0 / unit.picosecond,
    2.0 * unit.femtosecond,
)
simulation = Simulation(topology.topology, system, integrator)
simulation.context.setPositions(coordinates.positions)
if coordinates.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*coordinates.boxVectors)
initial = simulation.context.getState(getEnergy=True).getPotentialEnergy()
simulation.minimizeEnergy(maxIterations=int(sys.argv[3]))
final = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(json.dumps({
    "openmm_version": openmm.__version__,
    "platform": simulation.context.getPlatform().getName(),
    "particle_count": system.getNumParticles(),
    "force_count": system.getNumForces(),
    "initial_potential_energy_kj_mol": initial.value_in_unit(unit.kilojoule_per_mole),
    "minimized_potential_energy_kj_mol": final.value_in_unit(unit.kilojoule_per_mole),
    "minimization_iterations_requested": int(sys.argv[3]),
    "status": "passed",
}))
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


def write_openmm_runner(output_directory: str | Path) -> Path:
    """Write a standalone OpenMM minimization and dynamics runner."""
    path = Path(output_directory) / "run_openmm.py"
    path.write_text(_RUNNER, encoding="utf-8")
    path.chmod(0o755)
    return path


def validate_openmm_system(
    topology_path: str | Path,
    coordinate_path: str | Path,
    *,
    python_executable: str | Path | None = None,
    minimization_iterations: int = 10,
) -> dict[str, object]:
    """Load and briefly minimize an Amber-format system in OpenMM."""
    executable = _resolve_python(python_executable)
    result = subprocess.run(
        [
            executable,
            "-c",
            _VALIDATION_CODE,
            str(Path(topology_path).resolve()),
            str(Path(coordinate_path).resolve()),
            str(minimization_iterations),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "OpenMM failed to load or minimize the generated system "
            f"(exit code {result.returncode}):\n{result.stdout}\n{result.stderr}"
        )
    try:
        metadata = cast(
            dict[str, object], json.loads(result.stdout.strip().splitlines()[-1])
        )
    except (IndexError, json.JSONDecodeError) as error:
        raise RuntimeError(
            f"could not parse OpenMM validation output: {result.stdout}"
        ) from error
    for key in (
        "initial_potential_energy_kj_mol",
        "minimized_potential_energy_kj_mol",
    ):
        value = metadata.get(key)
        if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
            raise ValueError(f"OpenMM returned a non-finite {key}: {value}")
    metadata["python_executable"] = executable
    return metadata
