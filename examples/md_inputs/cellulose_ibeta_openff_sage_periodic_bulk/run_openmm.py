#!/usr/bin/env python3
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

CHAIN_AXIS = 'z'
PERIODIC_CHAINS = True
HAS_RESTRAINTS = False
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
        simulation.context.setParameter("cellulose_restraint_scale", 0.0)
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
