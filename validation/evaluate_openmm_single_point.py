#!/usr/bin/env python3
"""Evaluate one generated OpenMM validation system with canonical energy drivers."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import openmm
from openmm.app import GromacsTopFile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cellulose_forge.energy_drivers import (  # noqa: E402
    openmm_classification_from_topology,
    run_openmm_single_point,
)


def main() -> int:
    """Classify a cellulose system explicitly and print its canonical report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", type=Path, required=True)
    parser.add_argument("--coordinates", type=Path, required=True)
    parser.add_argument("--coordinate-reference", type=Path, required=True)
    parser.add_argument("--coordinate-sha256", required=True)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--include-directory", type=Path, required=True)
    parser.add_argument("--platform", default="Reference")
    arguments = parser.parse_args()

    topology = GromacsTopFile(
        str(arguments.topology), includeDir=str(arguments.include_directory)
    ).topology
    covalent_bonds = {(first.index, second.index) for first, second in topology.bonds()}
    system = openmm.XmlSerializer.deserialize(arguments.system.read_text(encoding="utf-8"))
    periodic_torsions = {
        index: ()
        for index, force in enumerate(system.getForces())
        if isinstance(force, openmm.PeriodicTorsionForce)
    }
    custom_roles = {
        index: "harmonic_improper"
        for index, force in enumerate(system.getForces())
        if isinstance(force, openmm.CustomTorsionForce)
    }
    classification = openmm_classification_from_topology(
        arguments.system,
        covalent_bonds=covalent_bonds,
        periodic_impropers=periodic_torsions,
        force_roles=custom_roles,
    )
    report = run_openmm_single_point(
        system_file=arguments.system,
        coordinate_file=arguments.coordinates,
        coordinate_reference_file=arguments.coordinate_reference,
        coordinate_sha256=arguments.coordinate_sha256,
        classification=classification,
        platform_name=arguments.platform,
    )
    print(json.dumps(report.to_dict(), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
