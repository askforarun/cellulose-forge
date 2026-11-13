"""Probe OpenFF Sage coverage of finite beta-cellobiose.

This is an experimental validation script, not a CelluloseForge backend.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tempfile
from importlib import metadata
from pathlib import Path

import numpy as np
import openmm
import parmed
from openff.interchange.drivers.gromacs import get_gromacs_energies
from openff.interchange.drivers.openmm import get_openmm_energies
from openff.toolkit import ForceField, Molecule
from openff.units import unit

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from cellulose_forge.crystallography import MillerIndex, load_allomorph  # noqa: E402
from cellulose_forge.geometry import build_chain_parallel_slab  # noqa: E402
from cellulose_forge.glycam import generate_glycam_gromacs_reference  # noqa: E402
from cellulose_forge.glycam_mapping import map_crystal_to_finite_glycam  # noqa: E402

CELLODEXTRIN_SMILES = {
    2: (
        "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O[C@@H]2"
        "[C@H](O[C@H]([C@@H]([C@H]2O)O)O)CO)O)O)O)O"
    ),
    6: (
        "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O[C@@H]2[C@H](O[C@H]"
        "([C@@H]([C@H]2O)O)O[C@@H]3[C@H](O[C@H]([C@@H]([C@H]3O)O)O"
        "[C@@H]4[C@H](O[C@H]([C@@H]([C@H]4O)O)O[C@@H]5[C@H](O[C@H]"
        "([C@@H]([C@H]5O)O)O[C@@H]6[C@H](O[C@H]([C@@H]([C@H]6O)O)O)"
        "CO)CO)CO)CO)CO)O)O)O)O"
    ),
}
CELLODEXTRIN_METADATA = {
    2: ("beta-D-cellobiose", "PubChem CID 10712 (beta-cellobiose)"),
    6: ("beta-D-cellohexaose", "PubChem CID 440948 (beta-D-cellohexaose)"),
}
CHARGE_METHODS = {
    "am1bcc": {
        "openff_method": "am1bcc",
        "backend": "AmberToolsToolkitWrapper",
        "conformer_behavior": (
            "OpenFF generated one charge conformer internally; it was not the retained "
            "CelluloseForge crystallographic conformer and is not stored by the toolkit"
        ),
    },
    "nagl": {
        "openff_method": "openff-gnn-am1bcc-1.0.0.pt",
        "backend": "NAGLToolkitWrapper",
        "conformer_behavior": "Graph-based and independent of molecular conformation",
    },
}

ENERGY_TOLERANCES_KJ_MOL = {
    "Bond": 0.01,
    "Angle": 0.01,
    "Torsion": 0.01,
    "vdW 1-4": 0.01,
    "Electrostatics 1-4": 0.01,
}
ENERGY_TOLERANCES_KJ_MOL_PER_ATOM = {
    "vdW": 0.01,
    "Electrostatics": 0.0015,
}
MAXIMUM_NVT_BOND_CHANGE_ANGSTROM = 0.2


def energy_values(report: object) -> dict[str, float]:
    """Return the available energy terms in kJ/mol."""
    energies = report.energies  # type: ignore[attr-defined]
    return {
        name: float(value.m_as(unit.kilojoule_per_mole))
        for name, value in energies.items()
        if value is not None
    }


def bond_lengths(positions: np.ndarray, molecule: Molecule) -> np.ndarray:
    """Measure all molecular bond lengths in angstrom."""
    return np.asarray(
        [
            np.linalg.norm(positions[bond.atom1_index] - positions[bond.atom2_index])
            for bond in molecule.bonds
        ]
    )


def run_nvt_smoke_test(
    system: openmm.System,
    molecule: Molecule,
    minimized_state: openmm.State,
    minimized_positions: np.ndarray,
) -> dict[str, float | int | bool]:
    """Run 1 ps NVT and verify that molecular bonds remain stable."""
    minimized_bond_lengths = bond_lengths(minimized_positions, molecule)
    integrator = openmm.LangevinMiddleIntegrator(
        300.0 * openmm.unit.kelvin,
        1.0 / openmm.unit.picosecond,
        0.001 * openmm.unit.picoseconds,
    )
    context = openmm.Context(system, integrator)
    context.setPositions(minimized_state.getPositions())
    context.setVelocitiesToTemperature(300.0 * openmm.unit.kelvin, 20261112)
    integrator.step(1000)
    state = context.getState(getEnergy=True, getPositions=True, getVelocities=True)
    positions = state.getPositions(asNumpy=True).value_in_unit(openmm.unit.angstrom)
    maximum_bond_change = float(
        np.max(np.abs(bond_lengths(positions, molecule) - minimized_bond_lengths))
    )
    return {
        "duration_ps": 1.0,
        "temperature_k": 300.0,
        "time_step_fs": 1.0,
        "random_seed": 20261112,
        "final_potential_energy_kj_mol": state.getPotentialEnergy().value_in_unit(
            openmm.unit.kilojoule_per_mole
        ),
        "final_kinetic_energy_kj_mol": state.getKineticEnergy().value_in_unit(
            openmm.unit.kilojoule_per_mole
        ),
        "maximum_bond_length_change_angstrom": maximum_bond_change,
        "passed": bool(
            np.isfinite(positions).all()
            and np.isfinite(maximum_bond_change)
            and maximum_bond_change < MAXIMUM_NVT_BOND_CHANGE_ANGSTROM
        ),
    }


def crystallographic_cellodextrin(directory: Path, degree_polymerization: int) -> Molecule:
    """Map one finite I-beta chain into an explicitly graphed OpenFF molecule."""
    crystal = build_chain_parallel_slab(
        load_allomorph("i_beta"),
        MillerIndex(1, 1, 0),
        layers=1,
        chains_per_layer=2,
        degree_polymerization=degree_polymerization,
        vacuum=20.0,
        padding_z=12.0,
        chain_boundary="finite",
    )
    files, _ = generate_glycam_gromacs_reference(
        degree_polymerization=degree_polymerization,
        chain_boundary="finite",
        directory=directory,
        parmed_python=sys.executable,
    )
    mapped = map_crystal_to_finite_glycam(crystal, files.pdb)
    reference = parmed.load_file(str(files.prmtop))
    reference.coordinates = mapped.cartesian_coordinates[: len(reference.atoms)]
    pdb = directory / f"crystallographic_dp{degree_polymerization}.pdb"
    reference.save(str(pdb), overwrite=True)
    # Topology.from_pdb cannot yet recognize the GLYCAM ROH/4GB/0GB
    # carbohydrate residue split, even when supplied the complete molecule as
    # a unique template. Use the explicit graph loader until OpenFF's supported
    # PDB route handles multi-residue carbohydrates.
    return Molecule.from_pdb_and_smiles(
        pdb,
        CELLODEXTRIN_SMILES[degree_polymerization],
        allow_undefined_stereo=False,
        name=f"crystallographic_cellodextrin_dp{degree_polymerization}",
    )


def main() -> None:  # noqa: PLR0915 - validation stages are intentionally explicit
    parser = argparse.ArgumentParser()
    parser.add_argument("--dp", type=int, choices=sorted(CELLODEXTRIN_SMILES), default=2)
    parser.add_argument(
        "--charge-method",
        choices=sorted(CHARGE_METHODS),
        default="am1bcc",
    )
    args = parser.parse_args()
    degree_polymerization = args.dp
    charge_method = args.charge_method
    cache_infix = "" if charge_method == "am1bcc" else f"_{charge_method}"
    charge_cache = Path(
        "validation/results/"
        f"openff_cellodextrin_dp{degree_polymerization}{cache_infix}_charged.json"
    )
    if charge_cache.exists():
        molecule = Molecule.from_json(charge_cache.read_text(encoding="utf-8"))
    else:
        with tempfile.TemporaryDirectory(
            prefix=f"cellulose_forge_crystal_dp{degree_polymerization}_"
        ) as temporary:
            molecule = crystallographic_cellodextrin(
                Path(temporary), degree_polymerization
            )
        openff_method = CHARGE_METHODS[charge_method]["openff_method"]
        if charge_method == "am1bcc":
            molecule.assign_partial_charges(openff_method, use_conformers=None)
        else:
            molecule.assign_partial_charges(openff_method)
        charge_cache.write_text(molecule.to_json(indent=2) + "\n", encoding="utf-8")

    graph_digest = hashlib.sha256(
        CELLODEXTRIN_SMILES[degree_polymerization].encode("utf-8")
    ).hexdigest()
    coordinate_digest = hashlib.sha256(
        np.asarray(molecule.conformers[0].m_as(unit.angstrom), dtype=np.float64).tobytes()
    ).hexdigest()
    cache_digest = hashlib.sha256(charge_cache.read_bytes()).hexdigest()
    charge_provenance = {
        "method": CHARGE_METHODS[charge_method]["openff_method"],
        "backend": CHARGE_METHODS[charge_method]["backend"],
        "use_conformers": None if charge_method == "am1bcc" else "not applicable",
        "conformer_behavior": CHARGE_METHODS[charge_method]["conformer_behavior"],
        "openff_toolkit_version": metadata.version("openff-toolkit"),
        "openff_nagl_version": (
            metadata.version("openff-nagl") if charge_method == "nagl" else "not used"
        ),
        "openff_nagl_models_version": (
            metadata.version("openff-nagl-models")
            if charge_method == "nagl"
            else "not used"
        ),
        "ambertools_version": (
            "26.0 (Conda environment pin)"
            if charge_method == "am1bcc"
            else "not used for NAGL charging"
        ),
        "graph_sha256": graph_digest,
        "crystallographic_coordinate_sha256": coordinate_digest,
        "charged_cache_sha256": cache_digest,
    }

    force_field = ForceField("openff-2.3.0.offxml")
    topology = molecule.to_topology()
    labels = force_field.label_molecules(topology)[0]
    system = force_field.create_openmm_system(topology)

    context = openmm.Context(system, openmm.VerletIntegrator(0.001))
    context.setPositions(molecule.conformers[0].to_openmm())
    initial_energy = context.getState(getEnergy=True).getPotentialEnergy()
    openmm.LocalEnergyMinimizer.minimize(context, tolerance=10.0, maxIterations=2000)
    minimized_state = context.getState(getEnergy=True, getPositions=True)
    minimized_energy = minimized_state.getPotentialEnergy()
    initial_positions = molecule.conformers[0].m_as(unit.angstrom)
    minimized_positions = minimized_state.getPositions(asNumpy=True).value_in_unit(
        openmm.unit.angstrom
    )
    displacement = minimized_positions - initial_positions
    del context
    nvt_result = run_nvt_smoke_test(system, molecule, minimized_state, minimized_positions)

    # Interchange/GROMACS requires box vectors even for an isolated molecule.
    # A 4 nm cube avoids interactions with periodic images in this export probe.
    topology.box_vectors = np.eye(3) * 4.0 * unit.nanometer
    interchange = force_field.create_interchange(
        topology,
        charge_from_molecules=[molecule],
    )
    interchange.positions = molecule.conformers[0]
    openmm_report = get_openmm_energies(
        interchange,
        detailed=True,
        combine_nonbonded_forces=False,
    )
    gromacs_report = get_gromacs_energies(interchange, detailed=True)
    openmm_energies = energy_values(openmm_report)
    gromacs_energies = energy_values(gromacs_report)
    common_terms = sorted(openmm_energies.keys() & gromacs_energies.keys())
    energy_tolerances = {
        **ENERGY_TOLERANCES_KJ_MOL,
        **{
            term: tolerance * molecule.n_atoms
            for term, tolerance in ENERGY_TOLERANCES_KJ_MOL_PER_ATOM.items()
        },
    }
    energy_comparison = {
        term: {
            "openmm_kj_mol": openmm_energies[term],
            "gromacs_kj_mol": gromacs_energies[term],
            "difference_kj_mol": openmm_energies[term] - gromacs_energies[term],
            "relative_difference": (
                abs(openmm_energies[term] - gromacs_energies[term])
                / max(abs(openmm_energies[term]), abs(gromacs_energies[term]), 1.0)
            ),
            "tolerance_kj_mol": energy_tolerances[term],
            "passed": abs(openmm_energies[term] - gromacs_energies[term])
            <= energy_tolerances[term],
        }
        for term in common_terms
    }
    missing_terms = sorted(energy_tolerances.keys() - set(common_terms))
    if missing_terms:
        raise RuntimeError(f"Energy comparison omitted required terms: {missing_terms}")
    failed_terms = [
        term for term, comparison in energy_comparison.items() if not comparison["passed"]
    ]
    with tempfile.TemporaryDirectory(prefix="cellulose_forge_openff_") as temporary:
        directory = Path(temporary)
        prefix = str(directory / f"cellodextrin_dp{degree_polymerization}")
        interchange.to_gromacs(prefix)
        interchange.to_amber(prefix)
        exported_files = {
            path.name: path.stat().st_size for path in sorted(directory.iterdir())
        }

    result = {
        "molecule": CELLODEXTRIN_METADATA[degree_polymerization][0],
        "degree_polymerization": degree_polymerization,
        "coordinates": f"CelluloseForge I-beta (110), finite DP{degree_polymerization} chain",
        "smiles_source": CELLODEXTRIN_METADATA[degree_polymerization][1],
        "openff_force_field": "openff-2.3.0.offxml",
        "charge_method": charge_method,
        "charge_provenance": charge_provenance,
        "atom_count": molecule.n_atoms,
        "bond_count": molecule.n_bonds,
        "formal_charge_e": molecule.total_charge.m_as(unit.elementary_charge),
        "partial_charge_e": molecule.partial_charges.sum().m_as(unit.elementary_charge),
        "parameter_counts": {name: len(values) for name, values in labels.items()},
        "openmm_forces": [type(force).__name__ for force in system.getForces()],
        "cross_engine_single_point": {
            "openmm": openmm_energies,
            "gromacs": gromacs_energies,
            "comparison": energy_comparison,
            "nonbonded_tolerance_policy_kj_mol_per_atom": (
                ENERGY_TOLERANCES_KJ_MOL_PER_ATOM
            ),
            "passed": not failed_terms,
            "failed_terms": failed_terms,
        },
        "initial_energy_kj_mol": initial_energy.value_in_unit(
            openmm.unit.kilojoule_per_mole
        ),
        "minimized_energy_kj_mol": minimized_energy.value_in_unit(
            openmm.unit.kilojoule_per_mole
        ),
        "minimization_rms_displacement_angstrom": float(
            np.sqrt(np.mean(np.sum(displacement**2, axis=1)))
        ),
        "minimization_max_displacement_angstrom": float(
            np.linalg.norm(displacement, axis=1).max()
        ),
        "openmm_nvt_smoke_test": nvt_result,
        "exported_files": exported_files,
    }
    destination = Path(
        "validation/results/"
        f"openff_cellodextrin_dp{degree_polymerization}{cache_infix}_probe.json"
    )
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
