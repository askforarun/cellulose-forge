"""Validate experimental OpenFF Sage support for a finite cellulose slab."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
from contextlib import chdir
from importlib import metadata
from pathlib import Path

import numpy as np
import openmm
from openff.interchange.drivers.gromacs import _process as process_gromacs_energies
from openff.interchange.drivers.gromacs import _run_gmx_energy
from openff.interchange.drivers.openmm import (
    _get_openmm_energies,
    _process,
)
from openff.toolkit import ForceField, Molecule, Topology
from openff.units import unit
from rdkit import Chem
from rdkit.Geometry import Point3D

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(REPOSITORY / "validation"))

from _provenance import provenance  # noqa: E402

from cellulose_forge.crystallography import MillerIndex, load_allomorph  # noqa: E402
from cellulose_forge.geometry import build_chain_parallel_slab  # noqa: E402
from cellulose_forge.glycam import generate_glycam_gromacs_reference  # noqa: E402
from cellulose_forge.glycam_mapping import map_crystal_to_finite_glycam  # noqa: E402

LAYERS = 2
CHAINS_PER_LAYER = 4
DEGREE_POLYMERIZATION = 6
CHAIN_COUNT = LAYERS * CHAINS_PER_LAYER
NVT_STEPS = 1000
MAXIMUM_BOND_CHANGE_ANGSTROM = 0.2
MINIMUM_XH_BOND_ANGSTROM = 0.8
MAXIMUM_XH_BOND_ANGSTROM = 1.2
MAXIMUM_HEAVY_BOND_DEVIATION_ANGSTROM = 0.12
CASES = (
    ("i_alpha", (0, 1, 0), "I-alpha (010)"),
    ("i_beta", (1, 1, 0), "I-beta (110)"),
    ("ii", (1, 1, 0), "II (110)"),
    ("iii_i", (1, 1, 0), "III-I (110)"),
)
REQUIRED_ENERGY_TERMS = frozenset(
    {"Bond", "Angle", "Torsion", "vdW", "vdW 1-4", "Electrostatics", "Electrostatics 1-4"}
)
ENERGY_RELATIVE_TOLERANCE = 1.0e-5


# Empirical mixed-precision acceptance criterion, not a general floating-point
# bound.  Mixed-precision GROMACS accumulates bonded energies in single
# precision, so its rounding error grows with the number of interactions
# summed (it also depends on term magnitudes, ordering and cancellation).
# validation/measure_openff_mixed_precision.py records the evidence: for the
# 36-chain I-alpha fibril, sequential single-precision summation errors are
# 2.25--3.22 eps32 sqrt(N).  The criterion uses 4 eps32 sqrt(N) and never falls
# below the fixed relative tolerance, which governs slab-sized systems.  If a
# case exceeds it, compare against a double-precision GROMACS build.
SINGLE_PRECISION_EPSILON = 2.0**-24
ACCUMULATION_ERROR_FACTOR = 4.0


def _tolerance(term: str, term_counts: dict[str, int] | None) -> float:
    count = (term_counts or {}).get(term, 0)
    return max(
        ENERGY_RELATIVE_TOLERANCE,
        ACCUMULATION_ERROR_FACTOR * SINGLE_PRECISION_EPSILON * float(np.sqrt(count)),
    )


def bonded_term_counts(system: openmm.System) -> dict[str, int]:
    """Count the interactions summed into each bonded energy category."""
    counts = {"Bond": 0, "Angle": 0, "Torsion": 0, "vdW 1-4": 0, "Electrostatics 1-4": 0}
    for force in system.getForces():
        if isinstance(force, openmm.HarmonicBondForce):
            counts["Bond"] += force.getNumBonds()
        elif isinstance(force, openmm.HarmonicAngleForce):
            counts["Angle"] += force.getNumAngles()
        elif isinstance(force, openmm.PeriodicTorsionForce):
            counts["Torsion"] += force.getNumTorsions()
        elif isinstance(force, openmm.CustomBondForce):
            # Interchange writes separate 1-4 Lennard-Jones and Coulomb forces
            # over the same pair list.
            counts["vdW 1-4"] = max(counts["vdW 1-4"], force.getNumBonds())
            counts["Electrostatics 1-4"] = counts["vdW 1-4"]
    return counts


def compare_energy_terms(
    openmm_values: dict[str, float],
    gromacs_values: dict[str, float],
    term_counts: dict[str, int] | None = None,
) -> dict[str, dict[str, float | bool]]:
    """Compare every required term with a documented relative tolerance."""
    missing_openmm = sorted(REQUIRED_ENERGY_TERMS - set(openmm_values))
    missing_gromacs = sorted(REQUIRED_ENERGY_TERMS - set(gromacs_values))
    if missing_openmm or missing_gromacs:
        raise RuntimeError(
            f"Missing required energy terms: OpenMM={missing_openmm}, GROMACS={missing_gromacs}"
        )
    comparison: dict[str, dict[str, float | bool]] = {}
    for term in sorted(REQUIRED_ENERGY_TERMS):
        relative = abs(openmm_values[term] - gromacs_values[term]) / max(
            abs(openmm_values[term]), abs(gromacs_values[term]), 1.0
        )
        tolerance = _tolerance(term, term_counts)
        comparison[term] = {
            "openmm_kj_mol": openmm_values[term],
            "gromacs_kj_mol": gromacs_values[term],
            "difference_kj_mol": openmm_values[term] - gromacs_values[term],
            "relative_difference": relative,
            "relative_tolerance": tolerance,
            "summed_interactions": (term_counts or {}).get(term, 0),
            "passed": relative <= tolerance,
        }
    return comparison

TIGHT_ENERGY_MDP = """integrator = md
nsteps = 0
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 0.9
fourierspacing = 0.05
ewald-rtol = 1e-6
pme-order = 8
vdwtype = Cut-off
vdw-modifier = Potential-switch
rvdw-switch = 0.8
rvdw = 0.9
DispCorr = no
constraints = none
periodic-molecules = no
nstenergy = 1
"""

MINIM_MDP = """integrator = steep
emtol = 1000
emstep = 0.01
nsteps = 5000
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 0.9
vdwtype = Cut-off
vdw-modifier = Potential-switch
rvdw-switch = 0.8
rvdw = 0.9
DispCorr = no
constraints = h-bonds
periodic-molecules = no
"""

NVT_MDP = """integrator = md
dt = 0.001
nsteps = 1000
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 0.9
vdwtype = Cut-off
vdw-modifier = Potential-switch
rvdw-switch = 0.8
rvdw = 0.9
DispCorr = no
constraints = h-bonds
constraint-algorithm = lincs
tcoupl = v-rescale
tc-grps = System
tau-t = 1.0
ref-t = 300
pcoupl = no
gen-vel = yes
gen-temp = 300
gen-seed = 20261112
periodic-molecules = no
nstlog = 100
nstenergy = 100
nstxout-compressed = 100
"""


def energies(report: object) -> dict[str, float]:
    """Convert an Interchange energy report to kJ/mol values."""
    return {
        name: float(value.m_as(unit.kilojoule_per_mole))
        for name, value in report.energies.items()  # type: ignore[attr-defined]
        if value is not None
    }


def disable_slab_dispersion_correction(system: openmm.System) -> None:
    """Disable homogeneous-fluid tail corrections for a vacuum slab."""
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.setUseDispersionCorrection(False)
        elif isinstance(force, openmm.CustomNonbondedForce):
            force.setUseLongRangeCorrection(False)


def openmm_slab_energies(interchange: object) -> dict[str, float]:
    """Evaluate detailed OpenMM energies without a slab-inappropriate tail term."""
    system = interchange.to_openmm_system(  # type: ignore[attr-defined]
        combine_nonbonded_forces=False,
        add_constrained_forces=True,
        ewald_tolerance=1.0e-6,
    )
    disable_slab_dispersion_correction(system)
    raw = _get_openmm_energies(
        system=system,
        box_vectors=interchange.box.to_openmm(),  # type: ignore[attr-defined]
        positions=interchange.positions.to_openmm(),  # type: ignore[attr-defined]
        round_positions=None,
        platform="Reference",
    )
    return energies(
        _process(
            raw,
            system=system,
            combine_nonbonded_forces=False,
            detailed=True,
        )
    )


def gromacs_slab_energies(interchange: object) -> dict[str, float]:
    """Evaluate detailed GROMACS energies with tight PME and no tail term."""
    with tempfile.TemporaryDirectory(prefix="cellulose_forge_openff_energy_") as tmp:
        directory = Path(tmp)
        prefix = str(directory / "structure")
        interchange.to_gromacs(prefix, decimal=8)  # type: ignore[attr-defined]
        mdp = directory / "point_energy.mdp"
        mdp.write_text(TIGHT_ENERGY_MDP, encoding="utf-8")
        normalized_mdp = re.sub(r"\s+", "", TIGHT_ENERGY_MDP.lower())
        if "dispcorr=no" not in normalized_mdp:
            raise RuntimeError("Vacuum-slab point energy must set DispCorr = no")
        with chdir(tmp):
            raw = _run_gmx_energy(
                top_file=f"{prefix}.top",
                gro_file=f"{prefix}.gro",
                mdp_file=mdp,
                maxwarn=0,
            )
        return energies(process_gromacs_energies(raw, detailed=True))


def input_geometry_checks(
    positions: np.ndarray,
    template: Molecule,
) -> dict[str, object]:
    """Check every chain's mapped covalent geometry before minimization."""
    atom_count = template.n_atoms
    heavy_distances: list[float] = []
    hydrogen_distances: list[float] = []
    for chain_index in range(CHAIN_COUNT):
        offset = chain_index * atom_count
        for bond in template.bonds:
            distance = float(
                np.linalg.norm(
                    positions[offset + bond.atom1_index]
                    - positions[offset + bond.atom2_index]
                )
            )
            elements = {
                template.atom(bond.atom1_index).symbol,
                template.atom(bond.atom2_index).symbol,
            }
            if "H" in elements:
                hydrogen_distances.append(distance)
            else:
                heavy_distances.append(distance)
    heavy_max_deviation = max(abs(value - 1.45) for value in heavy_distances)
    passed = bool(
        np.isfinite(hydrogen_distances).all()
        and np.isfinite(heavy_distances).all()
        and min(hydrogen_distances) >= MINIMUM_XH_BOND_ANGSTROM
        and max(hydrogen_distances) <= MAXIMUM_XH_BOND_ANGSTROM
        and heavy_max_deviation <= MAXIMUM_HEAVY_BOND_DEVIATION_ANGSTROM
    )
    return {
        "chains_checked": CHAIN_COUNT,
        "bonds_checked": CHAIN_COUNT * template.n_bonds,
        "hydrogen_bond_length_range_angstrom": [
            min(hydrogen_distances),
            max(hydrogen_distances),
        ],
        "heavy_bond_length_range_angstrom": [
            min(heavy_distances),
            max(heavy_distances),
        ],
        "maximum_heavy_bond_deviation_from_1_45_angstrom": heavy_max_deviation,
        "passed": passed,
    }


def stereochemistry_checks(
    positions: np.ndarray,
    template: Molecule,
) -> dict[str, object]:
    """Perceive and compare stereochemistry independently for every chain."""
    atom_count = template.n_atoms
    failures: list[int] = []

    def matches_template(chain_positions: np.ndarray) -> bool:
        rdmol = template.to_rdkit()
        rdmol.RemoveAllConformers()
        for atom in rdmol.GetAtoms():
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        conformer = Chem.Conformer(atom_count)
        for atom_index, xyz in enumerate(chain_positions):
            conformer.SetAtomPosition(atom_index, Point3D(*map(float, xyz)))
        rdmol.AddConformer(conformer, assignId=True)
        Chem.AssignStereochemistryFrom3D(
            rdmol,
            confId=0,
            replaceExistingTags=True,
        )
        perceived = Molecule.from_rdkit(
            rdmol,
            allow_undefined_stereo=False,
            hydrogens_are_explicit=True,
        )
        isomorphic, _ = Molecule.are_isomorphic(
            template,
            perceived,
            atom_stereochemistry_matching=True,
            bond_stereochemistry_matching=True,
            return_atom_map=True,
        )
        return bool(isomorphic)

    for chain_index in range(CHAIN_COUNT):
        chain_positions = positions[
            chain_index * atom_count : (chain_index + 1) * atom_count
        ]
        if not matches_template(chain_positions):
            failures.append(chain_index + 1)

    negative_control = positions[:atom_count].copy()
    negative_control[:, 0] = 2.0 * negative_control[:, 0].mean() - negative_control[:, 0]
    mirrored_chain_rejected = not matches_template(negative_control)
    return {
        "chains_checked": CHAIN_COUNT,
        "failed_chain_numbers": failures,
        "negative_control": "mirror first chain through its x centroid",
        "mirrored_chain_rejected": mirrored_chain_rejected,
        "passed": not failures and mirrored_chain_rejected,
    }


def run_gromacs_smoke_test(interchange: object) -> dict[str, object]:
    """Run GROMACS preprocessing, minimization, and 1 ps NVT."""
    with tempfile.TemporaryDirectory(prefix="cellulose_forge_openff_gromacs_") as tmp:
        directory = Path(tmp)
        version_result = subprocess.run(
            ["gmx", "--version"],
            text=True,
            capture_output=True,
            check=True,
        )
        version_match = re.search(
            r"^GROMACS version:\s+(.+)$",
            version_result.stdout + version_result.stderr,
            flags=re.MULTILINE,
        )
        version = version_match.group(1) if version_match else "unknown"
        prefix = str(directory / "structure")
        interchange.to_gromacs(prefix)  # type: ignore[attr-defined]
        minim_mdp = directory / "minim.mdp"
        nvt_mdp = directory / "nvt.mdp"
        minim_mdp.write_text(MINIM_MDP, encoding="utf-8")
        nvt_mdp.write_text(NVT_MDP, encoding="utf-8")

        commands = [
            [
                "gmx",
                "grompp",
                "-f",
                str(minim_mdp),
                "-c",
                f"{prefix}.gro",
                "-p",
                f"{prefix}.top",
                "-o",
                str(directory / "minim.tpr"),
                "-maxwarn",
                "0",
            ],
            [
                "gmx",
                "mdrun",
                "-deffnm",
                str(directory / "minim"),
                "-ntmpi",
                "1",
            ],
            [
                "gmx",
                "grompp",
                "-f",
                str(nvt_mdp),
                "-c",
                str(directory / "minim.gro"),
                "-p",
                f"{prefix}.top",
                "-o",
                str(directory / "nvt.tpr"),
                "-maxwarn",
                "0",
            ],
            [
                "gmx",
                "mdrun",
                "-deffnm",
                str(directory / "nvt"),
                "-ntmpi",
                "1",
            ],
        ]
        outputs = []
        for command in commands:
            completed = subprocess.run(
                command,
                cwd=directory,
                text=True,
                capture_output=True,
                check=True,
            )
            outputs.append(completed.stdout + completed.stderr)
        combined = "\n".join(outputs)
        forbidden = {
            "LINCS WARNING": r"LINCS WARNING",
            "nan": r"\bnan\b",
            "non-finite": r"non-finite",
        }
        failed_markers = [
            marker
            for marker, pattern in forbidden.items()
            if re.search(pattern, combined, flags=re.IGNORECASE)
        ]
        nvt_log = (directory / "nvt.log").read_text(encoding="utf-8")
        reached_final_step = "Statistics over 1001 steps" in nvt_log
        return {
            "gromacs_version": version,
            "grompp_maxwarn": 0,
            "minimization_completed": (directory / "minim.gro").exists(),
            "nvt_duration_ps": 1.0,
            "nvt_reached_final_step": reached_final_step,
            "failure_markers": failed_markers,
            "passed": bool(
                (directory / "minim.gro").exists()
                and reached_final_step
                and not failed_markers
            ),
        }


def system_bond_lengths(positions: np.ndarray, molecule: Molecule) -> np.ndarray:
    """Measure every bond in every identical chain."""
    atom_count = molecule.n_atoms
    return np.asarray(
        [
            np.linalg.norm(
                positions[offset + bond.atom1_index]
                - positions[offset + bond.atom2_index]
            )
            for chain_index in range(CHAIN_COUNT)
            for offset in [chain_index * atom_count]
            for bond in molecule.bonds
        ]
    )


def validate_case(
    allomorph: str,
    face: tuple[int, int, int],
    label: str,
    template: Molecule,
    charge_method: str,
) -> dict[str, object]:
    """Validate one finite reference slab in OpenMM and GROMACS."""
    crystal = build_chain_parallel_slab(
        load_allomorph(allomorph),
        MillerIndex(*face),
        layers=LAYERS,
        chains_per_layer=CHAINS_PER_LAYER,
        degree_polymerization=DEGREE_POLYMERIZATION,
        vacuum=20.0,
        padding_z=12.0,
        chain_boundary="finite",
    )
    with tempfile.TemporaryDirectory(prefix="cellulose_forge_openff_slab_") as temporary:
        files, _ = generate_glycam_gromacs_reference(
            degree_polymerization=DEGREE_POLYMERIZATION,
            chain_boundary="finite",
            directory=Path(temporary),
            parmed_python=sys.executable,
        )
        mapped = map_crystal_to_finite_glycam(crystal, files.pdb)

    expected_atoms = CHAIN_COUNT * template.n_atoms
    if len(mapped.atom_names) != expected_atoms:
        raise RuntimeError(
            f"Expected {expected_atoms} mapped atoms, found {len(mapped.atom_names)}"
        )

    topology = Topology.from_molecules([template] * CHAIN_COUNT)
    topology.box_vectors = crystal.basis_vectors * unit.angstrom
    force_field = ForceField("openff-2.3.0.offxml")
    interchange = force_field.create_interchange(
        topology,
        charge_from_molecules=[template],
    )
    interchange.positions = mapped.cartesian_coordinates * unit.angstrom
    geometry_checks = input_geometry_checks(mapped.cartesian_coordinates, template)
    if not geometry_checks["passed"]:
        raise RuntimeError(f"Mapped input geometry failed validation: {geometry_checks}")
    stereo_checks = stereochemistry_checks(mapped.cartesian_coordinates, template)
    if not stereo_checks["passed"]:
        raise RuntimeError(f"Mapped stereochemistry failed validation: {stereo_checks}")

    openmm_values = openmm_slab_energies(interchange)
    gromacs_values = gromacs_slab_energies(interchange)
    comparison = compare_energy_terms(
        openmm_values,
        gromacs_values,
        bonded_term_counts(
            interchange.to_openmm_system(combine_nonbonded_forces=False)  # type: ignore[attr-defined]
        ),
    )
    failed_energy_terms = [
        term for term, values in comparison.items() if not values["passed"]
    ]

    system = interchange.to_openmm_system(ewald_tolerance=1.0e-6)
    disable_slab_dispersion_correction(system)
    integrator = openmm.LangevinMiddleIntegrator(
        300.0 * openmm.unit.kelvin,
        1.0 / openmm.unit.picosecond,
        0.001 * openmm.unit.picoseconds,
    )
    # Seed the thermostat noise so stored NVT results are reproducible.
    integrator.setRandomNumberSeed(20261112)
    platform = openmm.Platform.getPlatformByName("CPU")
    platform_properties = {"Threads": "1"}
    context = openmm.Context(system, integrator, platform, platform_properties)
    context.setPeriodicBoxVectors(*interchange.box.to_openmm())
    context.setPositions(interchange.positions.to_openmm())
    initial_state = context.getState(getEnergy=True, getPositions=True)
    openmm.LocalEnergyMinimizer.minimize(context, tolerance=10.0, maxIterations=2000)
    minimized_state = context.getState(getEnergy=True, getPositions=True)
    minimized_positions = minimized_state.getPositions(asNumpy=True).value_in_unit(
        openmm.unit.angstrom
    )
    minimized_bonds = system_bond_lengths(minimized_positions, template)
    context.setVelocitiesToTemperature(300.0 * openmm.unit.kelvin, 20261112)
    integrator.step(NVT_STEPS)
    final_state = context.getState(getEnergy=True, getPositions=True)
    final_positions = final_state.getPositions(asNumpy=True).value_in_unit(
        openmm.unit.angstrom
    )
    maximum_bond_change = float(
        np.max(np.abs(system_bond_lengths(final_positions, template) - minimized_bonds))
    )

    return {
        "system": f"CelluloseForge {label} finite DP6 slab",
        "allomorph": allomorph,
        "face": list(face),
        "layers": LAYERS,
        "chains_per_layer": CHAINS_PER_LAYER,
        "chain_count": CHAIN_COUNT,
        "atom_count": expected_atoms,
        "openff_force_field": "openff-2.3.0.offxml",
        "charge_method": f"cached {charge_method} DP6 template",
        "input_geometry_checks": geometry_checks,
        "stereochemistry_checks": stereo_checks,
        "cross_engine_single_point": {
            "dispersion_correction": "disabled for vacuum slab in both engines",
            "pme_accuracy": {
                "openmm_ewald_error_tolerance": 1.0e-6,
                "gromacs_ewald_rtol": 1.0e-6,
                "gromacs_fourierspacing_nm": 0.05,
                "gromacs_pme_order": 8,
            },
            "comparison": comparison,
            "passed": not failed_energy_terms,
            "failed_terms": failed_energy_terms,
        },
        "initial_potential_energy_kj_mol": initial_state.getPotentialEnergy().value_in_unit(
            openmm.unit.kilojoule_per_mole
        ),
        "minimized_potential_energy_kj_mol": (
            minimized_state.getPotentialEnergy().value_in_unit(
                openmm.unit.kilojoule_per_mole
            )
        ),
        "openmm_nvt_smoke_test": {
            "duration_ps": 1.0,
            "platform": platform.getName(),
            "platform_properties": platform_properties,
            "maximum_bond_length_change_angstrom": maximum_bond_change,
            "final_potential_energy_kj_mol": final_state.getPotentialEnergy().value_in_unit(
                openmm.unit.kilojoule_per_mole
            ),
            "passed": bool(
                np.isfinite(final_positions).all()
                and maximum_bond_change < MAXIMUM_BOND_CHANGE_ANGSTROM
            ),
        },
        "gromacs_smoke_test": run_gromacs_smoke_test(interchange),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--charge-method", choices=("am1bcc", "nagl"), default="nagl")
    args = parser.parse_args()
    charge_method = args.charge_method
    cache_infix = "" if charge_method == "am1bcc" else f"_{charge_method}"
    cache = Path(
        f"validation/results/openff_cellodextrin_dp6{cache_infix}_charged.json"
    )
    if not cache.exists():
        raise FileNotFoundError(
            "Run probe_openff_cellobiose.py --dp 6 "
            f"--charge-method {charge_method} first"
        )
    template = Molecule.from_json(cache.read_text(encoding="utf-8"))
    cases = {
        allomorph: validate_case(
            allomorph,
            face,
            label,
            template,
            charge_method,
        )
        for allomorph, face, label in CASES
    }
    result = {
        "matrix": "finite DP6 cellulose reference slabs",
        "provenance": provenance(
            ("openff-toolkit", "openff-interchange", "openff-nagl", "openff-nagl-models", "openmm")
        ),
        "charge_method": charge_method,
        "case_count": len(cases),
        "software_versions": {
            "openff-toolkit": metadata.version("openff-toolkit"),
            "openff-interchange": metadata.version("openff-interchange"),
            "openmm": metadata.version("openmm"),
            "gromacs": cases["i_beta"]["gromacs_smoke_test"][  # type: ignore[index]
                "gromacs_version"
            ],
        },
        "cases": cases,
        "passed": all(
            case["input_geometry_checks"]["passed"]  # type: ignore[index]
            and case["stereochemistry_checks"]["passed"]  # type: ignore[index]
            and case["cross_engine_single_point"]["passed"]  # type: ignore[index]
            and case["openmm_nvt_smoke_test"]["passed"]  # type: ignore[index]
            and case["gromacs_smoke_test"]["passed"]  # type: ignore[index]
            for case in cases.values()
        ),
    }
    destination = Path(
        f"validation/results/openff_finite_slab{cache_infix}_matrix.json"
    )
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
