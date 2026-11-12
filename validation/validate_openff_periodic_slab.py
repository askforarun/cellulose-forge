"""Validate experimental periodic-chain OpenFF cellulose slabs."""

from __future__ import annotations

import hashlib
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
from openff.interchange.drivers.openmm import _get_openmm_energies, _process
from openff.toolkit import ForceField, Molecule, Topology
from openff.units import unit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "validation"))

from validate_openff_finite_slab import (  # noqa: E402
    MAXIMUM_BOND_CHANGE_ANGSTROM,
    MAXIMUM_HEAVY_BOND_DEVIATION_ANGSTROM,
    MAXIMUM_XH_BOND_ANGSTROM,
    MINIM_MDP,
    MINIMUM_XH_BOND_ANGSTROM,
    NVT_MDP,
    TIGHT_ENERGY_MDP,
    compare_energy_terms,
    disable_slab_dispersion_correction,
    energies,
)

from cellulose_forge.crystallography import MillerIndex, load_allomorph  # noqa: E402
from cellulose_forge.geometry import build_chain_parallel_slab  # noqa: E402
from cellulose_forge.glycam import generate_glycam_gromacs_reference  # noqa: E402
from cellulose_forge.glycam_mapping import map_crystal_to_periodic_glycam  # noqa: E402

CASES = (
    ("i_alpha", (0, 1, 0), "I-alpha (010)"),
    ("i_beta", (1, 1, 0), "I-beta (110)"),
    ("ii", (1, 1, 0), "II (110)"),
    ("iii_i", (1, 1, 0), "III-I (110)"),
)
LAYERS = 2
CHAINS_PER_LAYER = 4
CHAIN_COUNT = LAYERS * CHAINS_PER_LAYER
DP = 6
FIRST_RESIDUE = 2
LAST_RESIDUE = 7
REFERENCE_INTERNAL_RESIDUE = 3
MINIMUM_GLYCOSIDIC_BOND_ANGSTROM = 1.3
MAXIMUM_GLYCOSIDIC_BOND_ANGSTROM = 1.6


def enable_periodic_molecular_interactions(system: openmm.System) -> None:
    """Use minimum-image distances for every bonded term and exception."""
    for force in system.getForces():
        setter = getattr(force, "setUsesPeriodicBoundaryConditions", None)
        if setter is not None:
            setter(True)
        if isinstance(force, openmm.NonbondedForce):
            force.setExceptionsUsePeriodicBoundaryConditions(True)


def minimum_image(displacement: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Return a triclinic minimum-image displacement."""
    fractional = np.linalg.solve(box.T, displacement)
    fractional -= np.round(fractional)
    return fractional @ box


def periodic_bond_checks(
    positions: np.ndarray, box: np.ndarray, template: Molecule
) -> dict[str, object]:
    """Check every closed-loop glycosidic bond, including box crossings."""
    first_c1 = next(
        atom.molecule_atom_index
        for atom in template.atoms
        if atom.name == "C1" and int(atom.metadata["residue_number"]) == FIRST_RESIDUE
    )
    last_o4 = next(
        atom.molecule_atom_index
        for atom in template.atoms
        if atom.name == "O4" and int(atom.metadata["residue_number"]) == LAST_RESIDUE
    )
    lengths = []
    for chain in range(CHAIN_COUNT):
        offset = chain * template.n_atoms
        delta = positions[offset + first_c1] - positions[offset + last_o4]
        lengths.append(float(np.linalg.norm(minimum_image(delta, box))))
    return {
        "boundary_bonds_checked": CHAIN_COUNT,
        "length_range_angstrom": [min(lengths), max(lengths)],
        "passed": bool(np.isfinite(lengths).all())
        and min(lengths) >= MINIMUM_GLYCOSIDIC_BOND_ANGSTROM
        and max(lengths) <= MAXIMUM_GLYCOSIDIC_BOND_ANGSTROM,
    }


def unwrap_chain_for_stereochemistry(
    positions: np.ndarray, box: np.ndarray, template: Molecule
) -> np.ndarray:
    """Unwrap a chain along all bonds except the box-closing linkage."""
    first_c1 = next(
        atom.molecule_atom_index
        for atom in template.atoms
        if atom.name == "C1" and int(atom.metadata["residue_number"]) == FIRST_RESIDUE
    )
    last_o4 = next(
        atom.molecule_atom_index
        for atom in template.atoms
        if atom.name == "O4" and int(atom.metadata["residue_number"]) == LAST_RESIDUE
    )
    closing = frozenset((first_c1, last_o4))
    adjacency: dict[int, list[int]] = {index: [] for index in range(template.n_atoms)}
    for bond in template.bonds:
        pair = frozenset((bond.atom1_index, bond.atom2_index))
        if pair == closing:
            continue
        adjacency[bond.atom1_index].append(bond.atom2_index)
        adjacency[bond.atom2_index].append(bond.atom1_index)
    unwrapped = np.empty_like(positions)
    unwrapped[0] = positions[0]
    visited = {0}
    pending = [0]
    while pending:
        source = pending.pop()
        for target in adjacency[source]:
            if target in visited:
                continue
            delta = minimum_image(positions[target] - positions[source], box)
            unwrapped[target] = unwrapped[source] + delta
            visited.add(target)
            pending.append(target)
    if len(visited) != template.n_atoms:
        raise RuntimeError("Periodic chain did not unwrap into one connected graph")
    return unwrapped


STEREOCENTRE_NAMES = ("C1", "C2", "C3", "C4", "C5")


def centre_parity(
    positions: np.ndarray,
    atom: object,
    box: np.ndarray | None = None,
) -> int:
    """Return signed local parity using a deterministic neighbour order."""
    neighbours = sorted(atom.bonded_atoms, key=lambda item: (item.symbol, item.name))[:3]
    centre = positions[atom.molecule_atom_index]
    vectors = [positions[item.molecule_atom_index] - centre for item in neighbours]
    if box is not None:
        vectors = [minimum_image(vector, box) for vector in vectors]
    value = float(np.dot(vectors[0], np.cross(vectors[1], vectors[2])))
    if not np.isfinite(value) or value == 0.0:
        raise RuntimeError("Undefined stereocentre parity")
    return int(np.sign(value))


def independent_stereo_reference() -> tuple[dict[str, int], dict[str, str]]:
    """Read parity and absolute labels from finite PubChem internal residue 3."""
    finite = Molecule.from_json(
        (ROOT / "validation/results/openff_cellodextrin_dp6_nagl_charged.json").read_text(
            encoding="utf-8"
        )
    )
    positions = finite.conformers[0].m_as(unit.angstrom)
    atoms = [
        atom
        for atom in finite.atoms
        if atom.name in STEREOCENTRE_NAMES
        and int(atom.metadata["residue_number"]) == REFERENCE_INTERNAL_RESIDUE
    ]
    return (
        {atom.name: centre_parity(positions, atom) for atom in atoms},
        {atom.name: str(atom.stereochemistry) for atom in atoms},
    )


def independent_periodic_stereochemistry_checks(
    positions: np.ndarray,
    box: np.ndarray,
    template: Molecule,
    reference_parity: dict[str, int],
    reference_labels: dict[str, str],
) -> dict[str, object]:
    """Check every cyclic residue against independent finite-chain chemistry."""
    label_failures = [
        (int(atom.metadata["residue_number"]), atom.name, atom.stereochemistry)
        for atom in template.atoms
        if atom.name in STEREOCENTRE_NAMES and atom.stereochemistry != reference_labels[atom.name]
    ]
    parity_failures: list[tuple[int, int, str]] = []
    for chain in range(CHAIN_COUNT):
        chain_positions = positions[chain * template.n_atoms : (chain + 1) * template.n_atoms]
        for atom in template.atoms:
            if atom.name in STEREOCENTRE_NAMES:
                if centre_parity(chain_positions, atom, box) != reference_parity[atom.name]:
                    parity_failures.append(
                        (chain + 1, int(atom.metadata["residue_number"]), atom.name)
                    )

    controls: dict[str, bool] = {}
    first_chain = positions[: template.n_atoms]
    for name in ("C1", "C3"):
        target = next(
            atom
            for atom in template.atoms
            if atom.name == name
            and int(atom.metadata["residue_number"]) == REFERENCE_INTERNAL_RESIDUE
        )
        altered = first_chain.copy()
        neighbours = sorted(target.bonded_atoms, key=lambda item: (item.symbol, item.name))[:3]
        left, right = neighbours[0].molecule_atom_index, neighbours[1].molecule_atom_index
        altered[[left, right]] = altered[[right, left]]
        controls[f"{name}_local_inversion_rejected"] = (
            centre_parity(altered, target, box) != reference_parity[name]
        )
    mirrored = first_chain.copy()
    mirrored[:, 0] = 2.0 * mirrored[:, 0].mean() - mirrored[:, 0]
    controls["whole_chain_mirror_rejected"] = all(
        centre_parity(mirrored, atom, box) != reference_parity[atom.name]
        for atom in template.atoms
        if atom.name in STEREOCENTRE_NAMES
    )
    return {
        "reference": "PubChem-derived finite DP6 internal residue 3",
        "centres_checked": CHAIN_COUNT * DP * len(STEREOCENTRE_NAMES),
        "template_label_failures": label_failures,
        "coordinate_parity_failures": parity_failures,
        "negative_controls": controls,
        "passed": not label_failures and not parity_failures and all(controls.values()),
    }


def periodic_geometry_checks(
    positions: np.ndarray, box: np.ndarray, template: Molecule
) -> dict[str, object]:
    """Check all covalent distances using minimum-image displacements."""
    heavy: list[float] = []
    hydrogen: list[float] = []
    for chain in range(CHAIN_COUNT):
        offset = chain * template.n_atoms
        for bond in template.bonds:
            delta = positions[offset + bond.atom2_index] - positions[offset + bond.atom1_index]
            distance = float(np.linalg.norm(minimum_image(delta, box)))
            symbols = {
                template.atom(bond.atom1_index).symbol,
                template.atom(bond.atom2_index).symbol,
            }
            (hydrogen if "H" in symbols else heavy).append(distance)
    maximum_deviation = max(abs(value - 1.45) for value in heavy)
    return {
        "chains_checked": CHAIN_COUNT,
        "bonds_checked": CHAIN_COUNT * template.n_bonds,
        "hydrogen_bond_length_range_angstrom": [min(hydrogen), max(hydrogen)],
        "heavy_bond_length_range_angstrom": [min(heavy), max(heavy)],
        "maximum_heavy_bond_deviation_from_1_45_angstrom": maximum_deviation,
        "passed": bool(np.isfinite(hydrogen).all() and np.isfinite(heavy).all())
        and min(hydrogen) >= MINIMUM_XH_BOND_ANGSTROM
        and max(hydrogen) <= MAXIMUM_XH_BOND_ANGSTROM
        and maximum_deviation <= MAXIMUM_HEAVY_BOND_DEVIATION_ANGSTROM,
    }


def periodic_all_bond_lengths(
    positions: np.ndarray, box: np.ndarray, template: Molecule
) -> np.ndarray:
    """Measure every molecular bond by minimum image."""
    return np.asarray(
        [
            np.linalg.norm(
                minimum_image(
                    positions[offset + bond.atom2_index] - positions[offset + bond.atom1_index],
                    box,
                )
            )
            for chain in range(CHAIN_COUNT)
            for offset in [chain * template.n_atoms]
            for bond in template.bonds
        ]
    )


def openmm_values(interchange: object) -> tuple[dict[str, float], openmm.System]:
    system = interchange.to_openmm_system(  # type: ignore[attr-defined]
        combine_nonbonded_forces=False,
        add_constrained_forces=True,
        ewald_tolerance=1.0e-6,
    )
    enable_periodic_molecular_interactions(system)
    disable_slab_dispersion_correction(system)
    raw = _get_openmm_energies(
        system=system,
        box_vectors=interchange.box.to_openmm(),  # type: ignore[attr-defined]
        positions=interchange.positions.to_openmm(),  # type: ignore[attr-defined]
        round_positions=None,
        platform="Reference",
    )
    report = _process(raw, system=system, combine_nonbonded_forces=False, detailed=True)
    return energies(report), system


def gromacs_values(interchange: object) -> dict[str, float]:
    with tempfile.TemporaryDirectory(prefix="cellulose_forge_openff_periodic_") as tmp:
        prefix = str(Path(tmp) / "structure")
        interchange.to_gromacs(prefix, decimal=8)  # type: ignore[attr-defined]
        mdp = Path(tmp) / "point.mdp"
        mdp.write_text(
            TIGHT_ENERGY_MDP.replace("periodic-molecules = no", "periodic-molecules = yes"),
            encoding="utf-8",
        )
        with chdir(tmp):
            raw = _run_gmx_energy(
                top_file=f"{prefix}.top", gro_file=f"{prefix}.gro", mdp_file=mdp, maxwarn=0
            )
        return energies(process_gromacs_energies(raw, detailed=True))


def run_gromacs_periodic_smoke_test(interchange: object) -> dict[str, object]:
    """Run isolated GROMACS minimization and 1 ps NVT for a periodic slab."""
    with tempfile.TemporaryDirectory(prefix="cellulose_forge_openff_periodic_gmx_") as tmp:
        directory = Path(tmp)
        prefix = str(directory / "structure")
        interchange.to_gromacs(prefix, decimal=8)  # type: ignore[attr-defined]
        minim = MINIM_MDP.replace("periodic-molecules = no", "periodic-molecules = yes")
        nvt = NVT_MDP.replace("periodic-molecules = no", "periodic-molecules = yes")
        (directory / "minim.mdp").write_text(minim, encoding="utf-8")
        (directory / "nvt.mdp").write_text(nvt, encoding="utf-8")
        commands = (
            [
                "gmx",
                "grompp",
                "-f",
                "minim.mdp",
                "-c",
                "structure.gro",
                "-p",
                "structure.top",
                "-o",
                "minim.tpr",
                "-maxwarn",
                "0",
            ],
            ["gmx", "mdrun", "-deffnm", "minim", "-ntmpi", "1"],
            [
                "gmx",
                "grompp",
                "-f",
                "nvt.mdp",
                "-c",
                "minim.gro",
                "-p",
                "structure.top",
                "-o",
                "nvt.tpr",
                "-maxwarn",
                "0",
            ],
            ["gmx", "mdrun", "-deffnm", "nvt", "-ntmpi", "1"],
        )
        logs = []
        for command in commands:
            completed = subprocess.run(
                command, cwd=directory, text=True, capture_output=True, check=True
            )
            logs.append(completed.stdout + completed.stderr)
        combined = "\n".join(logs) + (directory / "nvt.log").read_text(encoding="utf-8")
        failures = [
            marker
            for marker, pattern in {
                "LINCS": r"LINCS WARNING",
                "NaN": r"\bnan\b",
                "non-finite": r"non-finite",
            }.items()
            if re.search(pattern, combined, flags=re.IGNORECASE)
        ]
        reached = "Statistics over 1001 steps" in combined
        return {
            "grompp_maxwarn": 0,
            "minimization_completed": (directory / "minim.gro").exists(),
            "nvt_duration_ps": 1.0,
            "nvt_reached_final_step": reached,
            "failure_markers": failures,
            "passed": (directory / "minim.gro").exists() and reached and not failures,
        }


def validate_case(  # noqa: PLR0913, PLR0917
    allomorph: str,
    face: tuple[int, int, int],
    label: str,
    template: Molecule,
    reference_parity: dict[str, int],
    reference_labels: dict[str, str],
) -> dict[str, object]:
    crystal = build_chain_parallel_slab(
        load_allomorph(allomorph),
        MillerIndex(*face),
        layers=LAYERS,
        chains_per_layer=CHAINS_PER_LAYER,
        degree_polymerization=DP,
        vacuum=20.0,
        chain_boundary="periodic",
    )
    with tempfile.TemporaryDirectory(prefix="cellulose_forge_openff_periodic_map_") as tmp:
        files, _ = generate_glycam_gromacs_reference(
            degree_polymerization=DP,
            chain_boundary="periodic",
            directory=tmp,
            parmed_python=sys.executable,
        )
        mapped = map_crystal_to_periodic_glycam(crystal, files.pdb)
    if mapped.atom_count != CHAIN_COUNT * template.n_atoms:
        raise RuntimeError("Periodic mapped atom count does not match OpenFF graph")
    expected_names = tuple(atom.name for atom in template.atoms) * CHAIN_COUNT
    if tuple(mapped.atom_names) != expected_names:
        first = next(
            index
            for index, (observed, expected) in enumerate(
                zip(mapped.atom_names, expected_names, strict=True)
            )
            if observed != expected
        )
        raise RuntimeError(
            "Periodic mapped atom ordering does not match OpenFF graph: "
            f"index {first}, mapped={mapped.atom_names[first]}, "
            f"OpenFF={expected_names[first]}"
        )

    topology = Topology.from_molecules([template] * CHAIN_COUNT)
    topology.box_vectors = crystal.basis_vectors * unit.angstrom
    interchange = ForceField("openff-2.3.0.offxml").create_interchange(
        topology, charge_from_molecules=[template]
    )
    interchange.positions = mapped.cartesian_coordinates * unit.angstrom
    geometry = periodic_geometry_checks(
        mapped.cartesian_coordinates, crystal.basis_vectors, template
    )
    stereo = independent_periodic_stereochemistry_checks(
        mapped.cartesian_coordinates,
        crystal.basis_vectors,
        template,
        reference_parity,
        reference_labels,
    )
    boundary = periodic_bond_checks(mapped.cartesian_coordinates, crystal.basis_vectors, template)
    omm, system = openmm_values(interchange)
    gmx = gromacs_values(interchange)
    gromacs_smoke = run_gromacs_periodic_smoke_test(interchange)
    comparison = compare_energy_terms(omm, gmx)
    platform = openmm.Platform.getPlatformByName("CPU")
    platform_properties = {"Threads": "1"}
    minimization_integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, minimization_integrator, platform, platform_properties)
    context.setPeriodicBoxVectors(*interchange.box.to_openmm())
    context.setPositions(interchange.positions.to_openmm())
    initial = (
        context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(openmm.unit.kilojoule_per_mole)
    )
    openmm.LocalEnergyMinimizer.minimize(context, tolerance=10.0, maxIterations=2000)
    minimized = context.getState(getEnergy=True, getPositions=True)
    minimized_positions = minimized.getPositions(asNumpy=True).value_in_unit(openmm.unit.angstrom)
    minimized_boundary = periodic_bond_checks(minimized_positions, crystal.basis_vectors, template)
    minimized_all_bonds = periodic_geometry_checks(
        minimized_positions, crystal.basis_vectors, template
    )
    del context
    integrator = openmm.LangevinMiddleIntegrator(
        300.0 * openmm.unit.kelvin,
        1.0 / openmm.unit.picosecond,
        0.001 * openmm.unit.picoseconds,
    )
    context = openmm.Context(system, integrator, platform, platform_properties)
    context.setPeriodicBoxVectors(*interchange.box.to_openmm())
    context.setPositions(minimized.getPositions())
    context.setVelocitiesToTemperature(300.0 * openmm.unit.kelvin, 20261112)
    integrator.step(1000)
    final = context.getState(getEnergy=True, getPositions=True)
    final_positions = final.getPositions(asNumpy=True).value_in_unit(openmm.unit.angstrom)
    final_boundary = periodic_bond_checks(final_positions, crystal.basis_vectors, template)
    final_all_bonds = periodic_geometry_checks(final_positions, crystal.basis_vectors, template)
    # The crystallographic absolute-range threshold is an input/minimized-
    # geometry check, not a thermal-frame stability criterion.  Keep the
    # final ranges as diagnostics and use the explicit minimized-to-final
    # all-bond change below to decide whether the NVT trajectory passed.
    final_all_bonds.pop("passed")
    minimized_lengths = periodic_all_bond_lengths(
        minimized_positions, crystal.basis_vectors, template
    )
    final_lengths = periodic_all_bond_lengths(final_positions, crystal.basis_vectors, template)
    maximum_bond_change = float(np.max(np.abs(final_lengths - minimized_lengths)))
    all_bond_stability = {
        "bonds_checked": len(final_lengths),
        "maximum_change_angstrom": maximum_bond_change,
        "all_finite": bool(
            np.isfinite(minimized_lengths).all() and np.isfinite(final_lengths).all()
        ),
        "passed": bool(
            np.isfinite(maximum_bond_change) and maximum_bond_change < MAXIMUM_BOND_CHANGE_ANGSTROM
        ),
    }
    return {
        "system": label,
        "atom_count": mapped.atom_count,
        "geometry": geometry,
        "stereochemistry": stereo,
        "periodic_boundary_bonds": boundary,
        "energy_comparison": comparison,
        "initial_energy_kj_mol": initial,
        "minimized_energy_kj_mol": minimized.getPotentialEnergy().value_in_unit(
            openmm.unit.kilojoule_per_mole
        ),
        "openmm_nvt": {
            "duration_ps": 1.0,
            "platform": platform.getName(),
            "platform_properties": platform_properties,
            "minimized_boundary_bonds": minimized_boundary,
            "minimized_all_bonds": minimized_all_bonds,
            "final_boundary_bonds": final_boundary,
            "final_all_bond_ranges": final_all_bonds,
            "all_bond_stability": all_bond_stability,
            "final_energy_kj_mol": final.getPotentialEnergy().value_in_unit(
                openmm.unit.kilojoule_per_mole
            ),
            "passed": minimized_boundary["passed"]
            and final_boundary["passed"]
            and minimized_all_bonds["passed"]
            and all_bond_stability["passed"],
        },
        "gromacs_smoke_test": gromacs_smoke,
        "passed": bool(
            geometry["passed"]
            and stereo["passed"]
            and boundary["passed"]
            and all(value["passed"] for value in comparison.values())
            and minimized_boundary["passed"]
            and final_boundary["passed"]
            and minimized_all_bonds["passed"]
            and all_bond_stability["passed"]
            and gromacs_smoke["passed"]
        ),
    }


def main() -> None:
    template = Molecule.from_json(
        (ROOT / "validation/results/openff_cellodextrin_dp6_periodic_nagl_charged.json").read_text(
            encoding="utf-8"
        )
    )
    reference_parity, reference_labels = independent_stereo_reference()
    cases = {
        name: validate_case(name, face, label, template, reference_parity, reference_labels)
        for name, face, label in CASES
    }
    template_path = ROOT / "validation/results/openff_cellodextrin_dp6_periodic_nagl_charged.json"
    result = {
        "matrix": "periodic DP6 cellulose reference slabs",
        "charge_method": "NAGL openff-gnn-am1bcc-1.0.0.pt",
        "template_sha256": hashlib.sha256(template_path.read_bytes()).hexdigest(),
        "software_versions": {
            "openff-toolkit": metadata.version("openff-toolkit"),
            "openff-interchange": metadata.version("openff-interchange"),
            "openff-nagl": metadata.version("openff-nagl"),
            "openff-nagl-models": metadata.version("openff-nagl-models"),
            "openmm": metadata.version("openmm"),
            "gromacs": "2024.2-Homebrew",
        },
        "pme_settings": {
            "openmm_ewald_tolerance": 1e-6,
            "gromacs_ewald_rtol": 1e-6,
            "gromacs_pme_order": 8,
            "gromacs_fourierspacing_nm": 0.05,
        },
        "cases": cases,
        "passed": all(case["passed"] for case in cases.values()),
    }
    output = ROOT / "validation/results/openff_periodic_slab_nagl_matrix.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
