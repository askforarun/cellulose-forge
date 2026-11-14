"""Compare GLYCAM chitin energies in GROMACS, OpenMM, and LAMMPS.

The matrix covers finite and covalently periodic (100) slabs and periodic bulk
crystals of alpha- and beta-chitin. Explicit amide 1-4 terms are evaluated
analytically so their LAMMPS tabulated pseudo-bond contribution can be separated
from harmonic bonds and compared with the corresponding GROMACS categories.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
from pathlib import Path
from typing import cast

import numpy as np
from _provenance import provenance
from numpy.typing import NDArray
from validate_chitin_glycam import _build_arguments, _run, _run_with_input

from cellulose_forge.cli import main
from cellulose_forge.glycam import find_parmed_python
from cellulose_forge.gromacs_pairs import prepare_explicit_pairs

_ROOT = Path(__file__).resolve().parents[1]
_KJ_PER_KCAL = 4.184
_LAMMPS_COULOMB = 332.06371
# PME (GROMACS/OpenMM) and PPPM (LAMMPS) use different reciprocal-space
# approximations. The established cross-engine matrix permits 2e-4 for the
# combined nonbonded energy while bonded and explicit 1-4 categories remain
# near output precision.
_MAXIMUM_RELATIVE_DIFFERENCE = 2.0e-4
_MAXIMUM_BONDED_RELATIVE_DIFFERENCE = 1.0e-5
_ORTHOGONAL_GRO_BOX_FIELDS = 3
_TRICLINIC_GRO_BOX_FIELDS = 9


def _gromacs_energies(
    gmx: str, directory: Path, *, periodic: bool, dispersion: bool
) -> dict[str, float]:
    mdp = directory / "cross_engine.mdp"
    mdp.write_text(
        "integrator = md\nnsteps = 0\ncontinuation = yes\nconstraints = none\n"
        "cutoff-scheme = Verlet\nnstlist = 20\nrlist = 1.0\n"
        "coulombtype = PME\nrcoulomb = 1.0\newald-rtol = 1e-5\n"
        "vdwtype = Cut-off\nrvdw = 1.0\nvdw-modifier = None\n"
        f"coulomb-modifier = None\nDispCorr = {'EnerPres' if dispersion else 'no'}\n"
        "pbc = xyz\nnstenergy = 1\n"
        + ("periodic-molecules = yes\n" if periodic else "periodic-molecules = no\n"),
        encoding="utf-8",
    )
    _run(
        [
            gmx, "grompp", "-f", mdp.name, "-p", "system.top", "-c", "structure.gro",
            "-o", "cross_engine.tpr", "-po", "cross_engine_mdout.mdp", "-maxwarn", "0",
        ],
        directory,
    )
    _run([gmx, "mdrun", "-s", "cross_engine.tpr", "-deffnm", "cross_engine"], directory)
    names = [
        "Bond", "Angle", "Proper Dih.", "Per. Imp. Dih.", "LJ-14", "Coulomb-14",
        "LJ (SR)", "Coulomb (SR)", "Coul. recip.", "Potential",
    ]
    if dispersion:
        names.insert(-1, "Disper. corr.")
    output = _run_with_input(
        [gmx, "energy", "-f", "cross_engine.edr", "-o", "cross_engine.xvg"],
        directory,
        "\n".join(name.replace(" ", "-") for name in names) + "\n0\n",
    )
    values: dict[str, float] = {"Disper. corr.": 0.0}
    for name in names:
        match = re.search(rf"^{re.escape(name)}\s+([-+0-9.eE]+)", output, re.MULTILINE)
        if match is None:
            raise RuntimeError(f"GROMACS did not report {name}")
        values[name] = float(match.group(1))
    return values


def _gro_coordinates(path: Path) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    atom_count = int(lines[1])
    coordinates = np.asarray(
        [
            [float(line[20:28]), float(line[28:36]), float(line[36:44])]
            for line in lines[2 : 2 + atom_count]
        ],
        dtype=np.float64,
    )
    values = [float(value) for value in lines[2 + atom_count].split()]
    if len(values) == _ORTHOGONAL_GRO_BOX_FIELDS:
        box = np.diag(values)
    elif len(values) == _TRICLINIC_GRO_BOX_FIELDS:
        box = np.asarray(
            [
                [values[0], values[3], values[4]],
                [values[5], values[1], values[6]],
                [values[7], values[8], values[2]],
            ],
            dtype=np.float64,
        )
    else:
        raise ValueError("unsupported GRO box-vector record")
    return coordinates, box


def _explicit_pair_energies(
    topology: Path, coordinates_path: Path, temporary_topology: Path
) -> dict[str, float]:
    _, pairs = prepare_explicit_pairs(topology, temporary_topology)
    coordinates, box = _gro_coordinates(coordinates_path)
    inverse = np.linalg.inv(box)
    lj = 0.0
    coulomb = 0.0
    for pair in pairs:
        delta = coordinates[pair.second] - coordinates[pair.first]
        delta -= np.rint(delta @ inverse) @ box
        distance_angstrom = float(np.linalg.norm(delta)) * 10.0
        sigma_angstrom = pair.sigma_nm * 10.0
        sr6 = (sigma_angstrom / distance_angstrom) ** 6
        lj += 4.0 * pair.epsilon_kj_mol / _KJ_PER_KCAL * (sr6 * sr6 - sr6)
        coulomb += _LAMMPS_COULOMB * pair.charge_product_e2 / distance_angstrom
    return {"lj_kcal_mol": lj, "coulomb_kcal_mol": coulomb, "pair_count": len(pairs)}


def _lammps_energies(lammps: str, directory: Path) -> dict[str, float]:
    output = _run([lammps, "-in", "in.chitin"], directory)
    rows = [line.split() for line in output.splitlines() if re.match(r"^\s*0\s+", line)]
    if not rows:
        raise RuntimeError("LAMMPS did not report its run-0 energy row")
    values = [float(value) for value in rows[-1]]
    return dict(
        zip(
            (
                "step", "temperature", "potential", "bonds", "angles", "torsions",
                "lj", "coulomb_short", "coulomb_long", "dispersion_tail",
            ),
            values,
            strict=True,
        )
    )


def _relative(first: float, second: float) -> float:
    return abs(first - second) / max(abs(first), 1.0)


def _case(  # noqa: PLR0913, PLR0917
    directory: Path,
    allomorph: str,
    shape: str,
    boundary: str,
    tleap: str,
    parmed_python: str,
    openmm_python: str,
    gmx_double: str,
    lammps: str,
) -> dict[str, object]:
    if directory.exists():
        shutil.rmtree(directory)
    arguments = _build_arguments(
        allomorph, shape, boundary, directory, tleap, parmed_python, "gromacs"
    )
    arguments.extend(
        (
            "--md-engine", "openmm", "--md-engine", "lammps",
            "--openmm-python", openmm_python,
        )
    )
    assert main(arguments) == 0
    report = json.loads((directory / "build_report.json").read_text(encoding="utf-8"))
    dispersion = shape == "bulk"
    gromacs = _gromacs_energies(
        gmx_double, directory, periodic=boundary == "periodic", dispersion=dispersion
    )
    openmm = cast(dict[str, float], report["openmm_validation"]["initial_energies_kj_mol"])
    lammps_values = _lammps_energies(lammps, directory)
    openmm_dynamics = _run(
        [
            openmm_python,
            "run_openmm.py",
            "--platform",
            "CPU",
            "--steps",
            "1000",
            "--report-interval",
            "1000",
        ],
        directory,
    )
    if "Final potential energy:" not in openmm_dynamics:
        raise RuntimeError("OpenMM did not complete its two-ps NVT smoke test")
    explicit = _explicit_pair_energies(
        directory / "system.top", directory / "structure.gro", directory / "pair_parser.top"
    )
    (directory / "pair_parser.top").unlink()
    explicit_lj = explicit["lj_kcal_mol"]
    explicit_coulomb = explicit["coulomb_kcal_mol"]

    gromacs_categories = {
        "bonds": gromacs["Bond"],
        "angles": gromacs["Angle"],
        "torsions": gromacs["Proper Dih."] + gromacs["Per. Imp. Dih."],
        "one_four": gromacs["LJ-14"] + gromacs["Coulomb-14"],
        "nonbonded": sum(
            gromacs[name]
            for name in (
                "LJ (SR)", "Coulomb (SR)", "Coul. recip.", "Disper. corr.",
                "LJ-14", "Coulomb-14",
            )
        ),
    }
    openmm_categories = {
        "bonds": openmm["HarmonicBondForce"],
        "angles": openmm["HarmonicAngleForce"],
        "torsions": openmm["PeriodicTorsionForce"],
        "one_four": openmm["NonbondedForce exceptions (1-4)"],
        "nonbonded": openmm["NonbondedForce"],
    }
    openmm_comparison = {
        name: {
            "gromacs_kj_mol": value,
            "openmm_kj_mol": openmm_categories[name],
            "relative_difference": _relative(value, openmm_categories[name]),
        }
        for name, value in gromacs_categories.items()
    }

    gromacs_kcal = {name: value / _KJ_PER_KCAL for name, value in gromacs.items()}
    lammps_reference = {
        "bonds_plus_explicit_pairs": gromacs_kcal["Bond"] + explicit_lj + explicit_coulomb,
        "angles": gromacs_kcal["Angle"],
        "torsions": gromacs_kcal["Proper Dih."] + gromacs_kcal["Per. Imp. Dih."],
        "remaining_lj": (
            gromacs_kcal["LJ (SR)"] + gromacs_kcal["LJ-14"] - explicit_lj
        ),
        "remaining_coulomb": (
            gromacs_kcal["Coulomb (SR)"] + gromacs_kcal["Coul. recip."]
            + gromacs_kcal["Coulomb-14"] - explicit_coulomb
        ),
        "potential": gromacs_kcal["Potential"],
        "dispersion_correction": gromacs_kcal["Disper. corr."],
    }
    lammps_observed = {
        "bonds_plus_explicit_pairs": lammps_values["bonds"],
        "angles": lammps_values["angles"],
        "torsions": lammps_values["torsions"],
        "remaining_lj": lammps_values["lj"] - lammps_values["dispersion_tail"],
        "remaining_coulomb": lammps_values["coulomb_short"] + lammps_values["coulomb_long"],
        "potential": lammps_values["potential"],
        "dispersion_correction": lammps_values["dispersion_tail"],
    }
    lammps_comparison = {
        name: {
            "gromacs_kcal_mol": value,
            "lammps_kcal_mol": lammps_observed[name],
            "relative_difference": _relative(value, lammps_observed[name]),
        }
        for name, value in lammps_reference.items()
    }
    maximum = max(
        row["relative_difference"]
        for comparison in (openmm_comparison, lammps_comparison)
        for name, row in comparison.items()
        if name != "dispersion_correction"
    )
    maximum_bonded = max(
        *(
            openmm_comparison[name]["relative_difference"]
            for name in ("bonds", "angles", "torsions", "one_four")
        ),
        *(
            lammps_comparison[name]["relative_difference"]
            for name in ("bonds_plus_explicit_pairs", "angles", "torsions", "remaining_lj")
        ),
    )
    if maximum_bonded > _MAXIMUM_BONDED_RELATIVE_DIFFERENCE:
        raise RuntimeError(
            f"bonded or explicit-pair energy difference {maximum_bonded:.3g} exceeds tolerance"
        )
    if maximum > _MAXIMUM_RELATIVE_DIFFERENCE:
        raise RuntimeError(f"cross-engine energy difference {maximum:.3g} exceeds tolerance")
    return {
        "allomorph": allomorph,
        "shape": shape,
        "chain_boundary": boundary,
        "atom_count": report["force_field_atom_count"],
        "explicit_pairs": explicit,
        "openmm_vs_gromacs": openmm_comparison,
        "lammps_vs_gromacs": lammps_comparison,
        "maximum_relative_difference": maximum,
        "maximum_bonded_and_explicit_pair_relative_difference": maximum_bonded,
        "openmm_nvt_completed_ps": 2.0,
        "status": "passed",
    }


def main_cli() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--openmm-python", required=True)
    parser.add_argument("--gmx-double", required=True)
    parser.add_argument("--lammps", default=shutil.which("lmp") or shutil.which("lmp_serial"))
    parser.add_argument(
        "--output",
        type=Path,
        default=_ROOT / "validation" / "results" / "chitin_glycam_cross_engines.json",
    )
    args = parser.parse_args()
    if args.lammps is None:
        raise FileNotFoundError("LAMMPS executable is required")
    tleap = shutil.which("tleap") or "tleap"
    parmed_python = find_parmed_python()
    requested_cases = [
        (allomorph, shape, boundary)
        for allomorph in ("alpha", "beta")
        for shape, boundary in (("slab", "finite"), ("slab", "periodic"), ("bulk", "periodic"))
    ]
    cases = [
        _case(
            args.work_directory / f"{allomorph}_{shape}_{boundary}",
            allomorph,
            shape,
            boundary,
            tleap,
            parmed_python,
            args.openmm_python,
            args.gmx_double,
            args.lammps,
        )
        for allomorph, shape, boundary in requested_cases
    ]
    result = {
        "provenance": provenance(("numpy",)),
        "scope": "GLYCAM06j-1 chitin term-by-term cross-engine validation",
        "force_field_accuracy_assessed": False,
        "cases": cases,
        "status": "passed",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main_cli()
