"""Validate serialized OpenMM systems against GROMACS and short OpenMM dynamics.

For each case the script builds GROMACS and OpenMM outputs from the same mapped
structure, compares OpenMM and GROMACS zero-step energies by interaction
category, runs the generated OpenMM runner for a short NVT simulation, and,
for covalently periodic slabs, a short NPT simulation with the chain axis held
fixed. Bond lengths are rechecked against their equilibrium values afterwards.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "validation"

SLAB = [
    "--shape", "slab", "--face", "1", "1", "0",
    "--layers", "3", "--chains-per-layer", "4",
    "--degree-polymerization", "8", "--vacuum", "30",
]
CASES: dict[str, dict[str, Any]] = {
    "glycam_finite_slab": {
        "force_field": "glycam", "chain_boundary": "finite",
        "geometry": [*SLAB, "--padding-longitudinal", "12"],
    },
    "glycam_periodic_slab": {
        "force_field": "glycam", "chain_boundary": "periodic", "geometry": SLAB,
        "npt_axes": "x",
    },
    "charmm36_finite_slab": {
        "force_field": "charmm36", "chain_boundary": "finite",
        "geometry": [*SLAB, "--padding-longitudinal", "12"],
    },
    "charmm36_periodic_slab": {
        "force_field": "charmm36", "chain_boundary": "periodic", "geometry": SLAB,
        "npt_axes": "x",
    },
    "glycam_periodic_triclinic_bulk": {
        "force_field": "glycam", "chain_boundary": "periodic",
        "geometry": ["--shape", "bulk", "--cells", "3", "3", "2"],
    },
}
GROMACS_TERMS = {
    "bonds_and_urey_bradley": ("Bond", "U-B"),
    "angles": ("Angle",),
    "torsions": ("Proper Dih.", "Improper Dih.", "Per. Imp. Dih.", "CMAP Dih."),
    "one_four": ("LJ-14", "Coulomb-14"),
}
OPENMM_TERMS = {
    "bonds_and_urey_bradley": ("HarmonicBondForce",),
    "angles": ("HarmonicAngleForce",),
    "torsions": ("PeriodicTorsionForce", "CustomTorsionForce", "CMAPTorsionForce"),
    "one_four": ("NonbondedForce exceptions (1-4)",),
}
# CHARMM36 Urey-Bradley terms are one GROMACS angle type but two OpenMM forces.
COMBINED_UREY_BRADLEY = ("bonds_and_urey_bradley", "angles")

BOND_CHECK = r"""
import json, sys
import numpy, openmm
from openmm import unit
system = openmm.XmlSerializer.deserialize(open("system.xml").read())
state = openmm.XmlSerializer.deserialize(open(sys.argv[1]).read())
positions = numpy.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
box = numpy.asarray([[c for c in v.value_in_unit(unit.nanometer)]
                     for v in state.getPeriodicBoxVectors()])
worst = 0.0
for force in system.getForces():
    if isinstance(force, openmm.HarmonicBondForce):
        for index in range(force.getNumBonds()):
            first, second, length, _ = force.getBondParameters(index)
            delta = positions[second] - positions[first]
            for axis in (2, 1, 0):
                delta = delta - box[axis] * round(delta[axis] / box[axis][axis])
            worst = max(worst, abs(float(numpy.linalg.norm(delta))
                                   - length.value_in_unit(unit.nanometer)))
print(json.dumps({"maximum_bond_length_deviation_nm": worst,
                  "box_lengths_nm": [float(box[i][i]) for i in range(3)]}))
"""


def _run(command: list[str], cwd: Path, input_text: str | None = None) -> str:
    completed = subprocess.run(
        command, cwd=cwd, input=input_text, check=False, capture_output=True, text=True
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed: {' '.join(command)}\n{completed.stdout}\n{completed.stderr}"
        )
    return completed.stdout


def _gromacs_energies(gmx: str, directory: Path, *, periodic: bool) -> dict[str, float]:
    mdp = directory / "single_point.mdp"
    mdp.write_text(
        "integrator = md\nnsteps = 0\nconstraints = none\ncutoff-scheme = Verlet\n"
        "coulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\nnstenergy = 1\n"
        + ("periodic-molecules = yes\n" if periodic else ""),
        encoding="utf-8",
    )
    _run([gmx, "grompp", "-f", mdp.name, "-c", "structure.gro", "-p", "system.top",
          "-o", "single_point.tpr"], directory)
    _run([gmx, "mdrun", "-s", "single_point.tpr", "-deffnm", "single_point",
          "-ntmpi", "1", "-ntomp", "1"], directory)
    names = [term for terms in GROMACS_TERMS.values() for term in terms]
    stdout = _run(
        [gmx, "energy", "-f", "single_point.edr", "-o", "single_point.xvg"],
        directory,
        "\n".join(name.replace(" ", "-") for name in names) + "\n0\n",
    )
    values: dict[str, float] = {}
    for name in names:
        match = re.search(rf"^{re.escape(name)}\s+([-+0-9.eE]+)", stdout, re.MULTILINE)
        if match:
            values[name] = float(match.group(1))
    return values


def _compare(gromacs: dict[str, float], openmm: dict[str, float]) -> dict[str, Any]:
    rows: dict[str, dict[str, float]] = {}
    for category, gromacs_names in GROMACS_TERMS.items():
        rows[category] = {
            "gromacs_kj_mol": sum(gromacs.get(name, 0.0) for name in gromacs_names),
            "openmm_kj_mol": sum(openmm.get(name, 0.0) for name in OPENMM_TERMS[category]),
        }
    if "U-B" in gromacs:
        merged = {
            key: sum(rows[category][key] for category in COMBINED_UREY_BRADLEY)
            for key in ("gromacs_kj_mol", "openmm_kj_mol")
        }
        for category in COMBINED_UREY_BRADLEY:
            del rows[category]
        rows["bonds_angles_and_urey_bradley"] = merged
    for row in rows.values():
        reference = abs(row["gromacs_kj_mol"]) or 1.0
        row["relative_difference"] = abs(row["openmm_kj_mol"] - row["gromacs_kj_mol"]) / reference
    return {
        "categories": rows,
        "maximum_relative_difference": max(row["relative_difference"] for row in rows.values()),
    }


def _runner(args: argparse.Namespace, directory: Path, *extra: str) -> dict[str, Any]:
    stdout = _run(
        [str(args.openmm_python), "run_openmm.py", "--platform", args.platform,
         "--minimize-iterations", "200", "--steps", str(args.steps),
         "--report-interval", str(args.steps // 5), *extra],
        directory,
    )
    match = re.search(r"Final potential energy: ([-+0-9.eE]+)", stdout)
    rows = [line.split(",") for line in (directory / "state_data.csv").read_text(
        encoding="utf-8").splitlines()[1:]]
    bonds = json.loads(_run([str(args.openmm_python), "-c", BOND_CHECK, "final_state.xml"],
                            directory))
    return {
        "arguments": list(extra),
        "steps": args.steps,
        "final_potential_energy_kj_mol": float(match.group(1)) if match else None,
        "final_temperature_k": float(rows[-1][3]),
        "final_volume_nm3": float(rows[-1][4]),
        **bonds,
    }


def _validate_case(args: argparse.Namespace, name: str, case: dict[str, Any]) -> dict[str, Any]:
    directory = args.work_directory / name
    build_arguments = [
        "--allomorph", "I-beta", *case["geometry"],
        "--chain-boundary", case["chain_boundary"],
        "--force-field", case["force_field"],
        "--md-engine", "gromacs", "--md-engine", "openmm",
    ]
    command = [
        sys.executable, "-m", "cellulose_forge.cli", "build", *build_arguments,
        "--openmm-python", str(args.openmm_python),
        "--tleap-executable", args.tleap_executable,
        "--output", str(directory),
    ]
    if args.parmed_python:
        command.extend(["--parmed-python", str(args.parmed_python)])
    _run(command, ROOT)
    report = json.loads((directory / "build_report.json").read_text(encoding="utf-8"))
    openmm = report["openmm_validation"]
    periodic = case["chain_boundary"] == "periodic"
    result: dict[str, Any] = {
        "force_field": case["force_field"],
        "chain_boundary": case["chain_boundary"],
        "build_arguments": build_arguments,
        "atom_count": openmm["particle_count"],
        "periodic_forces": openmm["periodic_forces"],
        "chain_axis": openmm["chain_axis"],
        "initial_maximum_bond_length_deviation_nm": openmm["maximum_bond_length_deviation_nm"],
        "energy_comparison": _compare(
            _gromacs_energies(args.gmx, directory, periodic=periodic),
            openmm["initial_energies_kj_mol"],
        ),
        "nvt": _runner(args, directory),
    }
    if "npt_axes" in case:
        result["npt_fixed_chain_axis"] = _runner(
            args, directory, "--pressure-bar", "1", "--pressure-axes", case["npt_axes"]
        )
    return result


def main() -> int:
    """Run the OpenMM validation matrix and write a JSON summary."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--openmm-python", type=Path, required=True)
    parser.add_argument("--tleap-executable", default="tleap")
    parser.add_argument("--parmed-python", type=Path)
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--steps", type=int, default=500)
    parser.add_argument("--work-directory", type=Path,
                        default=VALIDATION / "work" / "openmm_engine")
    parser.add_argument("--output", type=Path,
                        default=VALIDATION / "results" / "openmm_engine_validation.json")
    args = parser.parse_args()
    args.work_directory.mkdir(parents=True, exist_ok=True)

    versions = json.loads(_run([str(args.openmm_python), "-c",
                                "import json, openmm; print(json.dumps(openmm.__version__))"],
                               ROOT))
    gromacs_version = re.search(r"GROMACS version:\s+(\S+)", _run([args.gmx, "--version"], ROOT))
    results = {name: _validate_case(args, name, case) for name, case in CASES.items()}
    summary = {
        "openmm_version": versions,
        "gromacs_version": gromacs_version.group(1) if gromacs_version else None,
        "platform": args.platform,
        "timestep_fs": 2.0,
        "temperature_k": 300.0,
        "cases": results,
        "maximum_relative_energy_difference": max(
            case["energy_comparison"]["maximum_relative_difference"]
            for case in results.values()
        ),
    }
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({name: {
        "max_rel_energy_diff": case["energy_comparison"]["maximum_relative_difference"],
        "nvt_T": case["nvt"]["final_temperature_k"],
        "nvt_max_bond_dev_nm": case["nvt"]["maximum_bond_length_deviation_nm"],
        "npt_max_bond_dev_nm": case.get("npt_fixed_chain_axis", {}).get(
            "maximum_bond_length_deviation_nm"),
    } for name, case in results.items()}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
