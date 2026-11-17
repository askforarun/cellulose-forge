"""Validate serialized OpenMM systems against GROMACS and short OpenMM dynamics.

For each case the script builds GROMACS and OpenMM outputs from the same mapped
structure, compares OpenMM and GROMACS zero-step energies by interaction
category, runs the generated OpenMM runner for a short NVT simulation, and,
for covalently periodic I-beta slabs, a short NPT simulation with the chain
axis held fixed. Bond lengths are rechecked against their equilibrium values
afterwards.

Bonded and 1-4 categories are force-field definitions and should agree to
GROMACS output precision. The total nonbonded energy is compared with GROMACS
run with matching settings (plain cutoffs without potential shifts, the same
PME accuracy, and the same dispersion correction); it differs slightly because
the engines choose different PME grids.
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
sys.path.insert(0, str(VALIDATION))

from _provenance import provenance  # noqa: E402

from cellulose_forge.energy_drivers import (  # noqa: E402
    EnergyReport,
    compare_reports,
    run_gromacs_single_point,
    validate_potential_sum,
)
from cellulose_forge.energy_drivers.coordinates import (  # noqa: E402
    read_gro_frame,
    write_coordinate_reference,
)
from cellulose_forge.energy_drivers.report import energy_report_from_dict  # noqa: E402

SLAB = [
    "--shape",
    "slab",
    "--face",
    "1",
    "1",
    "0",
    "--layers",
    "3",
    "--chains-per-layer",
    "4",
    "--degree-polymerization",
    "8",
    "--vacuum",
    "30",
]
ALLOMORPH_SLAB_FACES = {"I-alpha": "0 1 0", "II": "1 1 0", "III-I": "1 1 0"}
CASES: dict[str, dict[str, Any]] = {
    "glycam_finite_slab": {
        "force_field": "glycam",
        "chain_boundary": "finite",
        "geometry": [*SLAB, "--padding-longitudinal", "12"],
    },
    "glycam_periodic_slab": {
        "force_field": "glycam",
        "chain_boundary": "periodic",
        "geometry": SLAB,
        "npt_axes": "x",
    },
    "charmm36_finite_slab": {
        "force_field": "charmm36",
        "chain_boundary": "finite",
        "geometry": [*SLAB, "--padding-longitudinal", "12"],
    },
    "charmm36_periodic_slab": {
        "force_field": "charmm36",
        "chain_boundary": "periodic",
        "geometry": SLAB,
        "npt_axes": "x",
    },
    "glycam_periodic_triclinic_bulk": {
        "force_field": "glycam",
        "chain_boundary": "periodic",
        "geometry": ["--shape", "bulk", "--cells", "4", "4", "4"],
        "dispersion_correction": True,
    },
}
for _allomorph, _face in ALLOMORPH_SLAB_FACES.items():
    for _force_field in ("glycam", "charmm36"):
        CASES[f"{_force_field}_periodic_{_allomorph.lower().replace('-', '_')}_slab"] = {
            "allomorph": _allomorph,
            "force_field": _force_field,
            "chain_boundary": "periodic",
            "geometry": [
                "--shape",
                "slab",
                "--face",
                *_face.split(),
                "--layers",
                "3",
                "--chains-per-layer",
                "6",
                "--degree-polymerization",
                "8",
                "--vacuum",
                "30",
            ],
        }
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


def _gromacs_energies(  # noqa: PLR0913 - scientific inputs remain explicit
    gmx: str,
    directory: Path,
    *,
    periodic: bool,
    dispersion_correction: bool,
    coordinate_reference: Path,
    coordinate_sha256: str,
) -> EnergyReport:
    mdp = directory / "single_point.mdp"
    mdp.write_text(
        "integrator = md\nnsteps = 0\nconstraints = none\ncutoff-scheme = Verlet\n"
        "coulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\nvdw-modifier = None\n"
        "coulomb-modifier = None\newald-rtol = 1e-5\n"
        f"DispCorr = {'EnerPres' if dispersion_correction else 'no'}\n"
        "pbc = xyz\nnstenergy = 1\n" + ("periodic-molecules = yes\n" if periodic else ""),
        encoding="utf-8",
    )
    return run_gromacs_single_point(
        directory,
        coordinate_file=directory / "structure.gro",
        coordinate_reference_file=coordinate_reference,
        topology_file=directory / "system.top",
        mdp_file=mdp,
        coordinate_sha256=coordinate_sha256,
        executable=gmx,
    )


def _openmm_energies(
    args: argparse.Namespace,
    directory: Path,
    coordinate_reference: Path,
    coordinate_sha256: str,
) -> EnergyReport:
    evaluator = VALIDATION / "evaluate_openmm_single_point.py"
    stdout = _run(
        [
            str(args.openmm_python),
            str(evaluator),
            "--system",
            "validation_system.xml",
            "--coordinates",
            "structure.gro",
            "--coordinate-reference",
            str(coordinate_reference),
            "--coordinate-sha256",
            coordinate_sha256,
            "--topology",
            "system.top",
            "--include-directory",
            str(directory),
            "--platform",
            args.platform,
        ],
        directory,
    )
    payload = json.loads(stdout.strip().splitlines()[-1])
    return energy_report_from_dict(
        payload,
        expected_engine="openmm",
        coordinate_sha256=coordinate_sha256,
        required_terms=set(payload["terms"]),
    )


def _runner(args: argparse.Namespace, directory: Path, *extra: str) -> dict[str, Any]:
    stdout = _run(
        [
            str(args.openmm_python),
            "run_openmm.py",
            "--platform",
            args.platform,
            "--minimize-iterations",
            "200",
            "--steps",
            str(args.steps),
            "--report-interval",
            str(args.steps // 5),
            *extra,
        ],
        directory,
    )
    match = re.search(r"Final potential energy: ([-+0-9.eE]+)", stdout)
    rows = [
        line.split(",")
        for line in (directory / "state_data.csv").read_text(encoding="utf-8").splitlines()[1:]
    ]
    bonds = json.loads(
        _run([str(args.openmm_python), "-c", BOND_CHECK, "final_state.xml"], directory)
    )
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
        "--allomorph",
        case.get("allomorph", "I-beta"),
        *case["geometry"],
        "--chain-boundary",
        case["chain_boundary"],
        "--force-field",
        case["force_field"],
        "--md-engine",
        "gromacs",
        "--md-engine",
        "openmm",
    ]
    command = [
        sys.executable,
        "-m",
        "cellulose_forge.cli",
        "build",
        *build_arguments,
        "--openmm-python",
        str(args.openmm_python),
        "--tleap-executable",
        args.tleap_executable,
        "--output",
        str(directory),
    ]
    if args.parmed_python:
        command.extend(["--parmed-python", str(args.parmed_python)])
    _run(command, ROOT)
    report = json.loads((directory / "build_report.json").read_text(encoding="utf-8"))
    openmm = report["openmm_validation"]
    periodic = case["chain_boundary"] == "periodic"
    positions, box = read_gro_frame(directory / "structure.gro")
    coordinate_reference = directory / "coordinate_reference.json"
    coordinate_sha256 = write_coordinate_reference(coordinate_reference, positions, box)
    gromacs_report = _gromacs_energies(
        args.gmx,
        directory,
        periodic=periodic,
        dispersion_correction=openmm["dispersion_correction"],
        coordinate_reference=coordinate_reference,
        coordinate_sha256=coordinate_sha256,
    )
    openmm_report = _openmm_energies(args, directory, coordinate_reference, coordinate_sha256)
    for energy_report in (gromacs_report, openmm_report):
        validate_potential_sum(energy_report, absolute_kj_mol=0.1, relative=1.0e-6)
    energy_comparison = compare_reports(
        gromacs_report,
        openmm_report,
        absolute_kj_mol=0.1,
        relative=2.0e-4,
    )
    result: dict[str, Any] = {
        "allomorph": case.get("allomorph", "I-beta"),
        "force_field": case["force_field"],
        "chain_boundary": case["chain_boundary"],
        "dispersion_correction": openmm["dispersion_correction"],
        "build_arguments": build_arguments,
        "atom_count": openmm["particle_count"],
        "periodic_forces": openmm["periodic_forces"],
        "chain_axis": openmm["chain_axis"],
        "initial_maximum_bond_length_deviation_nm": openmm["maximum_bond_length_deviation_nm"],
        "energy_reports": {
            "gromacs": gromacs_report.to_dict(),
            "openmm": openmm_report.to_dict(),
        },
        "energy_comparison": energy_comparison,
        "nvt": _runner(args, directory),
    }
    if "npt_axes" in case:
        result["npt_fixed_chain_axis"] = _runner(
            args, directory, "--pressure-bar", "1", "--pressure-axes", case["npt_axes"]
        )
    return result


def _energy_comparisons_passed(results: dict[str, dict[str, Any]]) -> bool:
    """Return whether every selected case passed its energy comparison."""
    return all(case["energy_comparison"]["passed"] for case in results.values())


def main() -> int:
    """Run the OpenMM validation matrix and write a JSON summary."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--openmm-python", type=Path, required=True)
    parser.add_argument("--tleap-executable", default="tleap")
    parser.add_argument("--parmed-python", type=Path)
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--steps", type=int, default=500)
    parser.add_argument(
        "--case",
        action="append",
        dest="selected_cases",
        choices=sorted(CASES),
        help="run only this validation case; repeat to select more than one",
    )
    parser.add_argument(
        "--work-directory", type=Path, default=VALIDATION / "work" / "openmm_engine"
    )
    parser.add_argument(
        "--output", type=Path, default=VALIDATION / "results" / "openmm_engine_validation.json"
    )
    args = parser.parse_args()
    args.work_directory.mkdir(parents=True, exist_ok=True)

    versions = json.loads(
        _run(
            [
                str(args.openmm_python),
                "-c",
                "import json, openmm; print(json.dumps(openmm.__version__))",
            ],
            ROOT,
        )
    )
    gromacs_version = re.search(r"GROMACS version:\s+(\S+)", _run([args.gmx, "--version"], ROOT))
    commit = _run(["git", "rev-parse", "HEAD"], ROOT).strip()
    dirty = bool(_run(["git", "status", "--porcelain", "--untracked-files=no"], ROOT).strip())
    selected = set(args.selected_cases or CASES)
    results = {
        name: _validate_case(args, name, case) for name, case in CASES.items() if name in selected
    }
    summary = {
        "provenance": provenance(("cellulose-forge", "numpy")),
        "cellulose_forge_commit": commit,
        "cellulose_forge_working_tree_modified": dirty,
        "openmm_version": versions,
        "gromacs_version": gromacs_version.group(1) if gromacs_version else None,
        "platform": args.platform,
        "timestep_fs": 2.0,
        "temperature_k": 300.0,
        "cases": results,
        "maximum_energy_relative_difference": max(
            term["relative_difference"]
            for case in results.values()
            for term in case["energy_comparison"]["common_terms"].values()
        ),
    }
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                name: {
                    "energy_comparison_passed": case["energy_comparison"]["passed"],
                    "maximum_energy_relative_difference": max(
                        term["relative_difference"]
                        for term in case["energy_comparison"]["common_terms"].values()
                    ),
                    "nvt_T": case["nvt"]["final_temperature_k"],
                    "nvt_max_bond_dev_nm": case["nvt"]["maximum_bond_length_deviation_nm"],
                    "npt_max_bond_dev_nm": case.get("npt_fixed_chain_axis", {}).get(
                        "maximum_bond_length_deviation_nm"
                    ),
                }
                for name, case in results.items()
            },
            indent=1,
        )
    )
    return 0 if _energy_comparisons_passed(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
