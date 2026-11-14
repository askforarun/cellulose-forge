"""Validate the supported GLYCAM06j-1 and CHARMM36 cellulose matrix.

The matrix contains periodic bulk crystals for all four allomorphs and finite
and covalently periodic slabs for the documented chain-parallel faces.  Every
case is built from the public CLI, preprocessed with ``grompp -maxwarn 0``,
minimized, and propagated for a 1 ps NVT smoke test.  Periodic builds must also
carry a successful mapped-geometry and written-topology O4--C1 bond check in
their build report.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from _provenance import provenance

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "validation" / "results" / "force_field_matrix.json"

BULK_CELLS: dict[str, tuple[int, int, int]] = {
    "I-alpha": (4, 4, 4),
    "I-beta": (4, 4, 4),
    "II": (4, 4, 4),
    # III-I has a short a axis and needs six repeats for the 1.0 nm cutoff.
    "III-I": (6, 4, 4),
}
SLAB_FACES: dict[str, tuple[int, int, int]] = {
    "I-alpha": (0, 1, 0),
    "I-beta": (1, 1, 0),
    "II": (1, 1, 0),
    "III-I": (1, 1, 0),
}
FORCE_FIELDS = ("glycam", "charmm36")
PERIODIC_BOND_LIMIT_ANGSTROM = 2.0


@dataclass(frozen=True)
class MatrixCase:
    """One force-field/geometry/boundary entry in the validation matrix."""

    force_field: str
    allomorph: str
    shape: str
    boundary: str


def _slug(value: str) -> str:
    return value.lower().replace("-", "_")


def _run(
    command: list[str | Path],
    *,
    cwd: Path,
    log_path: Path,
    stdin: str | None = None,
) -> str:
    result = subprocess.run(
        [str(item) for item in command],
        cwd=cwd,
        input=stdin,
        text=True,
        capture_output=True,
        check=False,
    )
    output = result.stdout + result.stderr
    log_path.write_text(output, encoding="utf-8")
    if result.returncode != 0:
        rendered = " ".join(str(item) for item in command)
        raise RuntimeError(f"command failed ({result.returncode}): {rendered}\n{output}")
    return output


def _engine_version(executable: str) -> str:
    result = subprocess.run([executable, "--version"], text=True, capture_output=True, check=True)
    first = (result.stdout + result.stderr).splitlines()[0]
    return first.strip()


def _external_python_package_version(executable: Path, package: str) -> str:
    result = subprocess.run(
        [
            str(executable),
            "-c",
            (f"from importlib.metadata import version; print(version({package!r}))"),
        ],
        text=True,
        capture_output=True,
        check=True,
    )
    return result.stdout.strip()


def _build_command(
    args: argparse.Namespace,
    *,
    case: MatrixCase,
    output: Path,
) -> list[str | Path]:
    command: list[str | Path] = [
        args.python,
        "-m",
        "cellulose_forge.cli",
        "build",
        "--polymer",
        "cellulose",
        "--allomorph",
        case.allomorph,
        "--shape",
        case.shape,
        "--chain-boundary",
        case.boundary,
        "--force-field",
        case.force_field,
        "--output",
        output,
    ]
    if case.shape == "bulk":
        command.extend(["--cells", *(str(value) for value in BULK_CELLS[case.allomorph])])
    else:
        command.extend(
            [
                "--face",
                *(str(value) for value in SLAB_FACES[case.allomorph]),
                "--layers",
                "3",
                "--chains-per-layer",
                "6",
                "--degree-polymerization",
                "6",
                "--vacuum",
                "24",
            ]
        )
        if case.boundary == "finite":
            command.extend(["--padding-longitudinal", "12"])
    if case.force_field == "glycam":
        command.extend(["--tleap-executable", args.tleap_executable])
        if args.parmed_python is not None:
            command.extend(["--parmed-python", args.parmed_python])
    else:
        command.extend(["--charmm-force-field-directory", args.charmm_force_field_directory])
    return command


def _periodic_bond_check(report: dict[str, Any], *, periodic: bool) -> dict[str, Any] | None:
    value = report.get("periodic_glycosidic_bond_validation")
    if not periodic:
        if value is not None:
            raise ValueError("finite build unexpectedly reports periodic glycosidic bonds")
        return None
    if not isinstance(value, dict):
        raise ValueError("periodic build lacks glycosidic-bond validation")
    for key in ("mapped_geometry", "written_topology"):
        check = value.get(key)
        if not isinstance(check, dict) or int(check.get("bond_count", 0)) <= 0:
            raise ValueError(f"periodic {key} bond validation is missing or empty")
        if float(check.get("maximum_distance_angstrom", 99.0)) >= PERIODIC_BOND_LIMIT_ANGSTROM:
            raise ValueError(f"periodic {key} O4--C1 bond reaches 2 angstrom")
    return value


def _parse_minimization(output: str) -> dict[str, float | int]:
    steps = re.search(r"converged to Fmax <\s*[-+0-9.eE]+ in (\d+) steps", output)
    force = re.search(r"Maximum force\s*=\s*([-+0-9.eE]+)", output)
    potential = re.search(r"Potential Energy\s*=\s*([-+0-9.eE]+)", output)
    if steps is None or force is None or potential is None:
        raise ValueError("could not parse converged GROMACS minimization")
    return {
        "steps": int(steps.group(1)),
        "maximum_force_kj_mol_nm": float(force.group(1)),
        "potential_energy_kj_mol": float(potential.group(1)),
    }


def _check_nvt(output: str, log: str, *, expected_steps: int) -> None:
    text = (output + "\n" + log).lower()
    found: list[str] = []
    for label, pattern in (
        ("LINCS warning", r"\blincs warning\b"),
        ("fatal error", r"\bfatal error\b"),
        ("NaN", r"\bnan\b"),
        ("infinite", r"\binfinite\b"),
    ):
        if re.search(pattern, text):
            found.append(label)
    if found:
        raise ValueError(f"NVT output contains failure marker(s): {', '.join(found)}")
    if f"{expected_steps} steps" not in text and "finished mdrun" not in text:
        raise ValueError("NVT output does not record a completed mdrun")


def _validate_case(
    args: argparse.Namespace,
    *,
    case: MatrixCase,
) -> dict[str, Any]:
    name = "_".join((_slug(case.force_field), _slug(case.allomorph), case.shape, case.boundary))
    directory = args.work_directory / name
    directory.mkdir(parents=True, exist_ok=True)
    _run(
        _build_command(
            args,
            case=case,
            output=directory,
        ),
        cwd=ROOT,
        log_path=directory / "build.log",
    )
    report = json.loads((directory / "build_report.json").read_text(encoding="utf-8"))
    bond_check = _periodic_bond_check(report, periodic=case.boundary == "periodic")
    _run(
        [
            args.gmx,
            "grompp",
            "-f",
            "minim.mdp",
            "-c",
            "structure.gro",
            "-p",
            "system.top",
            "-o",
            "em.tpr",
            "-po",
            "em_out.mdp",
            "-maxwarn",
            "0",
        ],
        cwd=directory,
        log_path=directory / "grompp_em.log",
    )
    em = _run(
        [args.gmx, "mdrun", "-deffnm", "em", "-nt", str(args.threads)],
        cwd=directory,
        log_path=directory / "mdrun_em.log",
    )
    minimization = _parse_minimization(em)
    _run(
        [
            args.gmx,
            "grompp",
            "-f",
            "nvt.mdp",
            "-c",
            "em.gro",
            "-r",
            "em.gro",
            "-p",
            "system.top",
            "-o",
            "nvt.tpr",
            "-po",
            "nvt_out.mdp",
            "-maxwarn",
            "0",
        ],
        cwd=directory,
        log_path=directory / "grompp_nvt.log",
    )
    nvt = _run(
        [args.gmx, "mdrun", "-deffnm", "nvt", "-nt", str(args.threads)],
        cwd=directory,
        log_path=directory / "mdrun_nvt.log",
    )
    nvt_log = (directory / "nvt.log").read_text(encoding="utf-8", errors="replace")
    _check_nvt(nvt, nvt_log, expected_steps=1_000)
    return {
        "id": name,
        "force_field": ("GLYCAM06j-1" if case.force_field == "glycam" else "CHARMM36"),
        "allomorph": case.allomorph,
        "shape": case.shape,
        "chain_boundary": case.boundary,
        "face": (list(SLAB_FACES[case.allomorph]) if case.shape == "slab" else None),
        "cells": (list(BULK_CELLS[case.allomorph]) if case.shape == "bulk" else None),
        "atom_count": report["force_field_atom_count"],
        "chain_count": report["chain_count"],
        "degree_polymerization": report["degree_polymerization"],
        "periodic_glycosidic_bonds": bond_check,
        "minimization": minimization,
        "nvt": {"duration_ps": 1.0, "steps": 1_000, "status": "completed"},
        "status": "passed",
    }


def main() -> int:
    """Run and record the 24-system force-field validation matrix."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--tleap-executable", default="tleap")
    parser.add_argument("--parmed-python", type=Path)
    parser.add_argument("--charmm-force-field-directory", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--work-directory",
        type=Path,
        default=ROOT / "validation" / "work" / "force_field_matrix",
    )
    parser.add_argument("--output", type=Path, default=RESULT)
    args = parser.parse_args()
    if args.threads <= 0:
        parser.error("--threads must be positive")
    args.work_directory.mkdir(parents=True, exist_ok=True)
    matrix = [
        MatrixCase(force_field, allomorph, "bulk", "periodic")
        for force_field in FORCE_FIELDS
        for allomorph in BULK_CELLS
    ]
    matrix.extend(
        MatrixCase(force_field, allomorph, "slab", boundary)
        for force_field in FORCE_FIELDS
        for allomorph in SLAB_FACES
        for boundary in ("finite", "periodic")
    )
    cases = [_validate_case(args, case=case) for case in matrix]
    result_provenance = provenance(("cellulose-forge", "numpy"))
    if args.parmed_python is None:
        raise ValueError("--parmed-python is required for provenance-bearing GLYCAM validation")
    software_versions = result_provenance["software_versions"]
    assert isinstance(software_versions, dict)
    software_versions["parmed"] = _external_python_package_version(args.parmed_python, "ParmEd")
    result = {
        "validation": "GLYCAM06j-1 and CHARMM36 cellulose force-field matrix",
        "provenance": result_provenance,
        "gromacs_version": _engine_version(args.gmx),
        "protocol": {
            "bulk": "four periodic allomorphs with each force field",
            "slabs": "four faces, finite and periodic chains, with each force field",
            "preprocessing": "grompp -maxwarn 0",
            "minimization": "generated force-field-consistent minim.mdp",
            "dynamics": "generated 1 ps PME/NVT smoke test at 300 K",
            "periodic_bond_limit_angstrom": PERIODIC_BOND_LIMIT_ANGSTROM,
        },
        "case_count": len(cases),
        "passed_count": sum(case["status"] == "passed" for case in cases),
        "cases": cases,
        "status": "passed" if all(case["status"] == "passed" for case in cases) else "failed",
        "interpretation_limit": (
            "These runs validate construction, preprocessing, topology integrity, "
            "and short-run numerical stability; they do not establish equilibrated "
            "material properties or force-field accuracy."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
