"""Validate GLYCAM06j-1 alpha- and beta-chitin output with real engines.

The matrix covers periodic bulk crystals and finite/periodic (100) slabs with
GROMACS, plus a finite beta-chitin slab with native AmberTools ``sander``.
The one-ps NVT runs are smoke tests, not force-field-accuracy validation.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path

from _provenance import provenance

from cellulose_forge.cli import main
from cellulose_forge.glycam import find_parmed_python

_ROOT = Path(__file__).resolve().parents[1]
_PROBLEM = re.compile(r"WARNING|LINCS WARNING|\bnan\b|not finite", re.IGNORECASE)
_MINIMUM_PAIR_FIELDS = 3
_KJ_PER_KCAL = 4.184
_GROMACS_COULOMB_KCAL_ANGSTROM = 332.0637
_AMBER_COULOMB_KCAL_ANGSTROM = 332.0522
_ENERGY_DIFFERENCE_TOLERANCE_KCAL_MOL = 5.0e-3


def _run(command: list[str], directory: Path) -> str:
    result = subprocess.run(
        command,
        cwd=directory,
        capture_output=True,
        text=True,
        check=False,
    )
    output = result.stdout + result.stderr
    if result.returncode != 0:
        raise RuntimeError(f"{' '.join(command)} failed:\n{output}")
    return output


def _run_with_input(command: list[str], directory: Path, stdin: str) -> str:
    result = subprocess.run(
        command,
        cwd=directory,
        input=stdin,
        capture_output=True,
        text=True,
        check=False,
    )
    output = result.stdout + result.stderr
    if result.returncode != 0:
        raise RuntimeError(f"{' '.join(command)} failed:\n{output}")
    return output


def _external_versions(
    gmx: str,
    parmed_python: str,
    amber_output: Path,
) -> dict[str, str]:
    gromacs = subprocess.run(
        [gmx, "--version"], capture_output=True, text=True, check=False
    ).stdout
    gromacs_match = re.search(r"GROMACS version:\s*(.+)", gromacs)
    parmed = subprocess.run(
        [
            parmed_python,
            "-c",
            "from importlib.metadata import version; print(version('parmed'))",
        ],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    amber = amber_output.read_text(encoding="utf-8")
    amber_match = re.search(r"Amber\s+(\d+)\s+SANDER", amber)
    return {
        "gromacs": gromacs_match.group(1).strip() if gromacs_match else "unknown",
        "parmed": parmed,
        "amber_sander": amber_match.group(1) if amber_match else "unknown",
    }


def _build_arguments(  # noqa: PLR0913, PLR0917
    allomorph: str,
    shape: str,
    boundary: str,
    output: Path,
    tleap: str,
    parmed_python: str,
    engine: str,
) -> list[str]:
    arguments = [
        "build", "--polymer", "chitin", "--allomorph", allomorph,
        "--shape", shape,
    ]
    if shape == "bulk":
        cells = ("6", "2", "3") if allomorph == "alpha" else ("5", "3", "3")
        arguments.extend(("--cells", *cells))
    else:
        arguments.extend(
            (
                "--face", "1", "0", "0", "--layers", "2",
                "--chains-per-layer", "4", "--degree-polymerization", "6",
                "--vacuum", "30",
            )
        )
        if boundary == "finite":
            arguments.extend(("--padding-longitudinal", "12"))
    arguments.extend(
        (
            "--chain-boundary", boundary,
            "--force-field", "glycam",
            "--md-engine", engine,
            "--tleap-executable", tleap,
            "--parmed-python", parmed_python,
            "--output", str(output),
        )
    )
    return arguments


def _explicit_pair_count(topology: Path) -> int:
    section = ""
    count = 0
    for line in topology.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped[1:-1].strip().lower()
        elif section == "pairs" and stripped and not stripped.startswith(";"):
            fields = stripped.split()
            count += len(fields) >= _MINIMUM_PAIR_FIELDS and fields[2] == "2"
    return count


def _gromacs_case(  # noqa: PLR0913, PLR0917
    directory: Path,
    allomorph: str,
    shape: str,
    boundary: str,
    tleap: str,
    parmed_python: str,
    gmx: str,
) -> dict[str, object]:
    if directory.exists():
        shutil.rmtree(directory)
    directory.mkdir(parents=True)
    assert main(
        _build_arguments(
            allomorph, shape, boundary, directory, tleap, parmed_python, "gromacs"
        )
    ) == 0
    report = json.loads((directory / "build_report.json").read_text(encoding="utf-8"))
    logs = []
    logs.append(
        _run(
            [
                gmx, "grompp", "-f", "minim.mdp", "-p", "system.top",
                "-c", "structure.gro", "-o", "em.tpr", "-po", "em.mdout.mdp",
                "-maxwarn", "0",
            ],
            directory,
        )
    )
    logs.append(_run([gmx, "mdrun", "-deffnm", "em", "-nt", "1"], directory))
    logs.append(
        _run(
            [
                gmx, "grompp", "-f", "nvt.mdp", "-p", "system.top",
                "-c", "em.gro", "-o", "nvt.tpr", "-po", "nvt.mdout.mdp",
                "-maxwarn", "0",
            ],
            directory,
        )
    )
    logs.append(_run([gmx, "mdrun", "-deffnm", "nvt", "-nt", "1"], directory))
    log_text = "\n".join(logs) + (directory / "nvt.log").read_text(encoding="utf-8")
    if _PROBLEM.search(log_text):
        raise RuntimeError(f"problem marker found in {directory.name} GROMACS output")
    if not re.search(r"\b1000\s+1\.00000\b", (directory / "nvt.log").read_text()):
        raise RuntimeError(f"{directory.name} did not complete the one-ps NVT run")
    return {
        "allomorph": allomorph,
        "shape": shape,
        "chain_boundary": boundary,
        "atom_count": report["force_field_atom_count"],
        "chain_count": report["chain_count"],
        "explicit_mixed_scaling_pairs_per_chain": _explicit_pair_count(
            directory / "system.top"
        ),
        "grompp_maxwarn": 0,
        "minimization_completed": True,
        "nvt_completed_ps": 1.0,
        "periodic_bond_validation": report.get("periodic_glycosidic_bond_validation"),
        "status": "passed",
    }


def _amber_case(
    directory: Path,
    tleap: str,
    parmed_python: str,
    sander: str,
) -> dict[str, object]:
    if directory.exists():
        shutil.rmtree(directory)
    directory.mkdir(parents=True)
    assert main(
        _build_arguments(
            "beta", "slab", "finite", directory, tleap, parmed_python, "amber"
        )
    ) == 0
    _run(
        [
            sander, "-O", "-i", "min.in", "-p", "structure.prmtop",
            "-c", "structure.inpcrd", "-o", "min.out", "-r", "min.rst7",
        ],
        directory,
    )
    _run(
        [
            sander, "-O", "-i", "nvt.in", "-p", "structure.prmtop",
            "-c", "min.rst7", "-o", "nvt.out", "-r", "nvt.rst7",
        ],
        directory,
    )
    text = (directory / "nvt.out").read_text(encoding="utf-8")
    if "NSTEP =     1000" not in text or _PROBLEM.search(text):
        raise RuntimeError("Amber finite beta-chitin smoke test did not complete cleanly")
    report = json.loads((directory / "build_report.json").read_text(encoding="utf-8"))
    return {
        "allomorph": "beta",
        "shape": "slab",
        "chain_boundary": "finite",
        "atom_count": report["force_field_atom_count"],
        "chain_count": report["chain_count"],
        "minimization_completed": True,
        "nvt_completed_ps": 1.0,
        "status": "passed",
    }


def _amber_gromacs_energy_comparison(
    directory: Path,
    tleap: str,
    parmed_python: str,
    sander: str,
    gmx_double: str,
) -> dict[str, object]:
    """Compare finite chitin bonded and 1-4 terms on identical coordinates."""
    if directory.exists():
        shutil.rmtree(directory)
    directory.mkdir(parents=True)
    arguments = _build_arguments(
        "beta", "slab", "finite", directory, tleap, parmed_python, "gromacs"
    )
    amber_index = arguments.index("gromacs")
    arguments[amber_index : amber_index + 1] = ("gromacs", "--md-engine", "amber")
    assert main(arguments) == 0
    coordinate_script = (
        "import parmed as p,sys; "
        "a=p.load_file(sys.argv[1]); g=p.load_file(sys.argv[2]); "
        "a.coordinates=g.coordinates; a.box=g.box; a.save(sys.argv[3],overwrite=True)"
    )
    _run(
        [
            parmed_python,
            "-c",
            coordinate_script,
            "structure.prmtop",
            "structure.gro",
            "identical.inpcrd",
        ],
        directory,
    )
    (directory / "amber_sp.in").write_text(
        """Single-point energy on GROMACS-rounded coordinates
&cntrl
  imin=1, maxcyc=0, ncyc=0, ntb=1, cut=10.0,
/
&ewald
  vdwmeth=0,
/
""",
        encoding="utf-8",
    )
    _run(
        [
            sander, "-O", "-i", "amber_sp.in", "-p", "structure.prmtop",
            "-c", "identical.inpcrd", "-o", "amber_sp.out", "-r", "amber_sp.rst7",
        ],
        directory,
    )
    (directory / "gromacs_sp.mdp").write_text(
        """integrator = md
nsteps = 0
constraints = none
cutoff-scheme = Verlet
nstlist = 10
rlist = 1.0
coulombtype = Cut-off
rcoulomb = 1.0
vdwtype = Cut-off
rvdw = 1.0
DispCorr = no
pbc = xyz
periodic-molecules = no
nstenergy = 1
""",
        encoding="utf-8",
    )
    _run(
        [
            gmx_double, "grompp", "-f", "gromacs_sp.mdp", "-p", "system.top",
            "-c", "structure.gro", "-o", "gromacs_sp.tpr", "-po", "gromacs_sp.mdout.mdp",
            "-maxwarn", "0",
        ],
        directory,
    )
    _run(
        [gmx_double, "mdrun", "-s", "gromacs_sp.tpr", "-deffnm", "gromacs_sp"],
        directory,
    )
    _run_with_input(
        [gmx_double, "energy", "-f", "gromacs_sp.edr", "-o", "gromacs_sp.xvg"],
        directory,
        "Bond\nAngle\nProper-Dih.\nPer.-Imp.-Dih.\nLJ-14\nCoulomb-14\n0\n",
    )
    xvg_lines = (directory / "gromacs_sp.xvg").read_text(encoding="utf-8").splitlines()
    legends = [
        match.group(1)
        for line in xvg_lines
        if (match := re.match(r'@ s\d+ legend "(.*)"', line))
    ]
    values = next(
        line.split()[1:]
        for line in reversed(xvg_lines)
        if line and not line.startswith(("#", "@"))
    )
    gromacs = {
        name: float(value) / _KJ_PER_KCAL
        for name, value in zip(legends, values, strict=True)
    }
    amber_text = (directory / "amber_sp.out").read_text(encoding="utf-8")
    amber_matches = {
        name: float(value)
        for name, value in re.findall(
            r"\b(BOND|ANGLE|DIHED|1-4 VDW|1-4 EEL)\s*=\s*([-+0-9.Ee]+)",
            amber_text,
        )
    }
    terms = {
        "bonds": (gromacs["Bond"], amber_matches["BOND"]),
        "angles": (gromacs["Angle"], amber_matches["ANGLE"]),
        "dihedrals": (
            gromacs["Proper Dih."] + gromacs["Per. Imp. Dih."],
            amber_matches["DIHED"],
        ),
        "1-4_vdw": (gromacs["LJ-14"], amber_matches["1-4 VDW"]),
        "1-4_electrostatics": (gromacs["Coulomb-14"], amber_matches["1-4 EEL"]),
    }
    coulomb_ratio = _GROMACS_COULOMB_KCAL_ANGSTROM / _AMBER_COULOMB_KCAL_ANGSTROM
    result_terms: dict[str, object] = {}
    for name, (gromacs_value, amber_value) in terms.items():
        corrected_gromacs = (
            gromacs_value / coulomb_ratio if name == "1-4_electrostatics" else gromacs_value
        )
        result_terms[name] = {
            "gromacs_kcal_mol": gromacs_value,
            "amber_kcal_mol": amber_value,
            "raw_difference_kcal_mol": gromacs_value - amber_value,
            "coulomb_constant_corrected_difference_kcal_mol": corrected_gromacs
            - amber_value,
        }
    maximum_corrected = max(
        abs(term["coulomb_constant_corrected_difference_kcal_mol"])
        for term in result_terms.values()
        if isinstance(term, dict)
    )
    if maximum_corrected > _ENERGY_DIFFERENCE_TOLERANCE_KCAL_MOL:
        raise RuntimeError(
            f"Amber/GROMACS bonded or 1-4 energy difference is too large: {maximum_corrected}"
        )
    return {
        "system": "finite beta-chitin (100) slab, 8 chains, DP6",
        "coordinates": "identical GROMACS-rounded coordinates",
        "gromacs_precision": "double",
        "coulomb_constants_kcal_angstrom_per_mol_e2": {
            "gromacs": _GROMACS_COULOMB_KCAL_ANGSTROM,
            "amber": _AMBER_COULOMB_KCAL_ANGSTROM,
            "ratio": coulomb_ratio,
        },
        "terms": result_terms,
        "maximum_corrected_absolute_difference_kcal_mol": maximum_corrected,
        "status": "passed",
    }


def main_cli() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--gmx-double", default=shutil.which("gmx_d"))
    parser.add_argument(
        "--output",
        type=Path,
        default=_ROOT / "validation" / "results" / "chitin_glycam_engines.json",
    )
    args = parser.parse_args()
    tleap = shutil.which("tleap") or "tleap"
    gmx = shutil.which("gmx") or "gmx"
    sander = shutil.which("sander") or "sander"
    parmed_python = find_parmed_python()
    cases = [
        (allomorph, shape, boundary)
        for allomorph in ("alpha", "beta")
        for shape, boundary in (
            ("bulk", "periodic"),
            ("slab", "periodic"),
            ("slab", "finite"),
        )
    ]
    gromacs = [
        _gromacs_case(
            args.work_directory / f"{allomorph}_{shape}_{boundary}",
            allomorph,
            shape,
            boundary,
            tleap,
            parmed_python,
            gmx,
        )
        for allomorph, shape, boundary in cases
    ]
    amber = _amber_case(
        args.work_directory / "beta_slab_finite_amber",
        tleap,
        parmed_python,
        sander,
    )
    if args.gmx_double is None:
        raise FileNotFoundError("double-precision gmx_d is required for the energy comparison")
    energy_comparison = _amber_gromacs_energy_comparison(
        args.work_directory / "beta_slab_finite_energy_comparison",
        tleap,
        parmed_python,
        sander,
        args.gmx_double,
    )
    result = {
        "provenance": provenance(("numpy",)),
        "engine_versions": _external_versions(
            gmx, parmed_python, args.work_directory / "beta_slab_finite_amber" / "min.out"
        ),
        "scope": "GLYCAM06j-1 alpha- and beta-chitin engine smoke tests",
        "force_field_accuracy_assessed": False,
        "gromacs_cases": gromacs,
        "amber_case": amber,
        "amber_gromacs_energy_comparison": energy_comparison,
        "status": "passed",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main_cli()
