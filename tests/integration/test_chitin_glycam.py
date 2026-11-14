"""Real-tool validation of GLYCAM06j-1 alpha- and beta-chitin output."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from cellulose_forge.cli import main
from cellulose_forge.glycam import find_parmed_python


def _tools() -> tuple[str, str, str]:
    tleap = shutil.which("tleap")
    gmx = shutil.which("gmx")
    if tleap is None or gmx is None:
        pytest.skip("tleap and gmx are required for chitin GLYCAM integration tests")
    try:
        parmed_python = find_parmed_python()
    except FileNotFoundError:
        pytest.skip("a Python environment containing ParmEd is required")
    return tleap, gmx, parmed_python


@pytest.mark.parametrize("allomorph", ["alpha", "beta"])
@pytest.mark.parametrize("boundary", ["finite", "periodic"])
def test_chitin_glycam_slab_preprocesses_without_warnings(
    tmp_path: Path,
    allomorph: str,
    boundary: str,
) -> None:
    tleap, gmx, parmed_python = _tools()
    output = tmp_path / f"{allomorph}_{boundary}"
    arguments = [
        "build",
        "--polymer", "chitin",
        "--allomorph", allomorph,
        "--shape", "slab",
        "--face", "1", "0", "0",
        "--layers", "2",
        "--chains-per-layer", "4",
        "--degree-polymerization", "6",
        "--chain-boundary", boundary,
        "--vacuum", "30",
        "--force-field", "glycam",
        "--md-engine", "gromacs",
        "--tleap-executable", tleap,
        "--parmed-python", parmed_python,
        "--output", str(output),
    ]
    if boundary == "finite":
        arguments.extend(("--padding-longitudinal", "12"))

    assert main(arguments) == 0
    topology = (output / "system.top").read_text(encoding="utf-8")
    assert "CHITIN_CHAIN  8" in topology
    assert topology.count(" 2 0.8333333333 ") == 6
    completed = subprocess.run(
        [
            gmx,
            "grompp",
            "-f", "minim.mdp",
            "-p", "system.top",
            "-c", "structure.gro",
            "-o", "minim.tpr",
            "-po", "minim.mdout.mdp",
            "-maxwarn", "0",
        ],
        cwd=output,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "WARNING" not in completed.stdout + completed.stderr


def _openmm_python() -> str:
    candidates = [sys.executable]
    configured = os.environ.get("CELLULOSE_FORGE_OPENMM_PYTHON")
    if configured:
        candidates.insert(0, configured)
    for candidate in candidates:
        result = subprocess.run(
            [candidate, "-c", "import openmm"], capture_output=True, check=False
        )
        if result.returncode == 0:
            return candidate
    pytest.skip("a Python interpreter containing OpenMM is required")


@pytest.mark.parametrize("boundary", ["finite", "periodic"])
def test_chitin_glycam_builds_openmm_and_lammps(tmp_path: Path, boundary: str) -> None:
    tleap, _, parmed_python = _tools()
    openmm_python = _openmm_python()
    output = tmp_path / boundary
    arguments = [
        "build", "--polymer", "chitin", "--allomorph", "beta",
        "--shape", "slab", "--face", "1", "0", "0", "--layers", "2",
        "--chains-per-layer", "4", "--degree-polymerization", "6",
        "--chain-boundary", boundary, "--vacuum", "30", "--force-field", "glycam",
        "--md-engine", "openmm", "--md-engine", "lammps",
        "--openmm-python", openmm_python, "--tleap-executable", tleap,
        "--parmed-python", parmed_python, "--output", str(output),
    ]
    if boundary == "finite":
        arguments.extend(("--padding-longitudinal", "12"))
    assert main(arguments) == 0
    assert (output / "system.xml").is_file()
    assert (output / "in.chitin").is_file()
    assert (output / "explicit_pairs.table").is_file()
    report = json.loads((output / "build_report.json").read_text(encoding="utf-8"))
    assert report["lammps_validation"]["explicit_pair_override_count"] == 48


def test_chitin_glycam_restraints_use_chitin_define(tmp_path: Path) -> None:
    tleap, gmx, parmed_python = _tools()
    output = tmp_path / "restrained"
    assert main(
        [
            "build",
            "--polymer", "chitin",
            "--allomorph", "beta",
            "--shape", "slab",
            "--face", "1", "0", "0",
            "--layers", "2",
            "--chains-per-layer", "4",
            "--degree-polymerization", "6",
            "--chain-boundary", "periodic",
            "--vacuum", "30",
            "--force-field", "glycam",
            "--md-engine", "gromacs",
            "--restrain-layers", "1",
            "--restraint-force", "1000", "1000", "1000",
            "--tleap-executable", tleap,
            "--parmed-python", parmed_python,
            "--output", str(output),
        ]
    ) == 0
    topology = (output / "system.top").read_text(encoding="utf-8")
    assert "CHITIN_CHAIN_RESTRAINED  4" in topology
    assert "#ifdef POSRES_CHITIN" in topology
    mdp = (output / "minim.mdp").read_text(encoding="utf-8")
    (output / "restrained.mdp").write_text(
        mdp.replace("; define                  = -DPOSRES_CHITIN", "define = -DPOSRES_CHITIN"),
        encoding="utf-8",
    )
    completed = subprocess.run(
        [
            gmx, "grompp", "-f", "restrained.mdp", "-p", "system.top",
            "-c", "structure.gro", "-r", "structure.gro", "-o", "restrained.tpr",
            "-po", "restrained.mdout.mdp", "-maxwarn", "0",
        ],
        cwd=output,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
