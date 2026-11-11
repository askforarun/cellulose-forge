"""Engine-level OpenMM tests compared against GROMACS single-point energies.

These tests run only when an interpreter that can import OpenMM is available,
either the test interpreter itself or ``CELLULOSE_FORGE_OPENMM_PYTHON``, and
when GROMACS is installed to provide the reference energies.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from cellulose_forge.cli import main
from cellulose_forge.forcefields import discover_charmm36_directory
from cellulose_forge.glycam import GlycamChainFiles

pytestmark = [pytest.mark.openmm, pytest.mark.gromacs, pytest.mark.scientific_integration]

_GROMACS_BONDED_TERMS = (
    "Bond",
    "U-B",
    "Angle",
    "Proper Dih.",
    "Improper Dih.",
    "Per. Imp. Dih.",
    "CMAP Dih.",
)
_OPENMM_BONDED_FORCES = (
    "HarmonicBondForce",
    "HarmonicAngleForce",
    "PeriodicTorsionForce",
    "CustomTorsionForce",
    "CMAPTorsionForce",
)
# Both engines evaluate the same GRO coordinates, so agreement is limited only
# by single-precision GROMACS output.
_RELATIVE_TOLERANCE = 2.0e-5


def _openmm_python() -> str:
    candidates = [os.environ.get("CELLULOSE_FORGE_OPENMM_PYTHON"), sys.executable]
    for candidate in candidates:
        if not candidate:
            continue
        probe = subprocess.run(
            [candidate, "-c", "import openmm"],
            check=False,
            capture_output=True,
        )
        if probe.returncode == 0:
            return candidate
    pytest.skip("no Python interpreter with OpenMM; set CELLULOSE_FORGE_OPENMM_PYTHON")


def _gmx() -> str:
    executable = shutil.which("gmx")
    if executable is None:
        pytest.skip("GROMACS executable 'gmx' is not installed")
    return executable


def _use_committed_glycam_reference(
    monkeypatch: pytest.MonkeyPatch, *, periodic: bool
) -> None:
    def prepare(**kwargs: object) -> tuple[GlycamChainFiles, dict[str, object]]:
        directory = Path(str(kwargs["directory"]))
        stem = "cellulose_periodic_chain" if periodic else "cellulose_chain"
        fixture = "dp12_periodic_reference" if periodic else "dp12_reference"
        files = GlycamChainFiles(directory, stem)
        source = Path(__file__).parents[2] / f"examples/glycam/{fixture}"
        shutil.copyfile(source / f"{stem}.pdb", files.pdb)
        shutil.copyfile(source / f"{stem}.top", files.gromacs_topology)
        return files, {"test_fixture": fixture}

    monkeypatch.setattr("cellulose_forge.cli.generate_glycam_gromacs_reference", prepare)


def _gromacs_energies(directory: Path, *, periodic: bool) -> dict[str, float]:
    """Return GROMACS zero-step bonded and 1-4 energies for a written system."""
    gmx = _gmx()
    mdp = directory / "single_point.mdp"
    mdp.write_text(
        "integrator = md\nnsteps = 0\nconstraints = none\ncutoff-scheme = Verlet\n"
        "coulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\nnstenergy = 1\n"
        + ("periodic-molecules = yes\n" if periodic else ""),
        encoding="utf-8",
    )
    for command in (
        [gmx, "grompp", "-f", mdp.name, "-c", "structure.gro", "-p", "system.top",
         "-o", "single_point.tpr"],
        [gmx, "mdrun", "-s", "single_point.tpr", "-deffnm", "single_point",
         "-ntmpi", "1", "-ntomp", "1"],
    ):
        completed = subprocess.run(
            command, cwd=directory, check=False, capture_output=True, text=True
        )
        assert completed.returncode == 0, completed.stdout + completed.stderr
    terms = (*_GROMACS_BONDED_TERMS, "LJ-14", "Coulomb-14")
    energy = subprocess.run(
        [gmx, "energy", "-f", "single_point.edr", "-o", "single_point.xvg"],
        cwd=directory,
        input="\n".join(term.replace(" ", "-") for term in terms) + "\n0\n",
        check=False,
        capture_output=True,
        text=True,
    )
    assert energy.returncode == 0, energy.stdout + energy.stderr
    values: dict[str, float] = {}
    for term in terms:
        match = re.search(
            rf"^{re.escape(term)}\s+([-+0-9.eE]+)", energy.stdout, flags=re.MULTILINE
        )
        if match is not None:
            values[term] = float(match.group(1))
    return values


def _assert_openmm_matches_gromacs(output: Path, *, periodic: bool) -> None:
    report = json.loads((output / "build_report.json").read_text(encoding="utf-8"))
    openmm_energies = report["openmm_validation"]["initial_energies_kj_mol"]
    gromacs = _gromacs_energies(output, periodic=periodic)

    openmm_bonded = sum(openmm_energies.get(name, 0.0) for name in _OPENMM_BONDED_FORCES)
    gromacs_bonded = sum(gromacs.get(name, 0.0) for name in _GROMACS_BONDED_TERMS)
    openmm_onefour = openmm_energies["NonbondedForce exceptions (1-4)"]
    gromacs_onefour = gromacs["LJ-14"] + gromacs["Coulomb-14"]

    assert openmm_bonded == pytest.approx(gromacs_bonded, rel=_RELATIVE_TOLERANCE)
    assert openmm_onefour == pytest.approx(gromacs_onefour, rel=_RELATIVE_TOLERANCE)
    assert report["openmm_validation"]["maximum_bond_length_deviation_nm"] < 0.05
    if periodic:
        assert "NonbondedForce exceptions" in report["openmm_validation"]["periodic_forces"]
    for name in ("system.xml", "state.xml", "topology.cif", "run_openmm.py"):
        assert (output / name).is_file()


def _slab_arguments(chain_boundary: str) -> list[str]:
    padding = ["--padding-longitudinal", "12"] if chain_boundary == "finite" else []
    return [
        "build",
        "--allomorph",
        "I-beta",
        "--shape",
        "slab",
        "--face",
        "1",
        "1",
        "0",
        "--layers",
        "2",
        "--chains-per-layer",
        "4",
        "--degree-polymerization",
        "12",
        "--chain-boundary",
        chain_boundary,
        "--vacuum",
        "30",
        *padding,
        "--md-engine",
        "gromacs",
        "--md-engine",
        "openmm",
    ]


@pytest.mark.parametrize("chain_boundary", ["finite", "periodic"])
def test_glycam_openmm_energies_match_gromacs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, chain_boundary: str
) -> None:
    python = _openmm_python()
    _gmx()
    periodic = chain_boundary == "periodic"
    _use_committed_glycam_reference(monkeypatch, periodic=periodic)
    output = tmp_path / f"glycam_{chain_boundary}"

    assert (
        main(
            [
                *_slab_arguments(chain_boundary),
                "--force-field",
                "glycam",
                "--openmm-python",
                python,
                "--output",
                str(output),
            ]
        )
        == 0
    )
    _assert_openmm_matches_gromacs(output, periodic=periodic)


@pytest.mark.parametrize("chain_boundary", ["finite", "periodic"])
def test_charmm36_openmm_energies_match_gromacs(
    tmp_path: Path, chain_boundary: str
) -> None:
    python = _openmm_python()
    _gmx()
    try:
        discover_charmm36_directory()
    except FileNotFoundError:
        pytest.skip("the pinned CHARMM36 release is not installed")
    output = tmp_path / f"charmm36_{chain_boundary}"

    assert (
        main(
            [
                *_slab_arguments(chain_boundary),
                "--force-field",
                "charmm36",
                "--openmm-python",
                python,
                "--output",
                str(output),
            ]
        )
        == 0
    )
    _assert_openmm_matches_gromacs(output, periodic=chain_boundary == "periodic")
