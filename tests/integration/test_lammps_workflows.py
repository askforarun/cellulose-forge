"""Cross-engine validation of experimental LAMMPS output."""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any

import pytest

from cellulose_forge.cli import main
from cellulose_forge.forcefields import discover_charmm36_directory
from cellulose_forge.glycam import GlycamChainFiles

pytestmark = [
    pytest.mark.lammps,
    pytest.mark.openmm,
    pytest.mark.gromacs,
    pytest.mark.scientific_integration,
]

_BUNDLE_ENVIRONMENT = "CELLULOSE_FORGE_LAMMPS_BUNDLE_DIR"


def _redact_value(value: Any, *, key: str = "") -> Any:
    """Remove machine-local paths from retained JSON metadata."""
    if isinstance(value, dict):
        return {name: _redact_value(item, key=name) for name, item in value.items()}
    if isinstance(value, list):
        return [_redact_value(item, key=key) for item in value]
    if isinstance(value, str) and value.startswith("/"):
        if key == "reference_force_field":
            return Path(value).name
        if key == "python_executable":
            return "python with OpenMM"
        return "<local-path>"
    return value


def _copy_redacted(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.suffix == ".json":
        payload = _redact_value(json.loads(source.read_text(encoding="utf-8")))
        destination.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        return
    text = source.read_text(encoding="utf-8", errors="replace")
    prefixes = {
        str(source.parent.resolve()),
        str(Path.home()),
        str(Path.cwd().resolve()),
    }
    for prefix in sorted(prefixes, key=len, reverse=True):
        text = text.replace(prefix, "<local-path>")
    text = re.sub(r"/(?:Users|private|tmp|home)/\S+", "<local-path>", text)
    destination.write_text(text, encoding="utf-8")


def _retain_validation_case(
    output: Path,
    case_name: str,
    comparison: dict[str, object],
    *,
    glycam: bool,
) -> None:
    """Retain only portable text evidence when explicitly requested."""
    root_value = os.environ.get(_BUNDLE_ENVIRONMENT)
    if not root_value:
        return
    case = Path(root_value) / case_name
    files = {
        "common": ("build_report.json",),
        "lammps": ("in.cellulose", "log.lammps"),
        "gromacs": (
            (
                "glycam_single_point.mdp",
                "glycam_single_point.log",
                "glycam_energy.xvg",
                "system.top",
            )
            if glycam
            else (
                "cellulose_charmm36.itp",
                "mdout.mdp",
                "single_point.mdp",
                "single_point.log",
                "single_point.xvg",
                "system.top",
            )
        ),
    }
    if glycam:
        files["openmm"] = ("run_openmm.py",)
    for engine, names in files.items():
        for name in names:
            _copy_redacted(output / name, case / engine / name)
    payload = {"generated_on": date.today().isoformat(), **comparison}
    (case / ".comparison.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _openmm_python() -> str:
    for candidate in (os.environ.get("CELLULOSE_FORGE_OPENMM_PYTHON"), sys.executable):
        if (
            candidate
            and subprocess.run(
                [candidate, "-c", "import openmm"], check=False, capture_output=True
            ).returncode
            == 0
        ):
            return candidate
    pytest.skip("no Python interpreter with OpenMM")


def _lammps() -> str:
    for name in (
        os.environ.get("CELLULOSE_FORGE_LAMMPS_EXECUTABLE"),
        "lmp",
        "lmp_serial",
        "lammps",
    ):
        if not name:
            continue
        if executable := shutil.which(name):
            return executable
    pytest.skip("no LAMMPS executable")


def _gmx() -> str:
    if executable := shutil.which("gmx"):
        return executable
    pytest.skip("no GROMACS executable")


def _lammps_energies(directory: Path) -> dict[str, float]:
    completed = subprocess.run(
        [_lammps(), "-in", "in.cellulose", "-log", "log.lammps"],
        cwd=directory,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    log = (directory / "log.lammps").read_text(encoding="utf-8")
    lines = log.splitlines()
    header_index = next(
        (
            index
            for index, line in enumerate(lines)
            if line.split()[:9]
            == [
                "Step",
                "Temp",
                "PotEng",
                "E_bond",
                "E_angle",
                "E_dihed",
                "E_vdwl",
                "E_coul",
                "E_long",
            ]
        ),
        None,
    )
    assert header_index is not None
    headings = lines[header_index].split()
    row = next(
        (line.split() for line in lines[header_index + 1 :] if line.split()[:1] == ["0"]),
        None,
    )
    assert row is not None
    values = dict(zip(headings, row, strict=True))
    names = {
        "temperature": "Temp",
        "potential": "PotEng",
        "bond": "E_bond",
        "angle": "E_angle",
        "dihedral": "E_dihed",
        "vdw": "E_vdwl",
        "coulomb": "E_coul",
        "long": "E_long",
    }
    return {
        name: float(values[column]) * (1.0 if name == "temperature" else 4.184)
        for name, column in names.items()
    }


def _gromacs_charmm36_energies(directory: Path, *, periodic: bool) -> dict[str, float]:
    mdp = directory / "single_point.mdp"
    mdp.write_text(
        "integrator = md\nnsteps = 0\nconstraints = none\ncutoff-scheme = Verlet\n"
        "coulombtype = PME\nrcoulomb = 1.2\n"
        "vdwtype = cutoff\nrvdw = 1.2\nvdw-modifier = force-switch\nrvdw-switch = 1.0\n"
        "coulomb-modifier = None\newald-rtol = 1e-5\nDispCorr = no\n"
        "pbc = xyz\nnstenergy = 1\n" + ("periodic-molecules = yes\n" if periodic else ""),
        encoding="utf-8",
    )
    for command in (
        [
            _gmx(),
            "grompp",
            "-f",
            mdp.name,
            "-c",
            "structure.gro",
            "-p",
            "system.top",
            "-o",
            "single_point.tpr",
        ],
        [
            _gmx(),
            "mdrun",
            "-s",
            "single_point.tpr",
            "-deffnm",
            "single_point",
            "-ntmpi",
            "1",
            "-ntomp",
            "1",
        ],
    ):
        completed = subprocess.run(
            command, cwd=directory, check=False, capture_output=True, text=True
        )
        assert completed.returncode == 0, completed.stdout + completed.stderr
    terms = (
        "Bond",
        "U-B",
        "Angle",
        "Proper Dih.",
        "Improper Dih.",
        "LJ-14",
        "Coulomb-14",
        "LJ (SR)",
        "Coulomb (SR)",
        "Coul. recip.",
        "Potential",
    )
    completed = subprocess.run(
        [_gmx(), "energy", "-f", "single_point.edr", "-o", "single_point.xvg"],
        cwd=directory,
        input="\n".join(term.replace(" ", "-") for term in terms) + "\n0\n",
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    values: dict[str, float] = {}
    for term in terms:
        match = re.search(
            rf"^{re.escape(term)}\s+([-+0-9.eE]+)",
            completed.stdout,
            flags=re.MULTILINE,
        )
        if match is not None:
            values[term] = float(match.group(1))
    return values


def _gromacs_glycam_energies(directory: Path) -> dict[str, float]:
    mdp = directory / "glycam_single_point.mdp"
    mdp.write_text(
        "integrator = md\nnsteps = 0\ncontinuation = yes\nconstraints = none\n"
        "cutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\n"
        "vdwtype = cutoff\nrvdw = 1.0\nvdw-modifier = None\n"
        "coulomb-modifier = None\newald-rtol = 1e-5\nDispCorr = no\n"
        "pbc = xyz\nperiodic-molecules = yes\nnstenergy = 1\n",
        encoding="utf-8",
    )
    for command in (
        [
            _gmx(),
            "grompp",
            "-f",
            mdp.name,
            "-c",
            "structure.gro",
            "-p",
            "system.top",
            "-o",
            "glycam_single_point.tpr",
        ],
        [
            _gmx(),
            "mdrun",
            "-s",
            "glycam_single_point.tpr",
            "-deffnm",
            "glycam_single_point",
            "-ntmpi",
            "1",
            "-ntomp",
            "1",
        ],
    ):
        completed = subprocess.run(
            command, cwd=directory, check=False, capture_output=True, text=True
        )
        assert completed.returncode == 0, completed.stdout + completed.stderr
    terms = (
        "Bond",
        "Angle",
        "Proper Dih.",
        "LJ-14",
        "Coulomb-14",
        "LJ (SR)",
        "Coulomb (SR)",
        "Coul. recip.",
        "Potential",
    )
    completed = subprocess.run(
        [_gmx(), "energy", "-f", "glycam_single_point.edr", "-o", "glycam_energy.xvg"],
        cwd=directory,
        input="\n".join(term.replace(" ", "-") for term in terms) + "\n0\n",
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    values: dict[str, float] = {}
    for term in terms:
        match = re.search(
            rf"^{re.escape(term)}\s+([-+0-9.eE]+)",
            completed.stdout,
            flags=re.MULTILINE,
        )
        if match is not None:
            values[term] = float(match.group(1))
    return values


@pytest.mark.parametrize("allomorph,face", [("I-beta", (1, 1, 0)), ("II", (1, 1, 0))])
def test_periodic_glycam_lammps_energy_matches_openmm_and_gromacs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    allomorph: str,
    face: tuple[int, int, int],
) -> None:
    def prepare(**kwargs: object) -> tuple[GlycamChainFiles, dict[str, object]]:
        directory = Path(str(kwargs["directory"]))
        files = GlycamChainFiles(directory, "cellulose_periodic_chain")
        source = Path(__file__).parents[2] / "examples/glycam/dp12_periodic_reference"
        shutil.copyfile(source / "cellulose_periodic_chain.pdb", files.pdb)
        shutil.copyfile(source / "cellulose_periodic_chain.top", files.gromacs_topology)
        return files, {"test_fixture": "dp12_periodic_reference"}

    monkeypatch.setattr("cellulose_forge.cli.generate_glycam_gromacs_reference", prepare)
    output = tmp_path / allomorph.lower().replace("-", "_")
    assert (
        main(
            [
                "build",
                "--allomorph",
                allomorph,
                "--shape",
                "slab",
                "--face",
                *(str(value) for value in face),
                "--layers",
                "2",
                "--chains-per-layer",
                "4",
                "--degree-polymerization",
                "12",
                "--chain-boundary",
                "periodic",
                "--vacuum",
                "30",
                "--force-field",
                "glycam",
                "--md-engine",
                "gromacs",
                "--md-engine",
                "openmm",
                "--md-engine",
                "lammps",
                "--openmm-python",
                _openmm_python(),
                "--output",
                str(output),
            ]
        )
        == 0
    )
    lammps = _lammps_energies(output)
    lammps_kj_mol = lammps["potential"]
    report = json.loads((output / "build_report.json").read_text(encoding="utf-8"))
    openmm_kj_mol = report["openmm_validation"]["initial_potential_energy_kj_mol"]
    assert lammps_kj_mol == pytest.approx(openmm_kj_mol, rel=5.0e-4)
    gromacs = _gromacs_glycam_energies(output)
    assert lammps["bond"] == pytest.approx(gromacs["Bond"], rel=2.0e-5)
    assert lammps["angle"] == pytest.approx(gromacs["Angle"], rel=2.0e-5)
    assert lammps["dihedral"] == pytest.approx(gromacs["Proper Dih."], rel=2.0e-5)
    assert lammps["vdw"] == pytest.approx(gromacs["LJ-14"] + gromacs["LJ (SR)"], abs=0.1)
    assert lammps["coulomb"] + lammps["long"] == pytest.approx(
        gromacs["Coulomb-14"] + gromacs["Coulomb (SR)"] + gromacs["Coul. recip."],
        rel=5.0e-4,
    )
    case_name = f"glycam_periodic_{allomorph.lower().replace('-', '')}"
    _retain_validation_case(
        output,
        case_name,
        {
            "reference_engine": "OpenMM",
            "openmm_potential_kj_mol": openmm_kj_mol,
            "lammps_potential_kj_mol": lammps_kj_mol,
            "relative_difference": abs(lammps_kj_mol - openmm_kj_mol) / abs(openmm_kj_mol),
        },
        glycam=True,
    )


@pytest.mark.parametrize(
    "allomorph,face,chain_boundary",
    [
        ("I-beta", (1, 1, 0), "finite"),
        ("I-beta", (1, 1, 0), "periodic"),
        ("II", (1, 1, 0), "periodic"),
    ],
)
def test_charmm36_lammps_energy_matches_gromacs(
    tmp_path: Path,
    allomorph: str,
    face: tuple[int, int, int],
    chain_boundary: str,
) -> None:
    try:
        force_field_directory = discover_charmm36_directory()
    except (FileNotFoundError, ValueError) as error:
        pytest.skip(str(error))
    output = tmp_path / f"{allomorph.lower().replace('-', '_')}_{chain_boundary}"
    padding = ["--padding-longitudinal", "12"] if chain_boundary == "finite" else []
    assert (
        main(
            [
                "build",
                "--allomorph",
                allomorph,
                "--shape",
                "slab",
                "--face",
                *(str(value) for value in face),
                "--layers",
                "2",
                "--chains-per-layer",
                "6",
                "--degree-polymerization",
                "6",
                "--chain-boundary",
                chain_boundary,
                "--vacuum",
                "30",
                *padding,
                "--force-field",
                "charmm36",
                "--md-engine",
                "gromacs",
                "--md-engine",
                "lammps",
                "--openmm-python",
                _openmm_python(),
                "--charmm-force-field-directory",
                str(force_field_directory),
                "--output",
                str(output),
            ]
        )
        == 0
    )
    lammps = _lammps_energies(output)
    gromacs = _gromacs_charmm36_energies(output, periodic=chain_boundary == "periodic")
    assert lammps["bond"] == pytest.approx(gromacs["Bond"], rel=2.0e-5)
    assert lammps["angle"] == pytest.approx(
        gromacs.get("U-B", 0.0) + gromacs.get("Angle", 0.0), rel=2.0e-5
    )
    assert lammps["dihedral"] == pytest.approx(
        gromacs.get("Proper Dih.", 0.0) + gromacs.get("Improper Dih.", 0.0),
        rel=2.0e-5,
    )
    assert lammps["potential"] == pytest.approx(gromacs["Potential"], rel=5.0e-4)
    report = json.loads((output / "build_report.json").read_text(encoding="utf-8"))
    assert report["lammps_validation"]["urey_bradley_count"] > 0
    assert report["lammps_validation"]["weighted_one_four_count"] > 0
    case_name = f"charmm36_{chain_boundary}_{allomorph.lower().replace('-', '')}"
    _retain_validation_case(
        output,
        case_name,
        {
            "reference_engine": "GROMACS",
            "gromacs_potential_kj_mol": gromacs["Potential"],
            "lammps_potential_kj_mol": lammps["potential"],
            "relative_difference": abs(lammps["potential"] - gromacs["Potential"])
            / abs(gromacs["Potential"]),
        },
        glycam=False,
    )
