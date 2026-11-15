"""Real-topology integration checks for experimental NAMD output."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from cellulose_forge.cli import main
from cellulose_forge.forcefields import discover_charmm36_directory


def _force_field() -> Path:
    try:
        return discover_charmm36_directory()
    except FileNotFoundError:
        pytest.skip("the pinned CHARMM36 release is not installed")


def _section_rows(text: str, section: str) -> list[list[str]]:
    active = False
    rows: list[list[str]] = []
    for raw in text.splitlines():
        stripped = raw.split(";", maxsplit=1)[0].strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            active = stripped[1:-1].strip() == section
            continue
        if active and stripped and not stripped.startswith("#"):
            rows.append(stripped.split())
    return rows


def _psf_section(path: Path, marker: str, arity: int) -> list[tuple[int, ...]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = next(index for index, line in enumerate(lines) if marker in line)
    count = int(lines[header].split()[0])
    values: list[int] = []
    row = header + 1
    while len(values) < count * arity:
        values.extend(int(value) for value in lines[row].split())
        row += 1
    return [
        tuple(value - 1 for value in values[start : start + arity])
        for start in range(0, len(values), arity)
    ]


def test_real_periodic_cellulose_psf_matches_gromacs_topology(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    force_field = _force_field()
    parameters = tmp_path / "par_all36_carb.prm"
    parameters.write_text("* checksum-controlled test placeholder\n", encoding="utf-8")
    checksum = hashlib.sha256(parameters.read_bytes()).hexdigest()
    monkeypatch.setattr("cellulose_forge.namd.CHARMM36_CARB_PARAMETER_SHA256", checksum)
    output = tmp_path / "build"

    assert (
        main(
            [
                "build", "--polymer", "cellulose", "--allomorph", "I-beta",
                "--shape", "slab", "--face", "1", "1", "0", "--layers", "2",
                "--chains-per-layer", "8", "--degree-polymerization", "6",
                "--chain-boundary", "periodic", "--vacuum", "30",
                "--force-field", "charmm36", "--md-engine", "gromacs",
                "--md-engine", "namd", "--charmm-force-field-directory", str(force_field),
                "--charmm-parameter-file", str(parameters), "--output", str(output),
            ]
        )
        == 0
    )  # fmt: skip

    report = json.loads((output / "build_report.json").read_text(encoding="utf-8"))
    chain_count = report["chain_count"]
    expected = report["interaction_counts_per_chain"]
    psf = output / "system.psf"
    assert len(_psf_section(psf, "!NBOND", 2)) == chain_count * expected["bonds"]
    assert len(_psf_section(psf, "!NTHETA", 3)) == chain_count * expected["angles"]
    proper = _psf_section(psf, "!NPHI", 4)
    improper = _psf_section(psf, "!NIMPHI", 4)
    assert len(proper) + len(improper) == chain_count * expected["dihedrals"]

    itp_atoms = _section_rows((output / "cellulose_charmm36.itp").read_text(), "atoms")
    psf_atom_lines = (output / "system.psf").read_text().split("!NATOM\n", maxsplit=1)[1]
    psf_atom_lines = psf_atom_lines.split("\n\n", maxsplit=1)[0].splitlines()
    assert len(psf_atom_lines) == chain_count * len(itp_atoms)
    for itp, psf_row in zip(itp_atoms, psf_atom_lines[: len(itp_atoms)], strict=True):
        psf_fields = psf_row.split()
        assert psf_fields[4] == itp[4]
        assert psf_fields[5] == itp[1]
        assert float(psf_fields[6]) == pytest.approx(float(itp[6]), abs=5e-7)
        assert float(psf_fields[7]) == pytest.approx(float(itp[7]), abs=5e-5)

    assert report["namd_validation"]["parameter_sha256"] == checksum

    gmx = shutil.which("gmx")
    if gmx is None:
        pytest.skip("GROMACS is required for the single-point adapter integration check")
    subprocess.run(
        [
            gmx,
            "grompp",
            "-f",
            "single_point_gromacs.mdp",
            "-c",
            "single_point_coordinates.gro",
            "-p",
            "system.top",
            "-o",
            "single_point.tpr",
            "-maxwarn",
            "0",
        ],
        cwd=output,
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run(
        [gmx, "mdrun", "-s", "single_point.tpr", "-deffnm", "single_point", "-nsteps", "0"],
        cwd=output,
        check=True,
        capture_output=True,
        text=True,
    )
    selections = "\n".join(
        (
            "Bond",
            "U-B",
            "Proper-Dih.",
            "Improper-Dih.",
            "LJ-14",
            "Coulomb-14",
            "LJ-(SR)",
            "Coulomb-(SR)",
            "Coul.-recip.",
            "Potential",
            "0",
        )
    )
    subprocess.run(
        [gmx, "energy", "-f", "single_point.edr", "-o", "single_point.xvg"],
        cwd=output,
        input=selections + "\n",
        check=True,
        capture_output=True,
        text=True,
    )
    xvg_lines = (output / "single_point.xvg").read_text(encoding="utf-8").splitlines()
    legends = [line.split('"')[1] for line in xvg_lines if line.startswith("@ s")]
    values = [
        float(value)
        for line in xvg_lines
        if line and not line.startswith(("#", "@"))
        for value in line.split()[1:]
    ]
    energies = dict(zip(legends, values, strict=True))
    grouped = {
        "bond": energies["Bond"],
        "angle": energies["U-B"],
        "dihedral": energies["Proper Dih."],
        "improper": energies.get("Improper Dih.", 0.0),
        "electrostatics": (
            energies["Coulomb-14"] + energies["Coulomb (SR)"] + energies["Coul. recip."]
        ),
        "vdw": energies["LJ-14"] + energies["LJ (SR)"],
        "potential": energies["Potential"],
    }
    kcal = {name: value / 4.184 for name, value in grouped.items()}
    namd_fields = [
        0.0,
        kcal["bond"],
        kcal["angle"],
        kcal["dihedral"],
        kcal["improper"],
        kcal["electrostatics"],
        kcal["vdw"],
        0.0,
        0.0,
        0.0,
        kcal["potential"],
        0.0,
        kcal["potential"],
    ]
    (output / "single_point_namd.log").write_text(
        "ETITLE: TS BOND ANGLE DIHED IMPRP ELECT VDW BOUNDARY MISC KINETIC "
        "TOTAL TEMP POTENTIAL\nENERGY: "
        + " ".join(f"{value:.12g}" for value in namd_fields)
        + "\n",
        encoding="utf-8",
    )
    subprocess.run(
        [
            sys.executable,
            "single_point_energy.py",
            "single_point_namd.log",
            "single_point.xvg",
            "single_point_coordinates.gro",
        ],
        cwd=output,
        check=True,
    )
    comparison = json.loads(
        (output / "single_point_energy_comparison.json").read_text(encoding="utf-8")
    )
    assert comparison["passed"] is True
