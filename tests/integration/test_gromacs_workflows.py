"""Executable GROMACS smoke tests for the supported force-field backends."""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from cellulose_forge.cli import main
from cellulose_forge.forcefields import discover_charmm36_directory
from cellulose_forge.glycam import GlycamChainFiles

pytestmark = [pytest.mark.gromacs, pytest.mark.scientific_integration]


def _gmx() -> str:
    executable = shutil.which("gmx")
    if executable is None:
        pytest.skip("GROMACS executable 'gmx' is not installed")
    return executable


def _grompp(directory: Path, mdp: Path) -> None:
    completed = subprocess.run(
        [
            _gmx(),
            "grompp",
            "-f",
            str(mdp.resolve()),
            "-c",
            "structure.gro",
            "-p",
            "system.top",
            "-o",
            "smoke.tpr",
            "-po",
            "processed.mdp",
        ],
        cwd=directory,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert (directory / "smoke.tpr").is_file()


def _single_point_bond_energy_per_atom(directory: Path) -> float:
    """Return the initial bonded energy normalized by the mapped atom count."""
    mdp = Path(__file__).parents[2] / "validation/single_point_pme.mdp"
    _grompp(directory, mdp)
    completed = subprocess.run(
        [
            _gmx(),
            "mdrun",
            "-s",
            "smoke.tpr",
            "-deffnm",
            "single_point",
            "-ntmpi",
            "1",
            "-ntomp",
            "1",
            "-rdd",
            "1.2",
        ],
        cwd=directory,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    energy = subprocess.run(
        [_gmx(), "energy", "-f", "single_point.edr", "-o", "bond_energy.xvg"],
        cwd=directory,
        input="Bond\n0\n",
        check=False,
        capture_output=True,
        text=True,
    )
    assert energy.returncode == 0, energy.stdout + energy.stderr
    match = re.search(r"^Bond\s+([-+0-9.eE]+)", energy.stdout, flags=re.MULTILINE)
    assert match is not None, energy.stdout
    atom_count = int((directory / "structure.gro").read_text(encoding="utf-8").splitlines()[1])
    return float(match.group(1)) / atom_count


def _gro_coordinates_and_box(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read Cartesian coordinates and an orthorhombic box from a GRO file."""
    lines = path.read_text(encoding="utf-8").splitlines()
    atom_count = int(lines[1])
    coordinates = np.asarray(
        [
            [float(line[20:28]), float(line[28:36]), float(line[36:44])]
            for line in lines[2 : atom_count + 2]
        ]
    )
    box = np.asarray([float(value) for value in lines[atom_count + 2].split()[:3]])
    return coordinates, box


def _use_committed_glycam_reference(
    monkeypatch: pytest.MonkeyPatch, *, periodic: bool = False
) -> None:
    def prepare(**kwargs: object) -> tuple[GlycamChainFiles, dict[str, object]]:
        directory = Path(str(kwargs["directory"]))
        stem = "cellulose_periodic_chain" if periodic else "cellulose_chain"
        fixture = "dp12_periodic_reference" if periodic else "dp12_reference"
        files = GlycamChainFiles(directory, stem)
        source = Path(__file__).parents[2] / f"examples/glycam/{fixture}"
        shutil.copyfile(source / f"{stem}.pdb", files.pdb)
        shutil.copyfile(source / f"{stem}.top", files.gromacs_topology)
        return files, {"test_fixture": "dp12_reference"}

    monkeypatch.setattr("cellulose_forge.cli.generate_glycam_gromacs_reference", prepare)


@pytest.mark.parametrize(
    ("allomorph", "face"),
    [
        ("I-alpha", (0, 1, 0)),
        ("I-beta", (1, 1, 0)),
        ("II", (1, 1, 0)),
        ("III-I", (1, 1, 0)),
    ],
)
def test_finite_glycam_slab_passes_grompp(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    allomorph: str,
    face: tuple[int, int, int],
) -> None:
    _gmx()
    _use_committed_glycam_reference(monkeypatch)
    output = tmp_path / f"glycam_{allomorph}_finite_slab"
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
                "finite",
                "--vacuum",
                "30",
                "--padding-longitudinal",
                "12",
                "--force-field",
                "glycam",
                "--output",
                str(output),
            ]
        )
        == 0
    )
    _grompp(output, Path(__file__).parents[2] / "validation/grompp_minim.mdp")
    assert _single_point_bond_energy_per_atom(output) < 50.0


@pytest.mark.parametrize(
    ("allomorph", "face"),
    [
        ("I-alpha", (0, 1, 0)),
        ("I-beta", (1, 1, 0)),
        ("II", (1, 1, 0)),
        ("III-I", (1, 1, 0)),
    ],
)
@pytest.mark.parametrize("chain_boundary", ["finite", "periodic"])
def test_charmm36_slab_passes_grompp(
    tmp_path: Path,
    chain_boundary: str,
    allomorph: str,
    face: tuple[int, int, int],
) -> None:
    _gmx()
    try:
        force_field = discover_charmm36_directory()
    except FileNotFoundError:
        pytest.skip("the pinned CHARMM36 release is not installed")
    output = tmp_path / f"charmm36_{allomorph}_{chain_boundary}_slab"
    padding_arguments = (
        ["--padding-longitudinal", "12"] if chain_boundary == "finite" else []
    )
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
                "8",
                "--chain-boundary",
                chain_boundary,
                "--vacuum",
                "30",
                *padding_arguments,
                "--force-field",
                "charmm36",
                "--charmm-force-field-directory",
                str(force_field),
                "--output",
                str(output),
            ]
        )
        == 0
    )
    _grompp(output, Path(__file__).parents[2] / "validation/grompp_minim.mdp")
    assert _single_point_bond_energy_per_atom(output) < 50.0


def test_periodic_i_beta_glycam_slab_retains_core_geometry_during_short_nvt(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Check that corrected I-beta (110) lattice sites remain stable for 1 ps."""
    _gmx()
    _use_committed_glycam_reference(monkeypatch, periodic=True)
    output = tmp_path / "glycam_i_beta_periodic_nvt"
    assert (
        main(
            [
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
                "4",
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
                "--output",
                str(output),
            ]
        )
        == 0
    )
    root = Path(__file__).parents[2]
    _grompp(output, root / "validation/periodic_pme_minim.mdp")
    minimization = subprocess.run(
        [_gmx(), "mdrun", "-s", "smoke.tpr", "-deffnm", "minim", "-ntmpi", "1", "-ntomp", "1"],
        cwd=output,
        check=False,
        capture_output=True,
        text=True,
    )
    assert minimization.returncode == 0, minimization.stdout + minimization.stderr
    preprocessing = subprocess.run(
        [
            _gmx(),
            "grompp",
            "-f",
            str((root / "validation/periodic_pme_nvt_1ps.mdp").resolve()),
            "-c",
            "minim.gro",
            "-p",
            "system.top",
            "-o",
            "nvt.tpr",
        ],
        cwd=output,
        check=False,
        capture_output=True,
        text=True,
    )
    assert preprocessing.returncode == 0, preprocessing.stdout + preprocessing.stderr
    dynamics = subprocess.run(
        [_gmx(), "mdrun", "-s", "nvt.tpr", "-deffnm", "nvt", "-ntmpi", "1", "-ntomp", "1"],
        cwd=output,
        check=False,
        capture_output=True,
        text=True,
    )
    assert dynamics.returncode == 0, dynamics.stdout + dynamics.stderr
    assert "LINCS WARNING" not in (dynamics.stdout + dynamics.stderr)

    reference, box = _gro_coordinates_and_box(output / "minim.gro")
    final, _ = _gro_coordinates_and_box(output / "nvt.gro")
    atoms_per_chain = 12 * 21
    core = slice(4 * atoms_per_chain, 12 * atoms_per_chain)
    displacement = final[core] - reference[core]
    displacement -= np.round(displacement / box) * box
    displacement -= np.mean(displacement, axis=0)
    core_rmsd_nm = float(np.sqrt(np.mean(np.sum(displacement**2, axis=1))))
    assert core_rmsd_nm < 0.1
