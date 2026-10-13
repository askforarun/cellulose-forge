"""Executable GROMACS smoke tests for the supported force-field backends."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

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


def _use_committed_glycam_reference(monkeypatch: pytest.MonkeyPatch) -> None:
    def prepare(**kwargs: object) -> tuple[GlycamChainFiles, dict[str, object]]:
        directory = Path(str(kwargs["directory"]))
        files = GlycamChainFiles(directory, "cellulose_chain")
        source = Path(__file__).parents[2] / "examples/glycam/dp12_reference"
        shutil.copyfile(source / "cellulose_chain.pdb", files.pdb)
        shutil.copyfile(source / "cellulose_chain.top", files.gromacs_topology)
        return files, {"test_fixture": "dp12_reference"}

    monkeypatch.setattr("cellulose_forge.cli.generate_glycam_gromacs_reference", prepare)


def test_finite_glycam_slab_passes_grompp(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _gmx()
    _use_committed_glycam_reference(monkeypatch)
    output = tmp_path / "glycam_finite_slab"
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
                "2",
                "--chains-per-layer",
                "4",
                "--degree-polymerization",
                "12",
                "--chain-boundary",
                "finite",
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
    _grompp(output, Path(__file__).parents[2] / "validation/grompp_minim.mdp")


@pytest.mark.parametrize("chain_boundary", ["finite", "periodic"])
def test_charmm36_slab_passes_grompp(tmp_path: Path, chain_boundary: str) -> None:
    _gmx()
    try:
        force_field = discover_charmm36_directory()
    except FileNotFoundError:
        pytest.skip("the pinned CHARMM36 release is not installed")
    output = tmp_path / f"charmm36_{chain_boundary}_slab"
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
                "2",
                "--chains-per-layer",
                "4",
                "--degree-polymerization",
                "8",
                "--chain-boundary",
                chain_boundary,
                "--vacuum",
                "30",
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
