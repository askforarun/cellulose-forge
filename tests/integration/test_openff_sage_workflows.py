"""Optional integration tests for experimental production OpenFF output."""

from __future__ import annotations

import importlib
import json
import re
from importlib.util import find_spec
from pathlib import Path

import pytest

from cellulose_forge.cli import main
from cellulose_forge.crystallography import MillerIndex, load_allomorph
from cellulose_forge.geometry import build_bulk_crystal, build_chain_parallel_slab
from cellulose_forge.openff_sage import build_sage_system, write_sage_outputs

pytestmark = pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")


@pytest.mark.parametrize(("boundary", "padding"), (("finite", 12.0), ("periodic", 0.0)))
def test_openff_slab_outputs(tmp_path: Path, boundary: str, padding: float) -> None:
    crystal = build_chain_parallel_slab(
        load_allomorph("i_beta"),
        MillerIndex(1, 1, 0),
        layers=2,
        chains_per_layer=4,
        degree_polymerization=6,
        vacuum=30.0,
        padding_z=padding,
        chain_boundary=boundary,
    )
    build = build_sage_system(crystal, boundary)
    report = write_sage_outputs(build, tmp_path, ("gromacs", "openmm"), dispersion_correction=False)
    assert report["dispersion_correction"] is False
    assert report["gromacs_coordinates"] == "structure.gro"
    assert (tmp_path / "system.top").is_file()
    assert re.search(r"^DispCorr\s*=\s*no$", (tmp_path / "minim.mdp").read_text(), re.MULTILINE)
    openmm = importlib.import_module("openmm")
    system = openmm.XmlSerializer.deserialize((tmp_path / "system.xml").read_text())
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            assert force.getUseDispersionCorrection() is False
        elif isinstance(force, openmm.CustomNonbondedForce):
            assert force.getUseLongRangeCorrection() is False


def test_openff_periodic_bulk_outputs(tmp_path: Path) -> None:
    crystal = build_bulk_crystal(load_allomorph("i_beta"), (3, 3, 2))
    build = build_sage_system(crystal, "periodic")
    report = write_sage_outputs(build, tmp_path, ("gromacs", "openmm"), dispersion_correction=True)
    assert report["dispersion_correction"] is True
    minim = (tmp_path / "minim.mdp").read_text()
    assert re.search(r"^DispCorr\s*=\s*EnerPres$", minim, re.MULTILINE)
    assert (tmp_path / "system.xml").is_file()
    assert (tmp_path / "topology.cif").is_file()
    openmm = importlib.import_module("openmm")
    system = openmm.XmlSerializer.deserialize((tmp_path / "system.xml").read_text())
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            assert force.getUseDispersionCorrection() is True
        elif isinstance(force, openmm.CustomNonbondedForce):
            assert force.getUseLongRangeCorrection() is True


@pytest.mark.parametrize(
    "arguments",
    (
        ["--fibril-model", "36-chain", "--chain-boundary", "periodic"],
        [
            "--fibril-profile", "2,3,3,2", "--chain-boundary", "finite",
            "--padding-longitudinal", "12",
        ],
    ),
)
def test_openff_fibril_cli_outputs(tmp_path: Path, arguments: list[str]) -> None:
    output = tmp_path / "fibril"
    assert (
        main(
            [
                "build", "--allomorph", "I-beta", "--shape", "fibril", *arguments,
                "--degree-polymerization", "6", "--vacuum", "20",
                "--force-field", "openff", "--md-engine", "gromacs", "--md-engine", "openmm",
                "--output", str(output),
            ]
        )
        == 0
    )
    report = json.loads((output / "build_report.json").read_text())
    assert report["shape"] == "fibril"
    assert report["dispersion_correction"] is False
    assert report["stereochemistry"]["passed"] is True
    assert re.search(r"^DispCorr\s*=\s*no$", (output / "minim.mdp").read_text(), re.MULTILINE)
    openmm = importlib.import_module("openmm")
    system = openmm.XmlSerializer.deserialize((output / "system.xml").read_text())
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            assert force.getUseDispersionCorrection() is False
        elif isinstance(force, openmm.CustomNonbondedForce):
            assert force.getUseLongRangeCorrection() is False


def test_openff_openmm_only_build_writes_no_gromacs_files(tmp_path: Path) -> None:
    output = tmp_path / "openmm_only"
    assert (
        main(
            [
                "build", "--allomorph", "I-beta", "--shape", "bulk", "--cells", "3", "3", "2",
                "--chain-boundary", "periodic", "--force-field", "openff",
                "--md-engine", "openmm", "--output", str(output),
            ]
        )
        == 0
    )
    for name in ("structure.gro", "system.top", "minim.mdp", "nvt.mdp"):
        assert not (output / name).exists(), name
    for name in ("system.xml", "state.xml", "topology.cif", "run_openmm.py"):
        assert (output / name).is_file(), name
    report = json.loads((output / "build_report.json").read_text())
    assert not any(key.startswith("gromacs") for key in report["outputs"])
