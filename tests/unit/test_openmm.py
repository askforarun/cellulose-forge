"""Tests for OpenMM output and validation helpers."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from cellulose_forge.openmm import validate_openmm_system, write_openmm_runner


def test_write_openmm_runner_loads_generated_amber_files(tmp_path: Path) -> None:
    runner = write_openmm_runner(tmp_path)

    text = runner.read_text(encoding="utf-8")
    assert "AmberPrmtopFile(\"structure.prmtop\")" in text
    assert "AmberInpcrdFile(\"structure.inpcrd\")" in text
    assert "nonbondedMethod=PME" in text
    assert "simulation.minimizeEnergy" in text


def test_validate_openmm_system_returns_engine_metadata(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    topology = tmp_path / "structure.prmtop"
    coordinates = tmp_path / "structure.inpcrd"
    topology.touch()
    coordinates.touch()

    def run(*args: object, **kwargs: object) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout=(
                '{"openmm_version":"8.5","platform":"CPU","particle_count":42,'
                '"force_count":6,"initial_potential_energy_kj_mol":100.0,'
                '"minimized_potential_energy_kj_mol":80.0,'
                '"minimization_iterations_requested":10,"status":"passed"}\n'
            ),
            stderr="",
        )

    monkeypatch.setattr("cellulose_forge.openmm.subprocess.run", run)
    metadata = validate_openmm_system(
        topology, coordinates, python_executable=sys.executable
    )

    assert metadata["status"] == "passed"
    assert metadata["platform"] == "CPU"
    assert metadata["particle_count"] == 42
    assert metadata["python_executable"]


def test_validate_openmm_system_rejects_nonfinite_energy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    topology = tmp_path / "structure.prmtop"
    coordinates = tmp_path / "structure.inpcrd"
    topology.touch()
    coordinates.touch()

    def run(*args: object, **kwargs: object) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout=(
                '{"initial_potential_energy_kj_mol":NaN,'
                '"minimized_potential_energy_kj_mol":80.0}\n'
            ),
            stderr="",
        )

    monkeypatch.setattr("cellulose_forge.openmm.subprocess.run", run)
    with pytest.raises(ValueError, match="non-finite"):
        validate_openmm_system(
            topology, coordinates, python_executable=sys.executable
        )
