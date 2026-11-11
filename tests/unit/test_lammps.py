"""Tests for native LAMMPS output."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from cellulose_forge.glycam_mapping import GlycamMappedSystem
from cellulose_forge.lammps import write_lammps_system


def _system() -> GlycamMappedSystem:
    return GlycamMappedSystem(
        atom_names=("C1", "O4", "C1", "O4"),
        residue_names=("4GB",) * 4,
        residue_numbers=(1, 1, 2, 2),
        chain_numbers=(1, 1, 2, 2),
        cartesian_coordinates=np.zeros((4, 3)),
        basis_vectors=np.diag([20.0, 20.0, 20.0]),
    )


def test_writer_records_molecules_and_emits_safe_input(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    topology = tmp_path / "system.top"
    coordinates = tmp_path / "structure.gro"
    topology.touch()
    coordinates.touch()
    captured: dict[str, object] = {}

    def run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        request = json.loads(Path(command[-1]).read_text(encoding="utf-8"))
        captured.update(request)
        Path(str(request["data_path"])).write_text("LAMMPS data\n", encoding="utf-8")
        metadata = {
            "atom_count": 4,
            "bond_count": 2,
            "angle_count": 0,
            "dihedral_count": 0,
            "special_bonds": "lj/coul 0 0 1",
        }
        return subprocess.CompletedProcess(
            args=command, returncode=0, stdout=json.dumps(metadata), stderr=""
        )

    monkeypatch.setattr("cellulose_forge.lammps.subprocess.run", run)
    files, metadata = write_lammps_system(
        _system(), topology, coordinates, tmp_path / "out", python_executable=sys.executable
    )

    assert captured["molecule_ids"] == [1, 1, 2, 2]
    assert files.data.is_file()
    text = files.input.read_text(encoding="utf-8")
    assert "special_bonds lj/coul 0.0 0.0 1.0" in text
    assert "dihedral_style charmm" in text
    assert "dihedral_style fourier" not in text
    assert "read_data system.data" in text
    assert "pair_modify mix arithmetic tail no" in text
    assert "not a dynamics protocol" in text
    assert metadata["atom_count"] == 4


def test_writer_surfaces_converter_failure(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(
            args=command, returncode=1, stdout="", stderr="unsupported exception"
        )

    monkeypatch.setattr("cellulose_forge.lammps.subprocess.run", run)
    with pytest.raises(RuntimeError, match="unsupported exception"):
        write_lammps_system(_system(), tmp_path / "a.top", tmp_path / "a.gro", tmp_path)


def test_charmm36_writer_uses_exact_lammps_styles(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    topology = tmp_path / "system.top"
    coordinates = tmp_path / "structure.gro"
    topology.touch()
    coordinates.touch()

    def run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        request = json.loads(Path(command[-1]).read_text(encoding="utf-8"))
        assert request["force_field"] == "charmm36"
        Path(str(request["data_path"])).write_text("LAMMPS data\n", encoding="utf-8")
        return subprocess.CompletedProcess(
            args=command,
            returncode=0,
            stdout=json.dumps(
                {
                    "atom_count": 4,
                    "bond_count": 2,
                    "angle_count": 1,
                    "dihedral_count": 1,
                    "special_bonds": "charmm",
                }
            ),
            stderr="",
        )

    monkeypatch.setattr("cellulose_forge.lammps.subprocess.run", run)
    files, _ = write_lammps_system(
        _system(),
        topology,
        coordinates,
        tmp_path / "out",
        force_field="charmm36",
    )
    text = files.input.read_text(encoding="utf-8")
    assert "angle_style charmm" in text
    assert "dihedral_style charmmfsw" in text
    assert "pair_style lj/charmmfsw/coul/long 10.0 12.0" in text
    assert "special_bonds charmm" in text


def test_glycam_bulk_enables_lennard_jones_tail_correction(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    topology = tmp_path / "system.top"
    coordinates = tmp_path / "structure.gro"
    topology.touch()
    coordinates.touch()

    def run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        request = json.loads(Path(command[-1]).read_text(encoding="utf-8"))
        Path(str(request["data_path"])).write_text("LAMMPS data\n", encoding="utf-8")
        return subprocess.CompletedProcess(
            args=command,
            returncode=0,
            stdout=json.dumps({"atom_count": 4, "bond_count": 2}),
            stderr="",
        )

    monkeypatch.setattr("cellulose_forge.lammps.subprocess.run", run)
    files, metadata = write_lammps_system(
        _system(),
        topology,
        coordinates,
        tmp_path / "out",
        dispersion_correction=True,
    )
    assert "pair_modify mix arithmetic tail yes" in files.input.read_text(encoding="utf-8")
    assert metadata["dispersion_correction"] is True
