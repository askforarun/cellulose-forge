"""Tests for serialized OpenMM system output."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from cellulose_forge.glycam_mapping import GlycamMappedSystem
from cellulose_forge.openmm import (
    POSITION_RESTRAINT_SCALE_PARAMETER,
    OpenMMPositionRestraints,
    chain_axis,
    heavy_atom_restraints,
    render_openmm_runner,
    write_openmm_system,
)


def _two_chain_system() -> GlycamMappedSystem:
    """Return two two-residue chains running along Cartesian y."""
    names = ("C1", "H1", "O4", "C1", "H1", "O4")
    coordinates = np.asarray(
        [[0.0, y, 0.0] for y in (0.0, 0.5, 1.0, 5.2, 5.7, 6.2)]
    )
    return GlycamMappedSystem(
        atom_names=names * 2,
        residue_names=("4GB",) * 12,
        residue_numbers=(1, 1, 1, 2, 2, 2) * 2,
        chain_numbers=(1,) * 6 + (2,) * 6,
        cartesian_coordinates=np.vstack(
            [coordinates, coordinates + np.asarray([8.0, 0.0, 0.0])]
        ),
        basis_vectors=np.diag([24.0, 25.0, 30.0]),
    )


def _metadata(**overrides: object) -> dict[str, object]:
    metadata: dict[str, object] = {
        "openmm_version": "8.5",
        "particle_count": 12,
        "periodic_forces": [
            "NonbondedForce exceptions",
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
        ],
        "initial_energies_kj_mol": {
            "HarmonicBondForce": 120.0,
            "HarmonicAngleForce": 60.0,
            "PeriodicTorsionForce": 30.0,
            "NonbondedForce": -400.0,
            "NonbondedForce exceptions (1-4)": 900.0,
        },
        "initial_potential_energy_kj_mol": -190.0,
        "minimized_potential_energy_kj_mol": -250.0,
        "maximum_bond_length_deviation_nm": 0.01,
    }
    metadata.update(overrides)
    return metadata


def _fake_builder(
    monkeypatch: pytest.MonkeyPatch, metadata: dict[str, object]
) -> dict[str, object]:
    captured: dict[str, object] = {}

    def run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        captured["request"] = json.loads(Path(command[-1]).read_text(encoding="utf-8"))
        return subprocess.CompletedProcess(
            args=command, returncode=0, stdout=json.dumps(metadata) + "\n", stderr=""
        )

    monkeypatch.setattr("cellulose_forge.openmm.subprocess.run", run)
    return captured


def test_chain_axis_follows_anomeric_carbons() -> None:
    assert chain_axis(_two_chain_system()) == "y"


def test_heavy_atom_restraints_select_only_requested_chain() -> None:
    restraints = heavy_atom_restraints(_two_chain_system(), (2,), (100.0, 200.0, 300.0))

    assert restraints.atom_indices == (6, 8, 9, 11)
    assert restraints.force_constants_kj_mol_nm2 == (100.0, 200.0, 300.0)


def test_runner_embeds_axis_periodicity_and_restraint_parameter() -> None:
    text = render_openmm_runner(chain_axis="y", periodic_chains=True, has_restraints=True)

    assert "CHAIN_AXIS = 'y'" in text
    assert "PERIODIC_CHAINS = True" in text
    assert "HAS_RESTRAINTS = True" in text
    assert POSITION_RESTRAINT_SCALE_PARAMETER in text
    assert "add --no-restraints for pressure coupling" in text
    assert 'XmlSerializer.deserialize' in text
    assert "__" not in text.replace("__name__", "").replace("__main__", "").replace(
        "__doc__", ""
    ).replace("__future__", "")


def test_write_openmm_system_sends_request_and_writes_runner(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    topology = tmp_path / "system.top"
    coordinates = tmp_path / "structure.gro"
    topology.touch()
    coordinates.touch()
    captured = _fake_builder(monkeypatch, _metadata())
    restraints = OpenMMPositionRestraints((0, 2), (1000.0, 1000.0, 1000.0))

    files, metadata = write_openmm_system(
        _two_chain_system(),
        topology,
        coordinates,
        tmp_path / "out",
        periodic_chains=True,
        dispersion_correction=False,
        restraints=restraints,
        python_executable=sys.executable,
    )

    request = captured["request"]
    assert isinstance(request, dict)
    assert request["periodic_chains"] is True
    assert request["expected_atom_count"] == 12
    assert request["dispersion_correction"] is False
    assert request["validation_ewald_error_tolerance"] < request["ewald_error_tolerance"]
    assert request["include_directory"] == str(tmp_path.resolve())
    assert request["position_restraints"] == {
        "atom_indices": [0, 2],
        "force_constants_kj_mol_nm2": [1000.0, 1000.0, 1000.0],
    }
    assert files.runner.is_file()
    assert "PERIODIC_CHAINS = True" in files.runner.read_text(encoding="utf-8")
    assert metadata["status"] == "passed"
    assert metadata["chain_axis"] == "y"
    assert metadata["position_restraint_atom_count"] == 2


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"minimized_potential_energy_kj_mol": float("nan")}, "non-finite"),
        ({"maximum_bond_length_deviation_nm": 1.04}, "periodic image"),
        (
            {"initial_energies_kj_mol": {"HarmonicBondForce": 1.0e6}},
            "per atom",
        ),
        ({"periodic_forces": ["HarmonicBondForce"]}, "exceptions"),
    ],
)
def test_write_openmm_system_rejects_invalid_systems(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    overrides: dict[str, object],
    message: str,
) -> None:
    topology = tmp_path / "system.top"
    topology.touch()
    _fake_builder(monkeypatch, _metadata(**overrides))

    with pytest.raises(ValueError, match=message):
        write_openmm_system(
            _two_chain_system(),
            topology,
            topology,
            tmp_path / "out",
            periodic_chains=True,
            dispersion_correction=False,
            python_executable=sys.executable,
        )


def test_write_openmm_system_rejects_box_narrower_than_twice_cutoff(
    tmp_path: Path,
) -> None:
    narrow = GlycamMappedSystem(
        **{
            **_two_chain_system().__dict__,
            "basis_vectors": np.diag([24.0, 10.4, 30.0]),
        }
    )

    with pytest.raises(ValueError, match="twice the 1 nm nonbonded cutoff"):
        write_openmm_system(
            narrow,
            tmp_path / "system.top",
            tmp_path / "structure.gro",
            tmp_path / "out",
            periodic_chains=True,
            dispersion_correction=False,
            python_executable=sys.executable,
        )


def test_write_openmm_system_reports_builder_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    topology = tmp_path / "system.top"
    topology.touch()

    def run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(
            args=command, returncode=1, stdout="", stderr="No module named 'openmm'"
        )

    monkeypatch.setattr("cellulose_forge.openmm.subprocess.run", run)
    with pytest.raises(RuntimeError, match="No module named 'openmm'"):
        write_openmm_system(
            _two_chain_system(),
            topology,
            topology,
            tmp_path / "out",
            periodic_chains=False,
            dispersion_correction=False,
            python_executable=sys.executable,
        )


def test_validation_platform_avoids_reference_for_tilted_third_vector() -> None:
    from cellulose_forge.openmm import validation_platform_for_box  # noqa: PLC0415

    upright = [[31.2, 0.0, 0.0], [-11.27, 24.39, 0.0], [0.0, 0.0, 21.0]]
    tilted = [[31.2, 0.0, 0.0], [-11.27, 24.39, 0.0], [-11.23, -0.79, 21.03]]
    assert validation_platform_for_box(upright) == "Reference"
    assert validation_platform_for_box(tilted) == "CPU"
