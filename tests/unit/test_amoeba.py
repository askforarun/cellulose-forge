"""Tests for the experimental externally supplied AMOEBA backend."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path

import pytest

from cellulose_forge import amoeba
from cellulose_forge.cli import main


def _source(tmp_path: Path, degree_polymerization: int = 2) -> tuple[Path, Path]:
    parameter = tmp_path / amoeba.AMOEBA_PARAMETER_FILE
    parameter.write_text("test parameters\n", encoding="utf-8")
    reference = tmp_path / f"Oligomers/Dynamics_1/Vacuum/B_{degree_polymerization // 2}/glc.xyz_2"
    reference.parent.mkdir(parents=True)
    reference.write_text("test reference\n", encoding="utf-8")
    return parameter, reference


def test_source_files_are_pinned_by_hash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    parameter, reference = _source(tmp_path)
    monkeypatch.setattr(
        amoeba,
        "AMOEBA_PARAMETER_SHA256",
        hashlib.sha256(parameter.read_bytes()).hexdigest(),
    )
    monkeypatch.setitem(
        amoeba._REFERENCE_SHA256,
        2,
        hashlib.sha256(reference.read_bytes()).hexdigest(),
    )

    assert amoeba.validate_amoeba_source_directory(tmp_path, 2) == (parameter, reference)
    reference.write_text("changed\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not the audited file"):
        amoeba.validate_amoeba_source_directory(tmp_path, 2)


def test_source_validation_refuses_missing_files_and_unsupported_dp(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="missing"):
        amoeba.validate_amoeba_source_directory(tmp_path, 2)
    with pytest.raises(ValueError, match="finite DP 2, 4, 6, 8, 10"):
        amoeba.validate_amoeba_source_directory(tmp_path, 3)


def test_openmm_runner_disables_vdw_dispersion_correction() -> None:
    runner = amoeba.render_openmm_runner(
        restraint_atom_indices=(0, 2, 5),
        restraint_force=(1000.0, 500.0, 0.0),
    )
    assert "expected one AmoebaVdwForce" in runner
    assert "setUseDispersionCorrection(False)" in runner
    assert "getPlatformByName" not in runner
    assert "RESTRAINT_ATOM_INDICES = [0, 2, 5]" in runner
    assert "RESTRAINT_FORCE_KJ_MOL_NM2 = (1000.0, 500.0, 0.0)" in runner
    assert "CellulosePositionRestraints" in runner
    compile(runner, "run_amoeba_openmm.py", "exec")


def test_cli_refuses_amoeba_without_external_source(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit):
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
                "3",
                "--chains-per-layer",
                "6",
                "--degree-polymerization",
                "6",
                "--chain-boundary",
                "finite",
                "--force-field",
                "amoeba",
                "--output",
                str(tmp_path / "out"),
            ]
        )
    assert "requires --amoeba-source-directory" in capsys.readouterr().err


_SOURCE = os.environ.get("CELLULOSE_FORGE_AMOEBA_SOURCE")


@pytest.mark.skipif(
    _SOURCE is None,
    reason="set CELLULOSE_FORGE_AMOEBA_SOURCE to the audited external repository",
)
@pytest.mark.parametrize(
    ("allomorph", "geometry", "chain_count"),
    [
        (
            "I-beta",
            [
                "--shape",
                "slab",
                "--face",
                "1",
                "1",
                "0",
                "--layers",
                "3",
                "--chains-per-layer",
                "6",
            ],
            18,
        ),
        (
            "II",
            [
                "--shape",
                "slab",
                "--face",
                "1",
                "1",
                "0",
                "--layers",
                "3",
                "--chains-per-layer",
                "6",
            ],
            18,
        ),
        (
            "III-I",
            [
                "--shape",
                "slab",
                "--face",
                "1",
                "1",
                "0",
                "--layers",
                "3",
                "--chains-per-layer",
                "6",
            ],
            18,
        ),
        ("I-alpha", ["--shape", "fibril", "--fibril-model", "36-chain"], 36),
        ("I-beta", ["--shape", "fibril", "--fibril-model", "18-chain"], 18),
        ("II", ["--shape", "fibril", "--fibril-model", "36-chain"], 36),
    ],
)
def test_cli_writes_openmm_readable_amoeba_bundle(
    tmp_path: Path, allomorph: str, geometry: list[str], chain_count: int
) -> None:
    assert _SOURCE is not None
    output = tmp_path / "amoeba"
    arguments = [
        "build",
        "--allomorph",
        allomorph,
        *geometry,
        "--vacuum",
        "30",
        "--chain-boundary",
        "finite",
        "--degree-polymerization",
        "2",
        "--padding-longitudinal",
        "20",
        "--force-field",
        "amoeba",
        "--md-engine",
        "openmm",
        "--amoeba-source-directory",
        _SOURCE,
        "--restrain-chains",
        "1",
        "--restraint-force",
        "1000",
        "500",
        "0",
        "--output",
        str(output),
    ]
    assert main(arguments) == 0
    report = json.loads((output / "build_report.json").read_text(encoding="utf-8"))
    assert report["source_commit"] == amoeba.AMOEBA_SOURCE_COMMIT
    assert report["allomorph"] == allomorph
    assert report["chain_count"] == chain_count
    assert report["md_engines"] == ["openmm"]
    assert report["position_restraints"]["selected_chain_ids"] == [1]
    assert report["position_restraints"]["force_constants_kj_mol_nm2"] == [
        1000.0,
        500.0,
        0.0,
    ]
    expected_particles = chain_count * 45
    assert (
        int((output / "structure.xyz").read_text(encoding="utf-8").split()[0]) == expected_particles
    )
    runner_path = output / "run_amoeba_openmm.py"
    runner = runner_path.read_text(encoding="utf-8")
    assert "MonteCarloBarostat" not in runner
    assert "LangevinMiddleIntegrator" in runner
    openmm = pytest.importorskip("openmm")
    namespace: dict[str, object] = {"__file__": str(runner_path), "__name__": "amoeba_test"}
    previous = Path.cwd()
    try:
        os.chdir(output)
        exec(runner.split("integrator =", maxsplit=1)[0], namespace)
    finally:
        os.chdir(previous)
    system = namespace["system"]
    assert system.getNumParticles() == expected_particles
    vdw_forces = [force for force in system.getForces() if isinstance(force, openmm.AmoebaVdwForce)]
    assert len(vdw_forces) == 1
    assert vdw_forces[0].getUseDispersionCorrection() is False
    restraint_forces = [
        force
        for force in system.getForces()
        if isinstance(force, openmm.CustomExternalForce)
        and force.getName() == "CellulosePositionRestraints"
    ]
    assert len(restraint_forces) == 1
    heavy_per_chain = sum(
        not name.upper().startswith("H") for name, _atom_type, _residue in amoeba._atom_layout(2)
    )
    assert restraint_forces[0].getNumParticles() == heavy_per_chain
    restrained_particles = [
        restraint_forces[0].getParticleParameters(index)[0] for index in range(heavy_per_chain)
    ]
    particles_per_chain = len(amoeba._atom_layout(2))
    assert len({index // particles_per_chain for index in restrained_particles}) == 1
    expected_local = [
        index
        for index, (name, _atom_type, _residue) in enumerate(amoeba._atom_layout(2))
        if not name.upper().startswith("H")
    ]
    assert [index % particles_per_chain for index in restrained_particles] == expected_local
    integrator = openmm.VerletIntegrator(0.001 * openmm.unit.picoseconds)
    context = openmm.Context(
        system,
        integrator,
        openmm.Platform.getPlatformByName("Reference"),
    )
    context.setPositions(namespace["tinker"].positions)
    energy = (
        context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(openmm.unit.kilojoule_per_mole)
    )
    assert math.isfinite(energy)
    del context, integrator


def test_cli_refuses_i_alpha_slab_before_writing(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    output = tmp_path / "amoeba"
    with pytest.raises(SystemExit):
        main(
            [
                "build",
                "--allomorph",
                "I-alpha",
                "--shape",
                "slab",
                "--face",
                "0",
                "1",
                "0",
                "--layers",
                "3",
                "--chains-per-layer",
                "6",
                "--degree-polymerization",
                "6",
                "--chain-boundary",
                "finite",
                "--vacuum",
                "30",
                "--padding-longitudinal",
                "20",
                "--force-field",
                "amoeba",
                "--amoeba-source-directory",
                str(tmp_path),
                "--output",
                str(output),
            ]
        )
    assert "temporarily unsupported" in capsys.readouterr().err
    assert not output.exists()
