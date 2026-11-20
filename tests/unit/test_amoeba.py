"""Tests for the experimental externally supplied AMOEBA backend."""

from __future__ import annotations

import hashlib
import json
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
    runner = amoeba.render_openmm_runner()
    assert "expected one AmoebaVdwForce" in runner
    assert "setUseDispersionCorrection(False)" in runner
    assert "getPlatformByName" not in runner


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
    "geometry",
    [
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
        ["--shape", "fibril", "--fibril-model", "18-chain"],
    ],
)
def test_cli_writes_openmm_readable_amoeba_bundle(tmp_path: Path, geometry: list[str]) -> None:
    assert _SOURCE is not None
    output = tmp_path / "amoeba"
    arguments = [
        "build",
        "--allomorph",
        "I-beta",
        *geometry,
        "--vacuum",
        "30",
        "--chain-boundary",
        "finite",
        "--degree-polymerization",
        "6",
        "--padding-longitudinal",
        "20",
        "--force-field",
        "amoeba",
        "--md-engine",
        "openmm",
        "--md-engine",
        "tinker",
        "--amoeba-source-directory",
        _SOURCE,
        "--output",
        str(output),
    ]
    assert main(arguments) == 0
    report = json.loads((output / "build_report.json").read_text(encoding="utf-8"))
    assert report["source_commit"] == amoeba.AMOEBA_SOURCE_COMMIT
    assert report["chain_count"] == 18
    assert report["md_engines"] == ["openmm", "tinker"]
    assert (output / "structure.xyz").read_text(encoding="utf-8").split()[0] == "2322"
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
    assert system.getNumParticles() == 2322
    vdw_forces = [force for force in system.getForces() if isinstance(force, openmm.AmoebaVdwForce)]
    assert len(vdw_forces) == 1
    assert vdw_forces[0].getUseDispersionCorrection() is False

    tinker_only = arguments.copy()
    openmm_index = tinker_only.index("openmm")
    del tinker_only[openmm_index - 1 : openmm_index + 1]
    assert main(tinker_only) == 0
    assert not runner_path.exists()
    rebuilt_report = json.loads((output / "build_report.json").read_text(encoding="utf-8"))
    assert rebuilt_report["outputs"]["openmm_runner"] is None
