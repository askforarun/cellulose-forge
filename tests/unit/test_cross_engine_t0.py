from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path


def _energy(engine: str, checksum: str, bond: float = 100.0) -> dict[str, object]:
    return {
        "schema_version": 1,
        "engine": engine,
        "units": "kJ/mol",
        "coordinate_sha256": checksum,
        "terms": {"bond": bond, "potential": 200.0},
        "corrections": [],
    }


def test_cross_engine_t0_runner_compares_subfolders(tmp_path: Path) -> None:
    source = Path("validation/cross_engine_t0/run_t0_validation.py").resolve()
    runner = tmp_path / "run_t0_validation.py"
    runner.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
    shared = tmp_path / "shared"
    shared.mkdir()
    coordinates = shared / "coordinates.g96"
    coordinates.write_text("identical coordinates\n", encoding="utf-8")
    checksum = hashlib.sha256(coordinates.read_bytes()).hexdigest()
    engines: dict[str, object] = {}
    for engine in ("gromacs", "openmm", "lammps", "namd", "amber"):
        directory = tmp_path / engine
        directory.mkdir()
        (directory / "energy.json").write_text(
            json.dumps(_energy(engine, checksum)) + "\n", encoding="utf-8"
        )
        engines[engine] = {"directory": engine, "command": ["true"]}
    manifest = {
        "schema_version": 1,
        "system": "test system",
        "coordinate_file": "shared/coordinates.g96",
        "coordinate_sha256": checksum,
        "reference_engine": "gromacs",
        "absolute_tolerance_kj_mol": 0.05,
        "relative_tolerance": 1.0e-5,
        "required_terms": ["bond", "potential"],
        "engines": engines,
    }
    (tmp_path / "manifest.json").write_text(json.dumps(manifest) + "\n", encoding="utf-8")

    completed = subprocess.run(
        [sys.executable, str(runner), "--no-run"],
        cwd=tmp_path,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    report = json.loads((tmp_path / "t0_energy_comparison.json").read_text(encoding="utf-8"))
    assert report["passed"] is True
    assert set(report["comparisons"]) == {"amber", "lammps", "namd", "openmm"}


def test_cross_engine_t0_runner_rejects_energy_difference(tmp_path: Path) -> None:
    source = Path("validation/cross_engine_t0/run_t0_validation.py").resolve()
    runner = tmp_path / "run_t0_validation.py"
    runner.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
    coordinates = tmp_path / "coordinates.g96"
    coordinates.write_text("identical coordinates\n", encoding="utf-8")
    checksum = hashlib.sha256(coordinates.read_bytes()).hexdigest()
    engines: dict[str, object] = {}
    for engine, bond in (("gromacs", 100.0), ("openmm", 101.0)):
        directory = tmp_path / engine
        directory.mkdir()
        (directory / "energy.json").write_text(
            json.dumps(_energy(engine, checksum, bond)) + "\n", encoding="utf-8"
        )
        engines[engine] = {"directory": engine, "command": ["true"]}
    manifest = {
        "schema_version": 1,
        "system": "test system",
        "coordinate_file": "coordinates.g96",
        "coordinate_sha256": checksum,
        "reference_engine": "gromacs",
        "absolute_tolerance_kj_mol": 0.05,
        "relative_tolerance": 1.0e-5,
        "required_terms": ["bond", "potential"],
        "engines": engines,
    }
    (tmp_path / "manifest.json").write_text(json.dumps(manifest) + "\n", encoding="utf-8")

    completed = subprocess.run(
        [sys.executable, str(runner), "--no-run"],
        cwd=tmp_path,
        check=False,
    )
    assert completed.returncode == 1
    report = json.loads((tmp_path / "t0_energy_comparison.json").read_text(encoding="utf-8"))
    assert report["passed"] is False
    assert report["comparisons"]["openmm"]["common_terms"]["bond"]["passed"] is False
