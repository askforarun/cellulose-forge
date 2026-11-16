from __future__ import annotations

import hashlib
import json
import os
import runpy
import subprocess
import sys
from argparse import Namespace
from pathlib import Path

import pytest

from cellulose_forge.energy_drivers import EnergyReport


def _energy(engine: str, checksum: str, bond: float = 100.0) -> dict[str, object]:
    return {
        "schema_version": 1,
        "engine": engine,
        "units": "kJ/mol",
        "coordinate_sha256": checksum,
        "terms": {"bond": bond, "potential": bond},
        "corrections": [],
    }


def test_cross_engine_single_point_runner_compares_subfolders(tmp_path: Path) -> None:
    source = Path("validation/cross_engine_single_point/run_single_point_validation.py").resolve()
    runner = tmp_path / "run_single_point_validation.py"
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
        env={**os.environ, "PYTHONPATH": str(Path("src").resolve())},
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    report = json.loads(
        (tmp_path / "single_point_energy_comparison.json").read_text(encoding="utf-8")
    )
    assert report["passed"] is True
    assert set(report["comparisons"]) == {"amber", "lammps", "namd", "openmm"}


def test_cross_engine_single_point_runner_rejects_energy_difference(tmp_path: Path) -> None:
    source = Path("validation/cross_engine_single_point/run_single_point_validation.py").resolve()
    runner = tmp_path / "run_single_point_validation.py"
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
        env={**os.environ, "PYTHONPATH": str(Path("src").resolve())},
        check=False,
    )
    assert completed.returncode == 1
    report = json.loads(
        (tmp_path / "single_point_energy_comparison.json").read_text(encoding="utf-8")
    )
    assert report["passed"] is False
    assert report["comparisons"]["openmm"]["common_terms"]["bond"]["passed"] is False


def test_selected_engines_cannot_omit_required_engine() -> None:
    namespace = runpy.run_path(
        "validation/cross_engine_single_point/run_single_point_validation.py",
        run_name="energy_runner_test",
    )
    arguments = Namespace(selected_engines=["gromacs", "openmm"])
    manifest = {
        "reference_engine": "gromacs",
        "engines": {
            "gromacs": {"optional": False},
            "openmm": {"optional": False},
            "lammps": {"optional": False},
            "amber": {"optional": True},
        },
    }

    with pytest.raises(ValueError, match=r"omit required.*lammps"):
        namespace["_selected_engines"](arguments, manifest)


def test_all_builtin_adapters_write_reports(  # noqa: PLR0915
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    namespace = runpy.run_path(
        "validation/cross_engine_single_point/run_single_point_validation.py",
        run_name="energy_runner_adapter_test",
    )
    run_adapter = namespace["_run_adapter"]
    globals_ = run_adapter.__globals__
    calls: dict[str, object] = {}

    def fake_gromacs(directory: Path, **kwargs: object) -> EnergyReport:
        calls["gromacs"] = (directory, kwargs)
        return EnergyReport("gromacs", "sha", {"bond": 1.0, "potential": 1.0})

    def fake_openmm(**kwargs: object) -> EnergyReport:
        calls["openmm"] = kwargs
        return EnergyReport("openmm", "sha", {"bond": 1.0, "potential": 1.0})

    def fake_namd(directory: Path, **kwargs: object) -> EnergyReport:
        calls["namd"] = (directory, kwargs)
        return EnergyReport("namd", "sha", {"bond": 1.0, "potential": 1.0})

    def fake_lammps(directory: Path, **kwargs: object) -> EnergyReport:
        calls["lammps"] = (directory, kwargs)
        return EnergyReport("lammps", "sha", {"bond": 1.0, "potential": 1.0})

    def fake_amber(directory: Path, **kwargs: object) -> EnergyReport:
        calls["amber"] = (directory, kwargs)
        return EnergyReport("amber", "sha", {"bond": 1.0, "potential": 1.0})

    monkeypatch.setitem(globals_, "run_gromacs_single_point", fake_gromacs)
    monkeypatch.setitem(globals_, "run_openmm_single_point", fake_openmm)
    monkeypatch.setitem(globals_, "run_namd_single_point", fake_namd)
    monkeypatch.setitem(globals_, "run_lammps_single_point", fake_lammps)
    monkeypatch.setitem(globals_, "run_glycam_amber_single_point", fake_amber)

    gromacs = tmp_path / "gromacs"
    openmm = tmp_path / "openmm"
    gromacs.mkdir()
    openmm.mkdir()
    namd = tmp_path / "namd"
    namd.mkdir()
    lammps = tmp_path / "lammps"
    lammps.mkdir()
    amber = tmp_path / "amber"
    amber.mkdir()
    run_adapter(
        "gromacs",
        {
            "adapter": "gromacs",
            "inputs": {
                "coordinate_file": "../shared.gro",
                "topology_file": "system.top",
                "mdp_file": "single_point.mdp",
            },
        },
        gromacs,
        "sha",
    )
    run_adapter(
        "openmm",
        {
            "adapter": "openmm",
            "inputs": {
                "system_file": "system.xml",
                "coordinate_file": "../shared.gro",
                "classification": {
                    "urey_bradley_bonds": {"0": [2]},
                    "periodic_impropers": {"2": []},
                    "rb_impropers": {},
                    "force_roles": {"4": "harmonic_improper"},
                },
            },
        },
        openmm,
        "sha",
    )
    run_adapter(
        "namd",
        {
            "adapter": "namd",
            "inputs": {
                "configuration_file": "single_point_namd.conf",
                "coordinate_file": "../shared.gro",
                "binary_coordinate_file": "single_point.coor",
                "threads": 2,
                "classified_terms": {"MISC": "restraint"},
            },
        },
        namd,
        "sha",
    )
    run_adapter(
        "lammps",
        {
            "adapter": "lammps",
            "inputs": {
                "input_file": "in.cellulose",
                "data_file": "system.data",
                "coordinate_file": "../shared.gro",
                "periodic_improper_count": 2,
                "periodic_improper_kj_mol": 1.25,
                "explicit_pair_kj_mol": {"vdw": 2.0, "electrostatics": 3.0},
            },
        },
        lammps,
        "sha",
    )
    run_adapter(
        "amber",
        {
            "adapter": "amber",
            "inputs": {
                "input_file": "single_point.in",
                "topology_file": "structure.prmtop",
                "amber_coordinate_file": "structure.inpcrd",
                "coordinate_file": "../shared.gro",
                "periodic_improper_kcal_mol": 0.0,
            },
        },
        amber,
        "sha",
    )

    assert json.loads((gromacs / "energy.json").read_text())["engine"] == "gromacs"
    assert json.loads((openmm / "energy.json").read_text())["engine"] == "openmm"
    assert json.loads((namd / "energy.json").read_text())["engine"] == "namd"
    assert json.loads((lammps / "energy.json").read_text())["engine"] == "lammps"
    assert json.loads((amber / "energy.json").read_text())["engine"] == "amber"
    openmm_call = calls["openmm"]
    assert isinstance(openmm_call, dict)
    classification = openmm_call["classification"]
    assert classification.urey_bradley_bonds == {0: (2,)}
    assert classification.periodic_impropers == {2: ()}
    assert classification.force_roles == {4: "harmonic_improper"}
    namd_call = calls["namd"]
    assert isinstance(namd_call, tuple)
    assert namd_call[1]["threads"] == 2
    assert namd_call[1]["classified_terms"] == {"MISC": "restraint"}
    lammps_call = calls["lammps"]
    assert isinstance(lammps_call, tuple)
    assert lammps_call[1]["periodic_improper_count"] == 2
    assert lammps_call[1]["explicit_pair_kj_mol"] == {"vdw": 2.0, "electrostatics": 3.0}
    amber_call = calls["amber"]
    assert isinstance(amber_call, tuple)
    assert amber_call[1]["periodic_improper"] == 0.0


def test_builtin_adapter_refuses_noncanonical_force_index(tmp_path: Path) -> None:
    namespace = runpy.run_path(
        "validation/cross_engine_single_point/run_single_point_validation.py",
        run_name="energy_runner_bad_adapter_test",
    )
    directory = tmp_path / "openmm"
    directory.mkdir()
    with pytest.raises(ValueError, match="force indices"):
        namespace["_run_adapter"](
            "openmm",
            {
                "adapter": "openmm",
                "inputs": {
                    "system_file": "system.xml",
                    "coordinate_file": "coordinates.gro",
                    "classification": {
                        "urey_bradley_bonds": {"01": []},
                        "periodic_impropers": {},
                        "rb_impropers": {},
                        "force_roles": {},
                    },
                },
            },
            directory,
            "sha",
        )


@pytest.mark.parametrize(
    "specification",
    [
        {"directory": "gromacs"},
        {
            "directory": "gromacs",
            "adapter": "gromacs",
            "command": ["true"],
            "inputs": {},
        },
    ],
)
def test_engine_execution_requires_exactly_one_mode(specification: dict[str, object]) -> None:
    namespace = runpy.run_path(
        "validation/cross_engine_single_point/run_single_point_validation.py",
        run_name="energy_runner_execution_mode_test",
    )
    with pytest.raises(ValueError, match="exactly one"):
        namespace["_execution_mode"]("gromacs", specification)


def test_engine_command_must_be_nonempty_strings() -> None:
    namespace = runpy.run_path(
        "validation/cross_engine_single_point/run_single_point_validation.py",
        run_name="energy_runner_command_test",
    )
    with pytest.raises(TypeError, match="non-empty list of strings"):
        namespace["_execution_mode"]("gromacs", {"directory": "gromacs", "command": []})
