from __future__ import annotations

import importlib.util
import struct
from pathlib import Path
from types import ModuleType


def _load_validator() -> ModuleType:
    path = Path(__file__).resolve().parents[2] / "validation" / "validate_drude_cross_engine.py"
    spec = importlib.util.spec_from_file_location("validate_drude_cross_engine", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_displace_drudes_changes_only_drude_coordinates(tmp_path: Path) -> None:
    validator = _load_validator()
    path = tmp_path / "coordinates.coor"
    original_values = [0.0] * 12
    original = struct.pack("=i", 4) + struct.pack("=12d", *original_values)
    path.write_bytes(original)

    returned = validator._displace_drudes(path, [1, 3])
    changed = struct.unpack("=12d", path.read_bytes()[4:])

    assert returned == original
    assert changed[3] == 0.01
    assert changed[10] == 0.01
    assert sum(value != 0.0 for value in changed) == 2


def test_displaced_state_compares_engine_bookkeeping_equivalents() -> None:
    validator = _load_validator()
    namd = {
        "bond": 120.0,
        "angle": 200.0,
        "dihedral": 300.0,
        "improper": 0.0,
        "cross": 0.0,
        "nonbonded": 80.0,
        "potential": 700.0,
    }
    openmm = {
        "terms_kj_mol": {
            "bond": 100.0,
            "angle": 200.0,
            "dihedral": 300.0,
            "improper": 0.0,
            "cross": 0.0,
            "nonbonded": 70.0,
            "drude": 30.0,
            "potential": 700.0,
        },
        "drude_spring_anisotropy_kj_mol": 20.0,
        "forces": [],
    }

    comparison = validator._compare_displaced_state(namd, openmm)

    assert comparison["comparisons"]["bond_plus_spring_anisotropy"]["passed"] is True
    assert comparison["comparisons"]["nonbonded_plus_thole"]["passed"] is True
    assert comparison["passed"] is True


def test_displacement_difference_gates_the_polarization_response() -> None:
    validator = _load_validator()
    reference_namd = {
        "bond": 100.0,
        "angle": 200.0,
        "dihedral": 300.0,
        "improper": 0.0,
        "cross": 0.0,
        "nonbonded": 100.0,
        "potential": 700.0,
    }
    displaced_namd = {
        "bond": 120.0,
        "angle": 200.0,
        "dihedral": 300.0,
        "improper": 0.0,
        "cross": 0.0,
        "nonbonded": 80.0,
        "potential": 700.0,
    }
    reference_openmm = {
        "terms_kj_mol": {**reference_namd, "drude": 0.0},
        "drude_spring_anisotropy_kj_mol": 0.0,
        "forces": [],
    }
    displaced_openmm = {
        "terms_kj_mol": {
            "bond": 100.0,
            "angle": 200.0,
            "dihedral": 300.0,
            "improper": 0.0,
            "cross": 0.0,
            "nonbonded": 70.0,
            "drude": 30.0,
            "potential": 700.0,
        },
        "drude_spring_anisotropy_kj_mol": 20.0,
        "forces": [],
    }

    comparison = validator._compare_polarization_response(
        reference_namd,
        reference_openmm,
        displaced_namd,
        displaced_openmm,
    )

    assert comparison["comparisons"]["spring_and_anisotropy"]["passed"] is True
    assert comparison["comparisons"]["total_polarization"]["passed"] is True
    assert comparison["passed"] is True


def test_displacement_difference_fails_on_response_scaled_gap() -> None:
    validator = _load_validator()
    reference_namd = {
        "bond": 100.0,
        "angle": 200.0,
        "dihedral": 300.0,
        "improper": 0.0,
        "cross": 0.0,
        "nonbonded": 100.0,
        "potential": 700.0,
    }
    displaced_namd = {**reference_namd, "bond": 120.0, "potential": 720.0}
    reference_openmm = {
        "terms_kj_mol": {**reference_namd, "drude": 0.0},
        "drude_spring_anisotropy_kj_mol": 0.0,
        "forces": [],
    }
    displaced_openmm = {
        "terms_kj_mol": {
            **reference_namd,
            "drude": 30.0,
            "potential": 730.0,
        },
        "drude_spring_anisotropy_kj_mol": 20.0,
        "forces": [],
    }

    comparison = validator._compare_polarization_response(
        reference_namd,
        reference_openmm,
        displaced_namd,
        displaced_openmm,
    )

    assert comparison["comparisons"]["total_polarization"]["passed"] is False
    assert comparison["passed"] is False


def test_displaced_state_requires_nonzero_drude_energy() -> None:
    validator = _load_validator()
    namd = {
        "bond": 100.0,
        "angle": 200.0,
        "dihedral": 300.0,
        "improper": 0.0,
        "cross": 0.0,
        "nonbonded": 80.0,
        "potential": 680.0,
    }
    openmm = {
        "terms_kj_mol": {**namd, "drude": 0.0},
        "drude_spring_anisotropy_kj_mol": 0.0,
        "forces": [],
    }

    comparison = validator._compare_displaced_state(namd, openmm)

    assert comparison["passed"] is False


def test_namd_single_point_disables_dynamics_only_conventions() -> None:
    validator = _load_validator()

    configuration = validator._namd_configuration(
        [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]]
    )

    assert "switching off" in configuration
    assert "LJcorrection off" in configuration
    assert "drudeHardWall off" in configuration
    assert "temperature 0" in configuration
    assert configuration.endswith("run 0\n")
