from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

from cellulose_forge.energy_drivers.report import EnergyReport


def _load_validator() -> ModuleType:
    path = Path(__file__).resolve().parents[2] / "validation" / "validate_openmm_engine.py"
    spec = importlib.util.spec_from_file_location("validate_openmm_engine", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_energy_comparisons_passed_requires_every_case_to_pass() -> None:
    validator = _load_validator()

    assert validator._energy_comparisons_passed(
        {
            "glycam": {"energy_comparison": {"passed": True}},
            "charmm36": {"energy_comparison": {"passed": True}},
        }
    )
    assert not validator._energy_comparisons_passed(
        {
            "glycam": {"energy_comparison": {"passed": True}},
            "charmm36": {"energy_comparison": {"passed": False}},
        }
    )


def test_dispersion_partition_uses_tight_combined_vdw_gate() -> None:
    validator = _load_validator()
    common = {
        "bond": 0.0,
        "angle": 0.0,
        "dihedral": 0.0,
        "improper": 0.0,
        "electrostatics": 0.0,
        "restraint": 0.0,
        "potential": -8000.0,
    }
    gromacs = EnergyReport(
        engine="gromacs",
        coordinate_sha256="abc",
        terms={**common, "vdw": -6700.0, "dispersion_correction": -1300.0},
    )
    openmm = EnergyReport(
        engine="openmm",
        coordinate_sha256="abc",
        terms={**common, "vdw": -6698.9, "dispersion_correction": -1301.1},
    )

    comparison = validator._compare_openmm_gromacs_reports(gromacs, openmm)

    assert comparison["common_terms"]["dispersion_correction"]["passed"] is True
    assert comparison["combined_vdw_and_dispersion_correction"]["passed"] is True
    assert comparison["passed"] is True
