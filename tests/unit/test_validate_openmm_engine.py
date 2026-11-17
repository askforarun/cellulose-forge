from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType


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
