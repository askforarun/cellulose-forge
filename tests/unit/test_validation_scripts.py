"""Validation entry points must remain importable after helper refactors."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
VALIDATION = ROOT / "validation"
OPTIONAL_IMPORTS = frozenset({"openff", "openmm", "parmed", "rdkit"})


@pytest.mark.parametrize("script", sorted(VALIDATION.glob("*.py")), ids=lambda p: p.name)
def test_validation_script_imports(script: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Import a validator without executing its command-line entry point."""
    monkeypatch.syspath_prepend(str(ROOT / "src"))
    monkeypatch.syspath_prepend(str(VALIDATION))
    module_name = f"_cellulose_forge_validation_import_{script.stem}"
    spec = importlib.util.spec_from_file_location(module_name, script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    except ModuleNotFoundError as error:
        missing_root = (error.name or "").split(".", maxsplit=1)[0]
        if missing_root in OPTIONAL_IMPORTS:
            pytest.skip(f"optional validation dependency is not installed: {error.name}")
        raise
    finally:
        sys.modules.pop(module_name, None)
