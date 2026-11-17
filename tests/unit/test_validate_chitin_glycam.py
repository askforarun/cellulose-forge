from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest


def _load_validator() -> ModuleType:
    path = Path(__file__).resolve().parents[2] / "validation" / "validate_chitin_glycam.py"
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location("validate_chitin_glycam", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_amber_periodic_improper_energy_comes_from_native_topology(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    validator = _load_validator()
    commands: list[list[str]] = []

    def fake_run(command: list[str], directory: Path) -> str:
        commands.append(command)
        assert directory == tmp_path
        return 'CELLULOSE_FORGE_AMBER_IMPROPERS={"count": 48, "energy_kcal_mol": 2.5}\n'

    monkeypatch.setattr(validator, "_run", fake_run)
    result = validator._amber_periodic_improper_energy(
        tmp_path,
        "parmed-python",
        tmp_path / "structure.prmtop",
        tmp_path / "identical.inpcrd",
    )

    assert result == {"count": 48, "energy_kcal_mol": 2.5}
    assert commands[0][0] == "parmed-python"
    assert commands[0][-2:] == [
        str(tmp_path / "structure.prmtop"),
        str(tmp_path / "identical.inpcrd"),
    ]


@pytest.mark.parametrize(
    "payload",
    [
        "",
        'CELLULOSE_FORGE_AMBER_IMPROPERS={"count": -1, "energy_kcal_mol": 2.5}\n',
        'CELLULOSE_FORGE_AMBER_IMPROPERS={"count": 48, "energy_kcal_mol": "2.5"}\n',
    ],
)
def test_amber_periodic_improper_classification_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, payload: str
) -> None:
    validator = _load_validator()
    monkeypatch.setattr(validator, "_run", lambda command, directory: payload)

    with pytest.raises(RuntimeError, match="Amber periodic-improper"):
        validator._amber_periodic_improper_energy(
            tmp_path,
            "parmed-python",
            tmp_path / "structure.prmtop",
            tmp_path / "identical.inpcrd",
        )
