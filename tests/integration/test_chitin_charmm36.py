"""CHARMM36 chitin builds through the CLI (GROMACS only)."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

from cellulose_forge.cli import main
from cellulose_forge.forcefields import discover_charmm36_directory

_ROOT = Path(__file__).parents[2]


def _force_field() -> Path:
    try:
        return discover_charmm36_directory()
    except FileNotFoundError:
        pytest.skip("the pinned CHARMM36 release is not installed")


def _slab_arguments(allomorph: str, boundary: str, output: Path, *extra: str) -> list[str]:
    padding = ["--padding-longitudinal", "12"] if boundary == "finite" else []
    return [
        "build", "--polymer", "chitin", "--allomorph", allomorph, "--shape", "slab",
        "--face", "1", "0", "0", "--layers", "2", "--chains-per-layer", "4",
        "--degree-polymerization", "6", "--chain-boundary", boundary, "--vacuum", "30",
        *padding, "--force-field", "charmm36", *extra, "--output", str(output),
    ]  # fmt: skip


@pytest.mark.parametrize(
    ("engine", "message"),
    [("openmm", "gromacs only"), ("lammps", "gromacs only")],
)
def test_cli_rejects_unvalidated_chitin_charmm36_engines(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], engine: str, message: str
) -> None:
    arguments = _slab_arguments("beta", "periodic", tmp_path / "out", "--md-engine", engine)
    with pytest.raises(SystemExit, match="2"):
        main(arguments)
    assert message in capsys.readouterr().err


def test_cli_rejects_chitin_charmm_template(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    template = tmp_path / "template.itp"
    template.write_text("", encoding="utf-8")
    arguments = _slab_arguments(
        "beta", "periodic", tmp_path / "out", "--charmm-template", str(template)
    )
    with pytest.raises(SystemExit, match="2"):
        main(arguments)
    assert "--charmm-template is not supported for chitin" in capsys.readouterr().err


@pytest.mark.parametrize("allomorph", ["alpha", "beta"])
@pytest.mark.parametrize("boundary", ["periodic", "finite"])
def test_cli_builds_chitin_charmm36_gromacs_slab(
    tmp_path: Path, allomorph: str, boundary: str
) -> None:
    force_field = _force_field()
    output = tmp_path / f"{allomorph}_{boundary}"
    arguments = _slab_arguments(
        allomorph, boundary, output, "--charmm-force-field-directory", str(force_field)
    )
    assert main(arguments) == 0

    assert (output / "chitin_charmm36.itp").is_file()
    assert not (output / "cellulose_charmm36.itp").exists()
    top = (output / "system.top").read_text(encoding="utf-8")
    assert '#include "chitin_charmm36.itp"' in top
    assert "cellulose" not in top

    report = json.loads((output / "build_report.json").read_text(encoding="utf-8"))
    assert report["polymer"] == "chitin"
    assert report["residue_name"] == "NAG"
    assert report["charmm36_residue"] == "BGLCNA"
    assert report["topology_source"] == f"native BGLCNA residue with {boundary} 14bb linkage"
    assert report["validation_scope"]["engines"] == ["gromacs"]
    assert report["md_engines"] == ["gromacs"]
    assert abs(report["generated_chain_charge"]) < 1e-6

    itp = (output / "chitin_charmm36.itp").read_text(encoding="utf-8")
    impropers = [
        row for row in itp.splitlines() if len(row.split()) == 5 and row.split()[-1] == "2"
    ]
    assert len(impropers) == 2 * 6

    gmx = shutil.which("gmx")
    if gmx is None:
        pytest.skip("GROMACS executable 'gmx' is not installed")
    prefix = "periodic" if boundary == "periodic" else "finite"
    completed = subprocess.run(
        [
            gmx, "grompp", "-f", str(_ROOT / f"validation/{prefix}_pme_minim.mdp"),
            "-c", "structure.gro", "-p", "system.top", "-o", "smoke.tpr",
            "-po", "processed.mdp", "-maxwarn", "0",
        ],
        cwd=output, capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "WARNING" not in completed.stdout + completed.stderr
