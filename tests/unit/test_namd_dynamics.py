"""Tests for the real-NAMD production dynamics acceptance analysis."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def _write_fixture(directory: Path, *, final_bond_angstrom: float = 1.1) -> Path:
    directory.mkdir()
    (directory / "system.psf").write_text(
        """PSF EXT

       2 !NATOM
       1 SEG 1 BGLC C1 CT 0.0 12.011 0
       2 SEG 1 BGLC H1 HA 0.0 1.008 0

       1 !NBOND: bonds
       1       2
""",
        encoding="utf-8",
    )
    (directory / "single_point_coordinates.gro").write_text(
        "fixture\n    2\n    1BGLC    C1    1   0.100   0.100   0.100\n"
        "    1BGLC    H1    2   0.210   0.100   0.100\n"
        "   3.00000   3.00000   3.00000\n",
        encoding="utf-8",
    )

    def pdb_atom(serial: int, name: str, x_coordinate: float, element: str) -> str:
        return (
            f"{'ATOM':<6}{serial:5d} {name:^4s} {'BGLC':>4s}{'A'}{1:4d}    "
            f"{x_coordinate:8.3f}{1.0:8.3f}{1.0:8.3f}{1.0:6.2f}{0.0:6.2f}"
            f"      {'SEG':<4s}{element:>2s}"
        )

    (directory / "namd.coor").write_text(
        "REMARK final\n"
        f"{pdb_atom(1, 'C1', 1.0, 'C')}\n"
        f"{pdb_atom(2, 'H1', 1.0 + final_bond_angstrom, 'H')}\n",
        encoding="utf-8",
    )
    (directory / "namd.xsc").write_text(
        "# box\n1500 30 0 0 0 30 0 0 0 30 15 15 15\n", encoding="utf-8"
    )
    energy_rows = "\n".join(
        f"ENERGY: {step} 1 1 1 0 -1 -1 0 0 1 2 300 1" for step in range(1001, 1501)
    )
    (directory / "production_dynamics.log").write_text(
        "ETITLE: TS BOND ANGLE DIHED IMPRP ELECT VDW BOUNDARY MISC KINETIC "
        "TOTAL TEMP POTENTIAL\n"
        f"{energy_rows}\n[Partition 0][Node 0] End of program\n",
        encoding="utf-8",
    )
    return directory


def _script() -> Path:
    return Path(__file__).resolve().parents[2] / "validation" / "namd" / "analyze_dynamics.py"


def test_namd_dynamics_analysis_accepts_completed_intact_run(tmp_path: Path) -> None:
    build = _write_fixture(tmp_path / "build")
    output = tmp_path / "result.json"
    subprocess.run(
        [sys.executable, str(_script()), str(build), "--output", str(output)],
        check=True,
        capture_output=True,
        text=True,
    )
    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["passed"] is True
    assert result["final_step"] == 1500
    assert result["dynamics_energy_rows"] == 500
    assert result["bond_length_range_angstrom"] == [1.1, 1.1]


def test_namd_dynamics_analysis_rejects_broken_bond(tmp_path: Path) -> None:
    build = _write_fixture(tmp_path / "build", final_bond_angstrom=2.0)
    completed = subprocess.run(
        [sys.executable, str(_script()), str(build)],
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode != 0
    result = json.loads(completed.stdout)
    assert result["criteria"]["all_bonds_are_0.85_to_1.75_angstrom"] is False
