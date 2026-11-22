"""Regression checks for committed chemically complete coordinate examples."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[2]


def _pdb_names(path: Path) -> list[str]:
    return [
        line[12:16].strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]


def _pdb_header(path: Path) -> list[str]:
    return [
        line for line in path.read_text(encoding="utf-8").splitlines() if line.startswith("REMARK")
    ]


@pytest.mark.parametrize(
    "relative_path",
    [
        "examples/md_inputs/cellulose_ibeta_110_charmm36_periodic/crystal_structure.pdb",
        "examples/md_inputs/cellulose_ibeta_110_glycam_finite/crystal_structure.pdb",
        "examples/fibrils/ibeta_36chain_dp12.pdb",
    ],
)
def test_shipped_i_beta_diagnostic_headers_are_crystallographic(relative_path: str) -> None:
    header = _pdb_header(_ROOT / relative_path)
    assert header[2] == "REMARK  chain and residue IDs follow the crystallographic model"
    assert header[3] == "REMARK  unit-cell chain A=origin chain B=center for cellulose I-beta"


def _gro_atoms(path: Path) -> list[tuple[int, str, np.ndarray]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    atom_count = int(lines[1])
    return [
        (
            int(line[:5]),
            line[10:15].strip(),
            np.array([float(line[20:28]), float(line[28:36]), float(line[36:44])]),
        )
        for line in lines[2 : atom_count + 2]
    ]


@pytest.mark.parametrize(
    "stem",
    [
        "ibeta_010_slab",
        "ibeta_020_slab",
        "ibeta_100_slab",
        "ibeta_110_slab",
        "ibeta_1m10_slab",
        "ibeta_200_slab",
        "ibeta_m1m10_slab",
    ],
)
def test_finite_face_examples_are_capped_and_longitudinally_separated(stem: str) -> None:
    directory = _ROOT / "examples/faces"
    names = _pdb_names(directory / f"{stem}.pdb")
    assert len(names) == 3096
    assert {name: names.count(name) for name in ("O1", "HO1", "HO4")} == {
        "O1": 24,
        "HO1": 24,
        "HO4": 24,
    }
    assert _pdb_header(directory / f"{stem}.pdb")[2] == (
        "REMARK  chain IDs follow the crystal; residues run reducing to non-reducing end"
    )
    report = json.loads((directory / f"{stem}.build_report.json").read_text(encoding="utf-8"))
    assert report["chemical_atom_count"] == 3096
    assert any("padding_z=12.000000" in item for item in report["geometry_annotations"])


@pytest.mark.parametrize(
    ("folder", "terminal_hydrogen", "expected_oh_nm"),
    [
        ("cellulose_ibeta_110_glycam_finite", "H4O", 0.097),
        ("chitin_alpha_100_charmm36_finite", "HO4", 0.096),
    ],
)
def test_finite_md_example_pdb_matches_simulation_atom_count(
    folder: str, terminal_hydrogen: str, expected_oh_nm: float
) -> None:
    directory = _ROOT / "examples/md_inputs" / folder
    pdb_names = _pdb_names(directory / "structure.pdb")
    gro_lines = (directory / "structure.gro").read_text(encoding="utf-8").splitlines()
    gro_count = int(gro_lines[1])
    assert len(pdb_names) == gro_count
    assert pdb_names == [line[10:15].strip() for line in gro_lines[2 : gro_count + 2]]
    assert {"O1", "HO1", terminal_hydrogen} <= set(pdb_names)
    header = _pdb_header(directory / "structure.pdb")
    assert header[2] == (
        "REMARK  chain IDs follow the crystal; residues run reducing to non-reducing end"
    )
    assert (directory / "crystal_structure.pdb").is_file()
    assert len(_pdb_names(directory / "crystal_structure.pdb")) < gro_count

    atoms = _gro_atoms(directory / "structure.gro")
    coordinates = {(residue, name): coordinate for residue, name, coordinate in atoms}
    terminal_residues = [residue for residue, name, _ in atoms if name == terminal_hydrogen]
    assert terminal_residues
    for residue in terminal_residues:
        distance_nm = np.linalg.norm(
            coordinates[residue, terminal_hydrogen] - coordinates[residue, "O4"]
        )
        assert distance_nm == pytest.approx(expected_oh_nm, abs=0.001)


@pytest.mark.parametrize(
    "folder",
    ["cellulose_ibeta_110_charmm36_periodic", "chitin_beta_100_charmm36_periodic"],
)
def test_periodic_charmm36_example_pdb_matches_simulation_atom_order(folder: str) -> None:
    directory = _ROOT / "examples/md_inputs" / folder
    pdb_names = _pdb_names(directory / "structure.pdb")
    gro_lines = (directory / "structure.gro").read_text(encoding="utf-8").splitlines()
    gro_count = int(gro_lines[1])
    assert len(pdb_names) == gro_count
    assert pdb_names == [line[10:15].strip() for line in gro_lines[2 : gro_count + 2]]
    assert not {"O1", "HO1", "HO4"} & set(pdb_names)
    header = _pdb_header(directory / "structure.pdb")
    assert header[2] == (
        "REMARK  chain IDs follow the crystal; residue numbers follow mapped topology"
    )
    crystal_path = directory / "crystal_structure.pdb"
    assert crystal_path.is_file()
    assert len(_pdb_names(crystal_path)) == gro_count
