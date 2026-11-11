"""Coordinate-only alpha- and beta-chitin crystals."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from cellulose_forge.cli import main
from cellulose_forge.crystallography import load_allomorph
from cellulose_forge.geometry import (
    BulkCrystal,
    build_bulk_crystal,
    build_unit_cell,
    validate_pyranose_chirality,
)
from cellulose_forge.mmcif import write_mmcif
from cellulose_forge.pdb import write_pdb

DEFINITIONS = ("alpha", "alpha_o6_minor", "beta")
RESIDUE_ATOMS = 27
RESIDUE_BONDS = (
    ("C1", "C2"), ("C2", "C3"), ("C3", "C4"), ("C4", "C5"), ("C5", "O5"), ("O5", "C1"),
    ("C5", "C6"), ("C6", "O6"), ("C3", "O3"), ("C4", "O4"), ("C2", "N2"), ("N2", "C7"),
    ("C7", "O7"), ("C7", "C8"),
)


def _residues(crystal: BulkCrystal) -> dict[tuple[int, int], dict[str, np.ndarray]]:
    residues: dict[tuple[int, int], dict[str, np.ndarray]] = {}
    for label, chain, residue, xyz in zip(
        crystal.atom_labels,
        crystal.chain_numbers,
        crystal.residue_numbers,
        crystal.cartesian_coordinates,
        strict=True,
    ):
        residues.setdefault((chain, residue), {})[label] = xyz
    return residues


def _minimum_image(delta: np.ndarray, basis: np.ndarray) -> np.ndarray:
    fractional = np.linalg.solve(basis.T, delta)
    return basis.T @ (fractional - np.round(fractional))


@pytest.mark.parametrize(
    ("name", "atoms", "chains"),
    [("alpha", 108, 2), ("alpha_o6_minor", 108, 2), ("beta", 54, 1)],
)
def test_unit_cell_counts_follow_space_group(name: str, atoms: int, chains: int) -> None:
    definition = load_allomorph(name, polymer="chitin")
    cell = build_unit_cell(definition)

    assert definition.polymer == "chitin"
    assert cell.atom_count == atoms
    assert len(set(cell.chain_numbers)) == chains
    assert cell.atom_count % RESIDUE_ATOMS == 0


@pytest.mark.parametrize("name", DEFINITIONS)
def test_bulk_geometry_bonds_and_chirality(name: str) -> None:
    crystal = build_bulk_crystal(load_allomorph(name, polymer="chitin"), (2, 2, 3))
    residues = _residues(crystal)

    lengths = [
        np.linalg.norm(atoms[first] - atoms[second])
        for atoms in residues.values()
        for first, second in RESIDUE_BONDS
    ]
    assert 1.2 < min(lengths)
    assert max(lengths) < 1.6

    for (chain, residue), atoms in residues.items():
        partner = min(
            np.linalg.norm(_minimum_image(other["C1"] - atoms["O4"], crystal.basis_vectors))
            for (other_chain, other_residue), other in residues.items()
            if other_chain == chain and other_residue != residue
        )
        assert 1.35 < partner < 1.5

    report = validate_pyranose_chirality(crystal, "N-acetyl-beta-D-glucosamine")
    assert report["status"] == "passed"
    assert report["residue_count"] == len(residues)


@pytest.mark.parametrize("name", DEFINITIONS)
def test_no_interchain_heavy_atom_clashes(name: str) -> None:
    crystal = build_bulk_crystal(load_allomorph(name, polymer="chitin"), (2, 2, 2))
    heavy = [i for i, label in enumerate(crystal.atom_labels) if not label.startswith("H")]
    coordinates = crystal.cartesian_coordinates[heavy]
    chains = np.asarray(crystal.chain_numbers)[heavy]
    basis = crystal.basis_vectors
    closest = min(
        np.linalg.norm(_minimum_image(coordinates[j] - coordinates[i], basis))
        for i in range(len(heavy))
        for j in range(i + 1, len(heavy))
        if chains[i] != chains[j]
    )
    # The closest interchain contact is the N-H...O=C amide hydrogen bond.
    assert closest > 2.5


def test_mirrored_chitin_is_rejected() -> None:
    crystal = build_unit_cell(load_allomorph("beta", polymer="chitin"))
    mirrored = replace(
        crystal,
        cartesian_coordinates=crystal.cartesian_coordinates * np.asarray([-1.0, 1.0, 1.0]),
    )

    with pytest.raises(ValueError, match="inverted C1"):
        validate_pyranose_chirality(mirrored, "N-acetyl-beta-D-glucosamine")


def test_definitions_record_primary_sources() -> None:
    assert load_allomorph("alpha", polymer="chitin").doi == "10.1021/bm801251e"
    assert load_allomorph("beta", polymer="chitin").doi == "10.1021/ma102240r"


def test_cli_writes_chitin_bulk_with_nag_residues(tmp_path: Path) -> None:
    output = tmp_path / "beta.pdb"

    assert main(
        [
            "build", "--polymer", "chitin", "--allomorph", "beta",
            "--shape", "bulk", "--cells", "2", "2", "2", "--output", str(output),
        ]
    ) == 0

    pdb = output.read_text(encoding="utf-8")
    report = json.loads(output.with_suffix(".build_report.json").read_text(encoding="utf-8"))
    assert " NAG " in pdb
    assert " BGC " not in pdb
    assert report["polymer"] == "chitin"
    assert report["chirality_validation"]["configuration"] == "N-acetyl-beta-D-glucosamine"
    assert report["chain_count"] == 4
    assert "P 1 1 21" in pdb
    assert report["residue_name"] == "NAG"
    assert report["crystallographic_source"]["doi"] == "10.1021/ma102240r"
    assert report["crystallographic_source"]["temperature_k"] == 293.0
    assert "o6_conformation" not in report


def test_default_alpha_report_records_major_o6(tmp_path: Path) -> None:
    output = tmp_path / "alpha.pdb"

    assert main(
        [
            "build", "--polymer", "chitin", "--allomorph", "alpha",
            "--shape", "unit-cell", "--output", str(output),
        ]
    ) == 0

    report = json.loads(output.with_suffix(".build_report.json").read_text(encoding="utf-8"))
    assert report["o6_conformation"] == "major"
    assert report["crystallographic_source"]["definition"] == "chitin/alpha"


def test_python_api_writes_nag_residues(tmp_path: Path) -> None:
    definition = load_allomorph("beta", polymer="chitin")
    crystal = build_bulk_crystal(definition, (1, 1, 2))
    write_pdb(crystal, tmp_path / "beta.pdb")
    write_mmcif(crystal, tmp_path / "beta.cif")

    assert build_unit_cell(definition).residue_name == "NAG"
    assert crystal.residue_name == "NAG"
    assert " NAG " in (tmp_path / "beta.pdb").read_text(encoding="utf-8")
    assert " NAG " in (tmp_path / "beta.cif").read_text(encoding="utf-8")
    assert build_unit_cell(load_allomorph("i_beta")).residue_name == "BGC"


def test_cli_selects_minor_alpha_o6_conformation(tmp_path: Path) -> None:
    output = tmp_path / "alpha_minor.pdb"

    assert main(
        [
            "build", "--polymer", "chitin", "--allomorph", "alpha",
            "--o6-conformation", "minor", "--shape", "unit-cell", "--output", str(output),
        ]
    ) == 0

    report = json.loads(output.with_suffix(".build_report.json").read_text(encoding="utf-8"))
    assert report["o6_conformation"] == "minor"
    assert "P 21 21 21" in output.read_text(encoding="utf-8")


@pytest.mark.parametrize(
    "arguments",
    [
        ["--polymer", "chitin", "--allomorph", "beta", "--shape", "slab", "--face", "1",
         "0", "0", "--layers", "1", "--chains-per-layer", "2", "--degree-polymerization", "4"],
        ["--polymer", "chitin", "--allomorph", "beta", "--shape", "bulk", "--cells", "1", "1",
         "1", "--force-field", "glycam", "--chain-boundary", "periodic"],
        ["--polymer", "chitin", "--allomorph", "I-beta", "--shape", "unit-cell"],
        ["--allomorph", "alpha", "--shape", "unit-cell"],
        ["--polymer", "chitin", "--allomorph", "beta", "--o6-conformation", "minor",
         "--shape", "unit-cell"],
    ],
)
def test_cli_rejects_unsupported_chitin_requests(tmp_path: Path, arguments: list[str]) -> None:
    with pytest.raises(SystemExit, match="2"):
        main(["build", *arguments, "--output", str(tmp_path / "invalid")])
