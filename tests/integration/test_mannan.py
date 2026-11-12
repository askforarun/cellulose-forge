"""Coordinate-only Mannan I crystals."""

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
    validate_mannose_chirality,
    validate_pyranose_chirality,
)

RESIDUE_ATOMS = 21
RESIDUE_BONDS = (
    ("C1", "C2"), ("C2", "C3"), ("C3", "C4"), ("C4", "C5"),
    ("C5", "O5"), ("O5", "C1"), ("C5", "C6"), ("C6", "O6"),
    ("C2", "O2"), ("C3", "O3"), ("C4", "O4"),
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


def test_mannan_i_unit_cell_matches_primary_model() -> None:
    definition = load_allomorph("i", polymer="mannan")
    crystal = build_unit_cell(definition)

    assert definition.polymer == "mannan"
    assert definition.space_group == "P 21 21 21"
    assert definition.doi == "10.1021/ma00176a014"
    assert definition.residue_name == "MAN"
    assert crystal.atom_count == 84
    assert len(set(crystal.chain_numbers)) == 2
    assert crystal.atom_count % RESIDUE_ATOMS == 0
    assert np.allclose(
        [definition.unit_cell.a, definition.unit_cell.b, definition.unit_cell.c],
        [8.92, 7.21, 10.27],
    )


def test_mannan_i_bulk_geometry_and_chirality() -> None:
    crystal = build_bulk_crystal(load_allomorph("i", polymer="mannan"), (2, 2, 3))
    residues = _residues(crystal)

    lengths = [
        np.linalg.norm(atoms[first] - atoms[second])
        for atoms in residues.values()
        for first, second in RESIDUE_BONDS
    ]
    assert 1.35 < min(lengths)
    assert max(lengths) < 1.6

    for (chain, residue), atoms in residues.items():
        partner = min(
            np.linalg.norm(_minimum_image(other["C1"] - atoms["O4"], crystal.basis_vectors))
            for (other_chain, other_residue), other in residues.items()
            if other_chain == chain and other_residue != residue
        )
        assert 1.35 < partner < 1.5

    report = validate_mannose_chirality(crystal)
    assert report["status"] == "passed"
    assert report["configuration"] == "beta-D-mannose"
    assert report["residue_count"] == len(residues)


def _angle(first: np.ndarray, centre: np.ndarray, third: np.ndarray) -> float:
    left = first - centre
    right = third - centre
    cosine = np.dot(left, right) / (np.linalg.norm(left) * np.linalg.norm(right))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def test_modelled_hydroxyl_hydrogens_have_tetrahedral_geometry() -> None:
    crystal = build_unit_cell(load_allomorph("i", polymer="mannan"))
    residues = _residues(crystal)
    atoms = residues[(1, 1)]

    for hydrogen, oxygen, carbon in (
        ("HO2", "O2", "C2"),
        ("HO3", "O3", "C3"),
        ("HO6", "O6", "C6"),
    ):
        assert np.linalg.norm(atoms[hydrogen] - atoms[oxygen]) == pytest.approx(0.96, abs=1e-6)
        assert _angle(atoms[carbon], atoms[oxygen], atoms[hydrogen]) == pytest.approx(
            109.5, abs=1e-5
        )


def test_modelled_hydrogens_point_to_reported_acceptors() -> None:
    crystal = build_unit_cell(load_allomorph("i", polymer="mannan"))
    residues = _residues(crystal)
    atoms = residues[(1, 1)]

    acceptors = {
        "HO2": min(
            (other["O5"] for (chain, _), other in residues.items() if chain != 1),
            key=lambda xyz: np.linalg.norm(
                _minimum_image(xyz - atoms["O2"], crystal.basis_vectors)
            ),
        ),
        "HO3": min(
            (
                other["O5"]
                for (chain, residue), other in residues.items()
                if chain == 1 and residue != 1
            ),
            key=lambda xyz: np.linalg.norm(
                _minimum_image(xyz - atoms["O3"], crystal.basis_vectors)
            ),
        ),
        "HO6": min(
            (
                other["O3"]
                for (chain, residue), other in residues.items()
                if chain == 1 and residue != 1
            ),
            key=lambda xyz: np.linalg.norm(
                _minimum_image(xyz - atoms["O6"], crystal.basis_vectors)
            ),
        ),
    }
    expected_oxygen_distances = {"HO2": 3.110, "HO3": 2.385, "HO6": 2.711}
    donors = {"HO2": "O2", "HO3": "O3", "HO6": "O6"}
    for hydrogen, acceptor in acceptors.items():
        oxygen = atoms[donors[hydrogen]]
        donor_acceptor = np.linalg.norm(
            _minimum_image(acceptor - oxygen, crystal.basis_vectors)
        )
        assert donor_acceptor == pytest.approx(expected_oxygen_distances[hydrogen], abs=0.002)
        hydrogen_acceptor = np.linalg.norm(
            _minimum_image(acceptor - atoms[hydrogen], crystal.basis_vectors)
        )
        assert hydrogen_acceptor > 1.3
        assert hydrogen_acceptor < donor_acceptor
        hydrogen_to_donor = _minimum_image(
            oxygen - atoms[hydrogen], crystal.basis_vectors
        )
        hydrogen_to_acceptor = _minimum_image(
            acceptor - atoms[hydrogen], crystal.basis_vectors
        )
        cosine = np.dot(hydrogen_to_donor, hydrogen_to_acceptor) / (
            np.linalg.norm(hydrogen_to_donor) * np.linalg.norm(hydrogen_to_acceptor)
        )
        hydrogen_bond_angle = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
        assert hydrogen_bond_angle >= 138.0


def test_modelled_hydroxyl_hydrogens_have_no_periodic_clashes() -> None:
    crystal = build_unit_cell(load_allomorph("i", polymer="mannan"))
    labels = np.asarray(crystal.atom_labels)
    coordinates = crystal.cartesian_coordinates
    hydroxyl = np.flatnonzero(np.isin(labels, ["HO2", "HO3", "HO6"]))
    hydrogens = np.flatnonzero(np.char.startswith(labels.astype(str), "H"))
    heavy = np.flatnonzero(~np.char.startswith(labels.astype(str), "H"))
    bonded_oxygen = {"HO2": "O2", "HO3": "O3", "HO6": "O6"}

    closest_hydrogen = np.inf
    closest_heavy = np.inf
    for atom in hydroxyl:
        hydrogen_distances = np.asarray(
            [
                np.linalg.norm(
                    _minimum_image(coordinates[other] - coordinates[atom], crystal.basis_vectors)
                )
                for other in hydrogens
                if other != atom
            ]
        )
        closest_hydrogen = min(closest_hydrogen, float(hydrogen_distances.min()))

        for other in heavy:
            if labels[other] == bonded_oxygen[str(labels[atom])]:
                distance_to_hydrogen = np.linalg.norm(
                    _minimum_image(
                        coordinates[other] - coordinates[atom], crystal.basis_vectors
                    )
                )
                if distance_to_hydrogen == pytest.approx(0.96, abs=1e-6):
                    continue
            closest_heavy = min(
                closest_heavy,
                float(
                    np.linalg.norm(
                        _minimum_image(
                            coordinates[other] - coordinates[atom], crystal.basis_vectors
                        )
                    )
                ),
            )

    assert closest_hydrogen >= 1.8
    # The primary Table I O3...O5 distance is only 2.385 A.  Its best
    # tetrahedral HO3 placement gives 1.58 A to O5, hence the explicit 1.55 A
    # floor rather than silently moving an experimental heavy atom.
    assert closest_heavy >= 1.55


def test_glucose_stereochemistry_rejects_mannan_c2_epimer() -> None:
    crystal = build_unit_cell(load_allomorph("i", polymer="mannan"))
    with pytest.raises(ValueError, match="inverted C2"):
        validate_pyranose_chirality(crystal, "beta-D-glucose")


def test_mannan_i_has_no_interchain_heavy_atom_clash() -> None:
    crystal = build_bulk_crystal(load_allomorph("i", polymer="mannan"), (2, 2, 2))
    labels = np.asarray(crystal.atom_labels)
    chains = np.asarray(crystal.chain_numbers)
    heavy = np.flatnonzero(np.asarray([not label.startswith("H") for label in labels]))
    inverse = np.linalg.inv(crystal.basis_vectors)
    closest = np.inf
    for atom in heavy:
        fractional = (
            crystal.cartesian_coordinates[heavy] - crystal.cartesian_coordinates[atom]
        ) @ inverse
        fractional -= np.round(fractional)
        distances = np.linalg.norm(fractional @ crystal.basis_vectors, axis=1)
        distances = distances[chains[heavy] != chains[atom]]
        closest = min(closest, float(distances.min()))
    assert closest > 3.0


def test_mirrored_mannan_is_rejected() -> None:
    crystal = build_unit_cell(load_allomorph("i", polymer="mannan"))
    mirrored = replace(
        crystal,
        cartesian_coordinates=crystal.cartesian_coordinates * np.asarray([-1.0, 1.0, 1.0]),
    )
    with pytest.raises(ValueError, match="inverted C1"):
        validate_mannose_chirality(mirrored)


def test_cli_writes_mannan_i_bulk_and_provenance(tmp_path: Path) -> None:
    output = tmp_path / "mannan_i.pdb"
    assert main(
        [
            "build", "--polymer", "mannan", "--allomorph", "I",
            "--shape", "bulk", "--cells", "2", "2", "2", "--output", str(output),
        ]
    ) == 0

    report = json.loads(output.with_suffix(".build_report.json").read_text(encoding="utf-8"))
    assert " MAN " in output.read_text(encoding="utf-8")
    assert report["polymer"] == "mannan"
    assert report["residue_name"] == "MAN"
    assert report["chain_count"] == 8
    assert report["chirality_validation"]["configuration"] == "beta-D-mannose"
    assert report["crystallographic_source"]["definition"] == "mannan/I"
    assert report["crystallographic_source"]["doi"] == "10.1021/ma00176a014"
    annotation = report["modelled_hydrogen_coordinates"]
    assert annotation["atoms"] == ["HO2", "HO3", "HO6"]
    assert annotation["requires_restrained_relaxation_before_force_field_use"] is True


@pytest.mark.parametrize(
    "arguments",
    [
        ["--polymer", "mannan", "--allomorph", "I", "--shape", "slab", "--face", "1", "0", "0"],
        [
            "--polymer", "mannan", "--allomorph", "I", "--shape", "fibril",
            "--fibril-model", "36-chain",
        ],
        [
            "--polymer", "mannan", "--allomorph", "I", "--shape", "bulk",
            "--cells", "1", "1", "1", "--force-field", "charmm36",
        ],
        ["--polymer", "mannan", "--allomorph", "II", "--shape", "unit-cell"],
    ],
)
def test_cli_refuses_unvalidated_mannan_requests(
    tmp_path: Path, arguments: list[str]
) -> None:
    with pytest.raises(SystemExit, match="2"):
        main(["build", *arguments, "--output", str(tmp_path / "invalid")])
