"""Tests for crystallographic-to-GLYCAM atom mapping."""

from pathlib import Path

import numpy as np
import pytest

from cellulose_forge.crystallography import MillerIndex, load_allomorph
from cellulose_forge.geometry import build_bulk_crystal, build_chain_parallel_slab
from cellulose_forge.glycam_mapping import (
    map_crystal_to_finite_glycam,
    map_crystal_to_periodic_glycam,
)

REFERENCE_PDB = Path("examples/glycam/dp12_reference/cellulose_chain.pdb")
PERIODIC_REFERENCE_PDB = Path(
    "examples/glycam/dp12_periodic_reference/cellulose_periodic_chain.pdb"
)


def _dp12_slab(chains: int = 2, *, periodic: bool = False):
    return build_chain_parallel_slab(
        load_allomorph("i_beta"),
        MillerIndex(1, 1, 0),
        layers=1,
        chains_per_layer=chains,
        degree_polymerization=12,
        vacuum=30.0,
        padding_z=0.0 if periodic else 12.0,
        chain_boundary="periodic" if periodic else "finite",
    )


def _angle(left: np.ndarray, centre: np.ndarray, right: np.ndarray) -> float:
    first = left - centre
    second = right - centre
    cosine = np.dot(first, second) / (np.linalg.norm(first) * np.linalg.norm(second))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def test_dp12_mapping_matches_validated_glycam_atom_and_residue_order() -> None:
    mapped = map_crystal_to_finite_glycam(_dp12_slab(), REFERENCE_PDB)

    assert mapped.atom_count == 2 * 255
    assert mapped.chain_count == 2
    assert mapped.atom_names[:5] == ("HO1", "O1", "C1", "H1", "O5")
    assert mapped.residue_names[:2] == ("ROH", "ROH")
    assert set(mapped.residue_numbers) == set(range(1, 14))
    assert mapped.residue_names[-22:] == ("0GB",) * 22
    assert mapped.atom_names[-22:] == (
        "C1",
        "H1",
        "O5",
        "C5",
        "H5",
        "C6",
        "H62",
        "H61",
        "O6",
        "H6O",
        "C4",
        "H4",
        "O4",
        "H4O",
        "C3",
        "H3",
        "O3",
        "H3O",
        "C2",
        "H2",
        "O2",
        "H2O",
    )


def test_mapping_reverses_crystal_residues_into_glycam_linkage_order() -> None:
    crystal = _dp12_slab()
    mapped = map_crystal_to_finite_glycam(crystal, REFERENCE_PDB)

    crystal_labels = np.asarray(crystal.atom_labels)
    crystal_residues = np.asarray(crystal.residue_numbers)
    expected_reducing_c1 = crystal.cartesian_coordinates[
        (crystal_residues == 12) & (crystal_labels == "C1")
    ][0]
    expected_nonreducing_c1 = crystal.cartesian_coordinates[
        (crystal_residues == 1) & (crystal_labels == "C1")
    ][0]
    mapped_names = np.asarray(mapped.atom_names)
    mapped_residues = np.asarray(mapped.residue_numbers)

    assert np.allclose(
        mapped.cartesian_coordinates[(mapped_residues == 2) & (mapped_names == "C1")][0],
        expected_reducing_c1,
    )
    assert np.allclose(
        mapped.cartesian_coordinates[(mapped_residues == 13) & (mapped_names == "C1")][0],
        expected_nonreducing_c1,
    )


def test_constructed_terminal_bonds_have_chemical_lengths() -> None:
    mapped = map_crystal_to_finite_glycam(_dp12_slab(), REFERENCE_PDB)
    names = np.asarray(mapped.atom_names)
    residues = np.asarray(mapped.residue_numbers)

    def coordinate(residue: int, name: str) -> np.ndarray:
        return mapped.cartesian_coordinates[(residues == residue) & (names == name)][0]

    assert np.linalg.norm(coordinate(1, "O1") - coordinate(2, "C1")) == pytest.approx(
        1.43, abs=0.08
    )
    assert np.linalg.norm(coordinate(1, "HO1") - coordinate(1, "O1")) == pytest.approx(
        0.96, abs=0.04
    )
    assert np.linalg.norm(coordinate(13, "H4O") - coordinate(13, "O4")) == pytest.approx(
        0.96, abs=0.04
    )


def test_constructed_finite_end_caps_have_no_subangstrom_nonbonded_clashes() -> None:
    mapped = map_crystal_to_finite_glycam(_dp12_slab(), REFERENCE_PDB)
    identities = list(
        zip(mapped.residue_numbers, mapped.atom_names, strict=True)
    )
    intended_terminal_bonds = {
        frozenset({(1, "HO1"), (1, "O1")}),
        frozenset({(1, "O1"), (2, "C1")}),
        frozenset({(13, "H4O"), (13, "O4")}),
    }
    terminal_atoms = {(1, "HO1"), (1, "O1"), (13, "H4O")}
    nonbonded_terminal_distances: list[float] = []
    for first, first_identity in enumerate(identities):
        if first_identity not in terminal_atoms:
            continue
        for second, second_identity in enumerate(identities):
            if first == second:
                continue
            if frozenset({first_identity, second_identity}) in intended_terminal_bonds:
                continue
            nonbonded_terminal_distances.append(
                float(
                    np.linalg.norm(
                        mapped.cartesian_coordinates[first]
                        - mapped.cartesian_coordinates[second]
                    )
                )
            )

    assert min(nonbonded_terminal_distances) > 0.8


def test_multiple_chains_repeat_the_same_255_atom_topology() -> None:
    mapped = map_crystal_to_finite_glycam(_dp12_slab(chains=4), REFERENCE_PDB)

    assert mapped.atom_count == 4 * 255
    assert mapped.chain_count == 4
    assert [mapped.chain_numbers.count(chain) for chain in (1, 2, 3, 4)] == [
        255,
        255,
        255,
        255,
    ]


def test_periodic_mapping_has_252_atoms_without_terminal_groups() -> None:
    mapped = map_crystal_to_periodic_glycam(
        _dp12_slab(chains=2, periodic=True), PERIODIC_REFERENCE_PDB
    )

    assert mapped.atom_count == 2 * 252
    assert mapped.chain_count == 2
    assert set(mapped.residue_names) == {"4GB"}
    assert set(mapped.residue_numbers) == set(range(1, 13))
    assert not {"O1", "HO1", "H4O"}.intersection(mapped.atom_names)


def test_iii_i_periodic_bulk_maps_to_complete_glycam_chain() -> None:
    crystal = build_bulk_crystal(load_allomorph("iii_i"), (1, 1, 6))

    mapped = map_crystal_to_periodic_glycam(crystal, PERIODIC_REFERENCE_PDB)

    assert mapped.chain_count == 1
    assert mapped.atom_count == 12 * 21
    assert set(mapped.residue_numbers) == set(range(1, 13))


def test_ii_mapping_generates_force_field_hydrogens_and_maps_antiparallel_chains() -> None:
    crystal = build_bulk_crystal(load_allomorph("ii"), (1, 1, 6))

    mapped = map_crystal_to_periodic_glycam(crystal, PERIODIC_REFERENCE_PDB)

    assert mapped.chain_count == 2
    assert mapped.atom_count == 2 * 12 * 21
    for chain_number in (1, 2):
        for residue_number in range(1, 13):
            residue = {
                name: coordinate
                for name, coordinate, chain, mapped_residue in zip(
                    mapped.atom_names,
                    mapped.cartesian_coordinates,
                    mapped.chain_numbers,
                    mapped.residue_numbers,
                    strict=True,
                )
                if chain == chain_number and mapped_residue == residue_number
            }
            for oxygen, hydrogen in (("O2", "H2O"), ("O3", "H3O"), ("O6", "H6O")):
                suffix = oxygen[1:]
                assert np.linalg.norm(residue[hydrogen] - residue[oxygen]) == pytest.approx(
                    0.96, abs=1.0e-6
                )
                assert _angle(residue[f"C{suffix}"], residue[oxygen], residue[hydrogen]) == (
                    pytest.approx(108.5, abs=1.0e-6)
                )
    inverse = np.linalg.inv(mapped.basis_vectors)
    closest_hydrogen = np.inf
    closest_heavy = np.inf
    for index, name in enumerate(mapped.atom_names):
        if name not in {"H2O", "H3O", "H6O"}:
            continue
        for other, other_name in enumerate(mapped.atom_names):
            if index == other:
                continue
            same_residue = (
                mapped.chain_numbers[index] == mapped.chain_numbers[other]
                and mapped.residue_numbers[index] == mapped.residue_numbers[other]
            )
            if same_residue and other_name == "O" + name[1]:
                continue
            delta = mapped.cartesian_coordinates[other] - mapped.cartesian_coordinates[index]
            fractional = delta @ inverse
            fractional -= np.round(fractional)
            distance = float(np.linalg.norm(fractional @ mapped.basis_vectors))
            if other_name.startswith("H"):
                closest_hydrogen = min(closest_hydrogen, distance)
            else:
                closest_heavy = min(closest_heavy, distance)
    assert closest_hydrogen >= 2.0
    assert closest_heavy >= 1.8
    for chain_number in (1, 2):
        chain_indices = [
            index for index, chain in enumerate(mapped.chain_numbers) if chain == chain_number
        ]
        for residue_number in range(1, 12):
            o4 = next(
                mapped.cartesian_coordinates[index]
                for index in chain_indices
                if mapped.residue_numbers[index] == residue_number
                and mapped.atom_names[index] == "O4"
            )
            c1 = next(
                mapped.cartesian_coordinates[index]
                for index in chain_indices
                if mapped.residue_numbers[index] == residue_number + 1
                and mapped.atom_names[index] == "C1"
            )
            displacement = c1 - o4
            fractional = np.linalg.solve(mapped.basis_vectors.T, displacement)
            fractional -= np.round(fractional)
            assert np.linalg.norm(fractional @ mapped.basis_vectors) < 1.7


def test_finite_ii_slab_maps_antiparallel_chains_and_constructs_hydrogens() -> None:
    crystal = build_chain_parallel_slab(
        load_allomorph("ii"),
        MillerIndex(1, 1, 0),
        layers=1,
        chains_per_layer=2,
        degree_polymerization=12,
        vacuum=30.0,
        padding_z=12.0,
        chain_boundary="finite",
    )

    mapped = map_crystal_to_finite_glycam(crystal, REFERENCE_PDB)

    assert mapped.chain_count == 2
    assert mapped.atom_count == 2 * 255
    for chain_number in (1, 2):
        for residue_number in range(2, 14):
            residue = {
                name: coordinate
                for name, coordinate, chain, mapped_residue in zip(
                    mapped.atom_names,
                    mapped.cartesian_coordinates,
                    mapped.chain_numbers,
                    mapped.residue_numbers,
                    strict=True,
                )
                if chain == chain_number and mapped_residue == residue_number
            }
            for oxygen, hydrogen in (("O2", "H2O"), ("O3", "H3O"), ("O6", "H6O")):
                assert np.linalg.norm(residue[hydrogen] - residue[oxygen]) == pytest.approx(
                    0.96,
                    abs=0.08,
                )
