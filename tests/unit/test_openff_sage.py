"""Regression tests for the optional experimental OpenFF Sage backend."""

from __future__ import annotations

import importlib
from collections import Counter
from importlib.util import find_spec

import numpy as np
import pytest

from cellulose_forge.crystallography import load_allomorph
from cellulose_forge.geometry import build_bulk_crystal
from cellulose_forge.openff_sage import (
    _REFERENCE_PARITY,
    _centre_parity,
    build_sage_molecule,
    build_sage_system,
    dispersion_correction_for_shape,
    reduced_periodic_box,
    require_openff,
)


def test_dispersion_policy_is_geometry_specific() -> None:
    assert dispersion_correction_for_shape("bulk") is True
    assert dispersion_correction_for_shape("crystal") is True
    assert dispersion_correction_for_shape("slab") is False
    with pytest.raises(ValueError, match="undefined"):
        dispersion_correction_for_shape("fibril")


def test_missing_optional_dependency_has_installation_message(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real_import = importlib.import_module

    def missing(name: str) -> object:
        if name == "openff.toolkit":
            raise ImportError("not installed")
        return real_import(name)

    monkeypatch.setattr(importlib, "import_module", missing)
    with pytest.raises(RuntimeError, match=r"environment-openff\.yml"):
        require_openff()


def test_triclinic_box_reduction_preserves_lattice_volume() -> None:
    box = np.asarray(((10.0, 0.0, 0.0), (-7.0, 9.0, 0.0), (6.0, 5.0, 8.0)))
    reduced = reduced_periodic_box(box)
    assert np.linalg.det(reduced) == pytest.approx(np.linalg.det(box))
    assert abs(reduced[1, 0]) <= reduced[0, 0] / 2
    assert abs(reduced[2, 0]) <= reduced[0, 0] / 2
    assert abs(reduced[2, 1]) <= reduced[1, 1] / 2


def _angle(left: np.ndarray, centre: np.ndarray, right: np.ndarray) -> float:
    first = left - centre
    second = right - centre
    cosine = np.dot(first, second) / (np.linalg.norm(first) * np.linalg.norm(second))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
def test_direct_graph_counts_and_periodic_scaling() -> None:
    toolkit = importlib.import_module("openff.toolkit")
    force_field_type, topology_type = toolkit.ForceField, toolkit.Topology
    force_field = force_field_type("openff-2.3.0.offxml")
    normalized: dict[int, dict[str, Counter[str]]] = {}
    for degree in (4, 6):
        molecule = build_sage_molecule(degree, "periodic")
        assert molecule.n_atoms == 21 * degree
        assert molecule.n_bonds == 22 * degree
        labels = force_field.label_molecules(topology_type.from_molecules([molecule]))[0]
        normalized[degree] = {}
        for handler in ("Bonds", "Angles", "ProperTorsions"):
            counts = Counter(parameter.id for parameter in labels[handler].values())
            assert all(count % degree == 0 for count in counts.values())
            normalized[degree][handler] = Counter(
                {parameter: count // degree for parameter, count in counts.items()}
            )
    assert normalized[4] == normalized[6]


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
@pytest.mark.parametrize(
    ("degree", "pubchem_cid", "isomeric_smiles"),
    (
        (
            2,
            10712,
            "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O[C@@H]2[C@H](O[C@H]"
            "([C@@H]([C@H]2O)O)O)CO)O)O)O)O",
        ),
        (
            6,
            440948,
            "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O[C@@H]2[C@H](O[C@H]"
            "([C@@H]([C@H]2O)O)O[C@@H]3[C@H](O[C@H]([C@@H]([C@H]3O)O)O"
            "[C@@H]4[C@H](O[C@H]([C@@H]([C@H]4O)O)O[C@@H]5[C@H](O[C@H]"
            "([C@@H]([C@H]5O)O)O[C@@H]6[C@H](O[C@H]([C@@H]([C@H]6O)O)O)"
            "CO)CO)CO)CO)CO)O)O)O)O",
        ),
    ),
)
def test_finite_graph_matches_pubchem_stereochemistry(
    degree: int, pubchem_cid: int, isomeric_smiles: str
) -> None:
    molecule_type = importlib.import_module("openff.toolkit").Molecule
    generated = build_sage_molecule(degree, "finite")
    reference = molecule_type.from_smiles(isomeric_smiles, allow_undefined_stereo=False)
    isomorphic, _ = molecule_type.are_isomorphic(
        generated,
        reference,
        atom_stereochemistry_matching=True,
        bond_stereochemistry_matching=True,
    )
    assert isomorphic, f"finite DP{degree} graph differs from PubChem CID {pubchem_cid}"


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
def test_cellulose_ii_modelled_hydrogens_are_tetrahedral() -> None:
    crystal = build_bulk_crystal(load_allomorph("ii"), (1, 1, 2))
    build = build_sage_system(crystal, "periodic")
    for chain in set(build.mapped.chain_numbers):
        for residue in set(build.mapped.residue_numbers):
            atoms = {
                name: coordinate
                for name, coordinate, atom_chain, atom_residue in zip(
                    build.mapped.atom_names,
                    build.mapped.cartesian_coordinates,
                    build.mapped.chain_numbers,
                    build.mapped.residue_numbers,
                    strict=True,
                )
                if atom_chain == chain and atom_residue == residue
            }
            for suffix in ("2", "3", "6"):
                assert np.linalg.norm(atoms[f"HO{suffix}"] - atoms[f"O{suffix}"]) == pytest.approx(
                    0.96, abs=1.0e-6
                )
                assert _angle(
                    atoms[f"C{suffix}"], atoms[f"O{suffix}"], atoms[f"HO{suffix}"]
                ) == pytest.approx(108.5, abs=1.0e-6)
    inverse = np.linalg.inv(build.mapped.basis_vectors)
    closest_hydrogen = np.inf
    closest_heavy = np.inf
    for index, name in enumerate(build.mapped.atom_names):
        if name not in {"HO2", "HO3", "HO6"}:
            continue
        for other, other_name in enumerate(build.mapped.atom_names):
            if index == other:
                continue
            same_residue = (
                build.mapped.chain_numbers[index] == build.mapped.chain_numbers[other]
                and build.mapped.residue_numbers[index] == build.mapped.residue_numbers[other]
            )
            if same_residue and other_name == "O" + name[2]:
                continue
            delta = (
                build.mapped.cartesian_coordinates[other]
                - build.mapped.cartesian_coordinates[index]
            )
            fractional = delta @ inverse
            fractional -= np.round(fractional)
            distance = float(np.linalg.norm(fractional @ build.mapped.basis_vectors))
            if other_name.startswith("H"):
                closest_hydrogen = min(closest_hydrogen, distance)
            else:
                closest_heavy = min(closest_heavy, distance)
    assert closest_hydrogen >= 2.0
    assert closest_heavy >= 1.8


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
def test_stereochemical_negative_controls_reject_inversions_and_mirror() -> None:
    molecule = build_sage_molecule(4, "periodic")
    # An idealized local tetrahedron is sufficient to prove that each control
    # reverses the independent signed-parity observable used at build time.
    positions = np.zeros((molecule.n_atoms, 3), dtype=float)
    for atom in molecule.atoms:
        if atom.name not in {"C1", "C3"} or int(atom.metadata["residue_number"]) != 2:
            continue
        centre = atom.molecule_atom_index
        positions[centre] = (4.0, 4.0, 4.0)
        neighbours = sorted(atom.bonded_atoms, key=lambda item: (item.symbol, item.name))[:3]
        offsets = np.eye(3)
        if int(np.sign(np.linalg.det(offsets))) != _REFERENCE_PARITY[atom.name]:
            offsets[[0, 1]] = offsets[[1, 0]]
        for neighbour, offset in zip(neighbours, offsets, strict=True):
            positions[neighbour.molecule_atom_index] = positions[centre] + offset
        original = _centre_parity(positions, atom, np.eye(3) * 20.0)
        inverted = positions.copy()
        left, right = (item.molecule_atom_index for item in neighbours[:2])
        inverted[[left, right]] = inverted[[right, left]]
        assert original == _REFERENCE_PARITY[atom.name]
        assert _centre_parity(inverted, atom, np.eye(3) * 20.0) == -original
        mirrored = positions.copy()
        mirrored[:, 0] *= -1
        assert _centre_parity(mirrored, atom, np.eye(3) * 20.0) == -original
