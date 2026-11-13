"""Regression tests for the optional experimental OpenFF Sage backend."""

from __future__ import annotations

import importlib
import time
from collections import Counter
from importlib.util import find_spec
from pathlib import Path

import numpy as np
import pytest
from hydroxyl_contacts import modelled_hydroxyl_contacts

from cellulose_forge import cli
from cellulose_forge.cli import main
from cellulose_forge.crystallography import load_allomorph
from cellulose_forge.geometry import build_bulk_crystal
from cellulose_forge.glycam_mapping import map_crystal_to_periodic_glycam
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
    assert dispersion_correction_for_shape("fibril") is False
    with pytest.raises(ValueError, match="undefined"):
        dispersion_correction_for_shape("cylinder")


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
    contacts = modelled_hydroxyl_contacts(
        build.mapped, {"HO2": "O2", "HO3": "O3", "HO6": "O6"}
    )
    assert contacts["hydrogen"] >= 1.9
    # Hydrogen-bond acceptors may approach to ~1.7 A, as in neutron structures.
    assert contacts["oxygen"] >= 1.6
    assert contacts["other_heavy"] >= 2.0
    assert contacts["hydrogen_bonded_fraction"] >= 0.95
    # The intrachain O3-H...O5 bond is present in every cellulose allomorph.
    assert contacts["o3_to_intrachain_o5_fraction"] == 1.0


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
def test_cellulose_ii_hydrogen_placement_scales_to_large_bulk() -> None:
    crystal = build_bulk_crystal(load_allomorph("ii"), (6, 5, 4))
    start = time.perf_counter()
    build = build_sage_system(crystal, "periodic")
    assert time.perf_counter() - start < 60.0
    contacts = modelled_hydroxyl_contacts(
        build.mapped, {"HO2": "O2", "HO3": "O3", "HO6": "O6"}
    )
    assert contacts["hydrogen_bonded_fraction"] >= 0.95
    # The intrachain O3-H...O5 bond is present in every cellulose allomorph.
    assert contacts["o3_to_intrachain_o5_fraction"] == 1.0


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


@pytest.mark.parametrize(
    "arguments",
    (
        ["--shape", "slab", "--face", "1", "1", "0", "--layers", "6", "--chains-per-layer", "4"],
        ["--shape", "fibril", "--fibril-model", "18-chain"],
    ),
)
def test_openff_vacuum_shapes_warn_below_cutoff_vacuum(
    tmp_path: Path, arguments: list[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    def stop(*_: object, **__: object) -> None:
        raise SystemExit(0)

    # Stop before any OpenFF work; only the CLI vacuum policy is under test.
    monkeypatch.setattr(cli, "build_sage_system", stop)
    with pytest.warns(UserWarning, match="less than 9 angstrom of transverse vacuum"):
        with pytest.raises(SystemExit):
            main(
                [
                    "build", "--allomorph", "I-beta", *arguments,
                    "--degree-polymerization", "6", "--chain-boundary", "periodic",
                    "--vacuum", "5", "--force-field", "openff", "--md-engine", "gromacs",
                    "--output", str(tmp_path / "out"),
                ]
            )


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
def test_glycam_and_openff_place_identical_cellulose_ii_hydrogens() -> None:
    crystal = build_bulk_crystal(load_allomorph("ii"), (2, 2, 6))
    glycam = map_crystal_to_periodic_glycam(
        crystal, Path("examples/glycam/dp12_periodic_reference/cellulose_periodic_chain.pdb")
    )
    openff = build_sage_system(crystal, "periodic").mapped

    def hydrogens_by_oxygen(
        mapped: object, names: dict[str, str]
    ) -> dict[tuple[float, float, float], np.ndarray]:
        atoms = list(
            zip(
                mapped.atom_names,  # type: ignore[attr-defined]
                mapped.chain_numbers,  # type: ignore[attr-defined]
                mapped.residue_numbers,  # type: ignore[attr-defined]
                np.asarray(mapped.cartesian_coordinates),  # type: ignore[attr-defined]
                strict=True,
            )
        )
        oxygen = {(name, chain, residue): xyz for name, chain, residue, xyz in atoms}
        return {
            tuple(np.round(oxygen[(names[name], chain, residue)], 4)): xyz
            for name, chain, residue, xyz in atoms
            if name in names
        }

    from_glycam = hydrogens_by_oxygen(glycam, {"H2O": "O2", "H3O": "O3", "H6O": "O6"})
    from_openff = hydrogens_by_oxygen(openff, {"HO2": "O2", "HO3": "O3", "HO6": "O6"})
    assert from_glycam.keys() == from_openff.keys()
    for oxygen_position, hydrogen in from_glycam.items():
        np.testing.assert_allclose(hydrogen, from_openff[oxygen_position], atol=1e-6)
