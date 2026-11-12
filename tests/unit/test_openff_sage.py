"""Regression tests for the optional experimental OpenFF Sage backend."""

from __future__ import annotations

import importlib
from collections import Counter
from importlib.util import find_spec

import numpy as np
import pytest

from cellulose_forge.openff_sage import (
    _REFERENCE_PARITY,
    _centre_parity,
    build_sage_molecule,
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
