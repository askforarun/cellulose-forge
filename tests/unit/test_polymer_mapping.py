"""Contract tests for the shared chain-template mapper.

OpenFF Sage and OPLS-AA both map crystals through ``polymer_mapping``. These
tests pin the contract they rely on:

* every crystal atom is copied unchanged from the residue that an independent
  O4->C1 walk (implemented here, not imported) places at that template
  position;
* atoms come back in template order for every chain;
* OpenFF and OPLS receive identical coordinates for every shared atom;
* stored fingerprints, recorded from the pre-refactor mapper (commit 311be50),
  still match;
* malformed templates and crystals are refused with descriptive errors.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from importlib.util import find_spec
from pathlib import Path

import numpy as np
import pytest

from cellulose_forge.crystallography import MillerIndex, load_allomorph
from cellulose_forge.geometry import BulkCrystal, build_bulk_crystal, build_chain_parallel_slab
from cellulose_forge.glycam_mapping import GlycamMappedSystem
from cellulose_forge.opls import build_chain_graph, build_opls_system
from cellulose_forge.polymer_mapping import (
    ChainTemplateAtom,
    chain_template,
    map_crystal_to_chain_template,
)

_DP = 4
_FINGERPRINTS = json.loads(
    (Path(__file__).resolve().parents[1] / "data" / "mapping_fingerprints.json").read_text(
        encoding="utf-8"
    )
)
_PLACED = frozenset({"O1", "HO1", "HO4"})
_LINK_ANGSTROM = 2.0
_CELLULOSE_II_PLACED = frozenset({"HO2", "HO3", "HO6"})
_REQUIRES_OPENFF = pytest.mark.skipif(
    find_spec("openff") is None, reason="optional OpenFF environment"
)


def _slab(allomorph: str, boundary: str, polymer: str = "cellulose") -> BulkCrystal:
    return build_chain_parallel_slab(
        load_allomorph(allomorph, polymer),
        MillerIndex(1, 0, 0),
        layers=2,
        chains_per_layer=2,
        degree_polymerization=_DP,
        vacuum=20.0,
        chain_boundary=boundary,
    )


_CASES = {
    "cellulose_i_beta_slab_finite": (lambda: _slab("i_beta", "finite"), "finite"),
    "cellulose_i_beta_slab_periodic": (lambda: _slab("i_beta", "periodic"), "periodic"),
    "cellulose_ii_bulk_periodic": (
        lambda: build_bulk_crystal(load_allomorph("ii"), (2, 2, 2)),
        "periodic",
    ),
    "chitin_beta_slab_finite": (lambda: _slab("beta", "finite", "chitin"), "finite"),
    "chitin_beta_slab_periodic": (lambda: _slab("beta", "periodic", "chitin"), "periodic"),
}


_CELLULOSE_CASES = sorted(case for case in _CASES if case.startswith("cellulose"))


def _opls_template(boundary: str) -> tuple[ChainTemplateAtom, ...]:
    graph = build_chain_graph(_DP, boundary)
    return chain_template((atom.name, atom.residue) for atom in graph.atoms)


def _fingerprint(mapped: GlycamMappedSystem) -> str:
    """Hash identities and coordinates rounded to 1e-6 A (platform-stable)."""
    digest = hashlib.sha256()
    for chain, residue, name, xyz in zip(
        mapped.chain_numbers,
        mapped.residue_numbers,
        mapped.atom_names,
        mapped.cartesian_coordinates,
        strict=True,
    ):
        digest.update(f"{chain} {residue} {name} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}\n".encode())
    return digest.hexdigest()


def _by_key(mapped: GlycamMappedSystem) -> dict[tuple[int, int, str], tuple[float, ...]]:
    return {
        (chain, residue, name): tuple(float(value) for value in xyz)
        for chain, residue, name, xyz in zip(
            mapped.chain_numbers,
            mapped.residue_numbers,
            mapped.atom_names,
            mapped.cartesian_coordinates,
            strict=True,
        )
    }


def _crystal_atoms(crystal: BulkCrystal) -> dict[tuple[int, int, str], np.ndarray]:
    return {
        (chain, residue, name): np.asarray(xyz)
        for chain, residue, name, xyz in zip(
            crystal.chain_numbers,
            crystal.residue_numbers,
            crystal.atom_labels,
            crystal.cartesian_coordinates,
            strict=True,
        )
    }


def _walk(crystal: BulkCrystal, chain: int, periodic: bool) -> list[int]:
    """Independently order one chain's residues along O4->C1 links."""
    atoms = _crystal_atoms(crystal)
    residues = sorted({residue for c, residue, _ in atoms if c == chain})
    inverse = np.linalg.inv(crystal.basis_vectors)

    def separation(donor: int, acceptor: int) -> float:
        delta = atoms[chain, acceptor, "C1"] - atoms[chain, donor, "O4"]
        if periodic:
            fractional = delta @ inverse
            delta = (fractional - np.round(fractional)) @ crystal.basis_vectors
        return float(np.linalg.norm(delta))

    following: dict[int, int] = {}
    for donor in residues:
        linked = [
            acceptor
            for acceptor in residues
            if acceptor != donor and separation(donor, acceptor) < _LINK_ANGSTROM
        ]
        if linked:
            following[donor] = linked[0]
    if periodic:
        order = [residues[0]]
    else:
        order = [next(residue for residue in residues if residue not in following.values())]
    while len(order) < len(residues):
        order.append(following[order[-1]])
    return order


def _assert_copied_from_walked_residues(
    crystal: BulkCrystal, mapped: GlycamMappedSystem, boundary: str
) -> None:
    atoms = _crystal_atoms(crystal)
    placed = _PLACED | (_CELLULOSE_II_PLACED if crystal.allomorph == "II" else frozenset())
    orders = {
        chain: _walk(crystal, chain, boundary == "periodic")
        for chain in dict.fromkeys(crystal.chain_numbers)
    }
    copied = 0
    for (chain, template_residue, name), xyz in _by_key(mapped).items():
        if name in placed:
            continue
        crystal_residue = orders[chain][template_residue - 1]
        assert xyz == tuple(atoms[chain, crystal_residue, name]), (chain, template_residue, name)
        copied += 1
    assert copied > 0


@pytest.mark.parametrize("case", _CELLULOSE_CASES)
def test_opls_atoms_are_copied_from_the_walked_crystal_residue(case: str) -> None:
    make, boundary = _CASES[case]
    crystal = make()
    mapped = build_opls_system(crystal, "opls-aa", boundary).mapped
    _assert_copied_from_walked_residues(crystal, mapped, boundary)


@_REQUIRES_OPENFF
@pytest.mark.parametrize("case", ["chitin_beta_slab_finite", "chitin_beta_slab_periodic"])
def test_sage_chitin_atoms_are_copied_from_the_walked_crystal_residue(case: str) -> None:
    from cellulose_forge.openff_sage import build_sage_system  # noqa: PLC0415

    make, boundary = _CASES[case]
    crystal = make()
    mapped = build_sage_system(crystal, boundary, "chitin").mapped
    _assert_copied_from_walked_residues(crystal, mapped, boundary)
    assert set(mapped.residue_names) == {"NAG"}


@pytest.mark.parametrize("boundary", ["finite", "periodic"])
def test_atoms_follow_template_order_for_every_chain(boundary: str) -> None:
    crystal = _slab("i_beta", boundary)
    template = _opls_template(boundary)
    mapped = map_crystal_to_chain_template(crystal, template, boundary, "BGC")
    chains = tuple(dict.fromkeys(crystal.chain_numbers))
    assert len(mapped.atom_names) == len(chains) * len(template)
    for index, chain in enumerate(chains):
        block = slice(index * len(template), (index + 1) * len(template))
        assert set(mapped.chain_numbers[block]) == {chain}
        assert mapped.atom_names[block] == tuple(atom.name for atom in template)
        assert mapped.residue_numbers[block] == tuple(atom.residue_number for atom in template)
    assert set(mapped.residue_names) == {"BGC"}


def test_finite_caps_bond_to_their_anchor_atoms() -> None:
    mapped = _by_key(
        map_crystal_to_chain_template(
            _slab("i_beta", "finite"), _opls_template("finite"), "finite", "BGC"
        )
    )
    for chain in {key[0] for key in mapped}:
        for donor, acceptor, residue in (("C1", "O1", 1), ("O1", "HO1", 1), ("O4", "HO4", _DP)):
            distance = np.linalg.norm(
                np.subtract(mapped[chain, residue, donor], mapped[chain, residue, acceptor])
            )
            assert 0.9 < distance < 1.5, (chain, donor, acceptor, distance)


def test_template_order_does_not_change_coordinates() -> None:
    crystal = _slab("i_beta", "finite")
    template = _opls_template("finite")
    forward = map_crystal_to_chain_template(crystal, template, "finite", "BGC")
    reverse = map_crystal_to_chain_template(crystal, template[::-1], "finite", "BGC")
    assert _by_key(forward) == _by_key(reverse)


@pytest.mark.parametrize("case", _CELLULOSE_CASES)
def test_opls_mapping_matches_stored_pre_refactor_fingerprint(case: str) -> None:
    make, boundary = _CASES[case]
    mapped = build_opls_system(make(), "opls-aa", boundary).mapped
    assert _fingerprint(mapped) == _FINGERPRINTS["cases"][case]["opls"]


@_REQUIRES_OPENFF
@pytest.mark.parametrize("case", sorted(_CASES))
def test_sage_mapping_matches_stored_pre_refactor_fingerprint(case: str) -> None:
    from cellulose_forge.openff_sage import build_sage_system  # noqa: PLC0415

    make, boundary = _CASES[case]
    polymer = "chitin" if case.startswith("chitin") else "cellulose"
    mapped = build_sage_system(make(), boundary, polymer).mapped
    assert _fingerprint(mapped) == _FINGERPRINTS["cases"][case]["sage"]


@_REQUIRES_OPENFF
@pytest.mark.parametrize("boundary", ["finite", "periodic"])
@pytest.mark.parametrize("allomorph", ["i_beta", "ii"])
def test_openff_and_opls_receive_identical_coordinates(allomorph: str, boundary: str) -> None:
    from cellulose_forge.openff_sage import build_sage_system  # noqa: PLC0415

    crystal = _slab(allomorph, boundary)
    sage = _by_key(build_sage_system(crystal, boundary, "cellulose").mapped)
    opls = _by_key(build_opls_system(crystal, "opls-aa", boundary).mapped)
    assert sage == opls


def test_opls_refuses_non_cellulose_crystals() -> None:
    mannan = build_bulk_crystal(load_allomorph("i", "mannan"), (1, 1, 2))
    with pytest.raises(ValueError, match=r"cellulose \(BGC\) residues only; this crystal has MAN"):
        build_opls_system(mannan, "opls-aa", "periodic")


def _without(crystal: BulkCrystal, chain: int, residue: int, name: str) -> BulkCrystal:
    keep = [
        index
        for index, key in enumerate(
            zip(crystal.chain_numbers, crystal.residue_numbers, crystal.atom_labels, strict=True)
        )
        if key != (chain, residue, name)
    ]

    def pick(values: tuple[object, ...]) -> tuple[object, ...]:
        return tuple(values[index] for index in keep) if values else values

    return replace(
        crystal,
        atom_labels=pick(crystal.atom_labels),
        cartesian_coordinates=crystal.cartesian_coordinates[keep],
        chain_ids=pick(crystal.chain_ids),
        residue_numbers=pick(crystal.residue_numbers),
        chain_numbers=pick(crystal.chain_numbers),
        layer_numbers=pick(crystal.layer_numbers),
    )


@pytest.mark.parametrize("name", ["C1", "O4"])
def test_missing_link_atom_is_reported_by_chain_and_residue(name: str) -> None:
    crystal = _slab("i_beta", "finite")
    chain = crystal.chain_numbers[0]
    residue = crystal.residue_numbers[0]
    with pytest.raises(
        ValueError, match=rf"chain {chain} residue {residue} lacks {name}; the glycosidic"
    ):
        map_crystal_to_chain_template(
            _without(crystal, chain, residue, name), _opls_template("finite"), "finite", "BGC"
        )


def test_second_c1_within_link_distance_is_refused_as_ambiguous() -> None:
    crystal = _slab("i_beta", "finite")
    chain = crystal.chain_numbers[0]
    order = _walk(crystal, chain, periodic=False)
    atoms = _crystal_atoms(crystal)
    donor, intruder = order[0], order[2]
    target = atoms[chain, donor, "O4"] + np.array((0.0, 0.0, 1.6))
    coordinates = np.array(crystal.cartesian_coordinates)
    index = next(
        index
        for index, key in enumerate(
            zip(crystal.chain_numbers, crystal.residue_numbers, crystal.atom_labels, strict=True)
        )
        if key == (chain, intruder, "C1")
    )
    coordinates[index] = target
    with pytest.raises(
        ValueError,
        match=rf"chain {chain} residue {donor} O4 has more than one C1 within 2\.0 A",
    ):
        map_crystal_to_chain_template(
            replace(crystal, cartesian_coordinates=coordinates),
            _opls_template("finite"),
            "finite",
            "BGC",
        )


def test_rejects_duplicate_template_atoms() -> None:
    with pytest.raises(ValueError, match="duplicate"):
        chain_template([("C1", 1), ("C1", 1)])


def test_rejects_empty_template() -> None:
    with pytest.raises(ValueError, match="at least one atom"):
        chain_template([])


def test_rejects_terminal_caps_on_periodic_chains() -> None:
    crystal = _slab("i_beta", "periodic")
    with pytest.raises(ValueError, match="terminal cap"):
        map_crystal_to_chain_template(crystal, _opls_template("finite"), "periodic", "BGC")


def test_rejects_degree_of_polymerization_mismatch() -> None:
    crystal = _slab("i_beta", "periodic")
    short = tuple(atom for atom in _opls_template("periodic") if atom.residue_number < _DP)
    with pytest.raises(ValueError, match="residues but the template has"):
        map_crystal_to_chain_template(crystal, short, "periodic", "BGC")


def test_rejects_template_atoms_absent_from_crystal() -> None:
    crystal = _slab("i_beta", "periodic")
    template = (*_opls_template("periodic"), ChainTemplateAtom("N2", 1))
    with pytest.raises(ValueError, match="absent from crystal"):
        map_crystal_to_chain_template(crystal, template, "periodic", "BGC")


def test_rejects_unknown_chain_boundary() -> None:
    crystal = _slab("i_beta", "finite")
    with pytest.raises(ValueError, match="finite or periodic"):
        map_crystal_to_chain_template(crystal, _opls_template("finite"), "cyclic", "BGC")
