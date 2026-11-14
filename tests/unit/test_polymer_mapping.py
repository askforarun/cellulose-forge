"""Contract tests for the force-field-neutral chain-template mapper.

OpenFF Sage and OPLS-AA both map crystals through ``polymer_mapping``. These
tests pin the contract they rely on: atoms come back in template order, every
shared (chain, residue, atom) receives the same coordinate regardless of the
backend that requested it, and malformed templates are refused.
"""

from __future__ import annotations

from importlib.util import find_spec

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


def _opls_template(boundary: str) -> tuple[ChainTemplateAtom, ...]:
    graph = build_chain_graph(_DP, boundary)
    return chain_template((atom.name, atom.residue) for atom in graph.atoms)


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


@pytest.mark.parametrize("boundary", ["finite", "periodic"])
@pytest.mark.parametrize("allomorph", ["i_beta", "ii"])
def test_atoms_follow_template_order_for_every_chain(allomorph: str, boundary: str) -> None:
    crystal = _slab(allomorph, boundary)
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
    assert np.isfinite(mapped.cartesian_coordinates).all()


@pytest.mark.parametrize("boundary", ["finite", "periodic"])
def test_crystal_atoms_are_copied_not_moved(boundary: str) -> None:
    crystal = _slab("i_beta", boundary)
    mapped = map_crystal_to_chain_template(crystal, _opls_template(boundary), boundary, "BGC")
    crystal_xyz = {
        (chain, name): set()
        for chain, name in zip(crystal.chain_numbers, crystal.atom_labels, strict=True)
    }
    for chain, name, xyz in zip(
        crystal.chain_numbers, crystal.atom_labels, crystal.cartesian_coordinates, strict=True
    ):
        crystal_xyz[chain, name].add(tuple(float(value) for value in xyz))
    for (chain, _, name), xyz in _by_key(mapped).items():
        if name not in {"O1", "HO1", "HO4"}:
            assert xyz in crystal_xyz[chain, name]


def test_template_order_does_not_change_coordinates() -> None:
    crystal = _slab("i_beta", "finite")
    template = _opls_template("finite")
    forward = map_crystal_to_chain_template(crystal, template, "finite", "BGC")
    reverse = map_crystal_to_chain_template(crystal, template[::-1], "finite", "BGC")
    assert _by_key(forward) == _by_key(reverse)


def test_opls_build_uses_the_shared_mapper() -> None:
    crystal = build_bulk_crystal(load_allomorph("i_beta"), (2, 2, _DP // 2))
    template = _opls_template("periodic")
    expected = map_crystal_to_chain_template(crystal, template, "periodic", "BGC")
    built = build_opls_system(crystal, "opls-aa", "periodic").mapped
    assert built.atom_names == expected.atom_names
    assert _by_key(built) == _by_key(expected)


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
@pytest.mark.parametrize("boundary", ["finite", "periodic"])
@pytest.mark.parametrize("allomorph", ["i_beta", "ii"])
def test_openff_and_opls_receive_identical_coordinates(allomorph: str, boundary: str) -> None:
    from cellulose_forge.openff_sage import build_sage_system  # noqa: PLC0415

    crystal = _slab(allomorph, boundary)
    sage = _by_key(build_sage_system(crystal, boundary, "cellulose").mapped)
    opls = _by_key(build_opls_system(crystal, "opls-aa", boundary).mapped)
    assert sage.keys() == opls.keys()
    assert sage == opls


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
@pytest.mark.parametrize("boundary", ["finite", "periodic"])
def test_openff_chitin_mapping_matches_direct_template_call(boundary: str) -> None:
    from cellulose_forge.openff_sage import build_sage_system  # noqa: PLC0415

    crystal = _slab("beta", boundary, "chitin")
    build = build_sage_system(crystal, boundary, "chitin")
    template = chain_template(
        (atom.name, int(atom.metadata["residue_number"])) for atom in build.molecule.atoms
    )
    direct = map_crystal_to_chain_template(crystal, template, boundary, "NAG")
    assert build.mapped.atom_names == direct.atom_names
    assert _by_key(build.mapped) == _by_key(direct)
    assert set(build.mapped.residue_names) == {"NAG"}


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
