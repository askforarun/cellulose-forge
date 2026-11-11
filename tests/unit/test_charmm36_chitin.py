"""CHARMM36 BGLCNA (chitin) topology built from the installed release."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pytest

from cellulose_forge.charmm36 import (
    _crystal_name,
    native_finite_charmm36_topology,
    render_finite_charmm36_itp,
)
from cellulose_forge.crystallography import load_allomorph
from cellulose_forge.forcefields import discover_charmm36_directory

_DP = 6


@pytest.fixture(scope="module")
def force_field() -> Path:
    try:
        return discover_charmm36_directory()
    except FileNotFoundError:
        pytest.skip("the pinned CHARMM36 release is not installed")


def _residue_summary(topology) -> dict[int, tuple[int, float]]:
    counts: dict[int, int] = defaultdict(int)
    charges: dict[int, float] = defaultdict(float)
    for atom in topology.topology_atoms:
        counts[atom.residue_number] += 1
        charges[atom.residue_number] += atom.charge
    return {residue: (counts[residue], round(charges[residue], 6)) for residue in counts}


def test_finite_chitin_dp6_atoms_charges_and_termini(force_field: Path) -> None:
    topology = native_finite_charmm36_topology(
        force_field, degree_polymerization=_DP, residue="BGLCNA"
    )
    assert topology.atom_count == 165
    assert abs(sum(atom.charge for atom in topology.topology_atoms)) < 1e-6
    summary = _residue_summary(topology)
    assert summary[1] == (29, -0.18)
    assert all(summary[residue] == (27, 0.0) for residue in range(2, _DP))
    assert summary[_DP] == (28, 0.18)
    names = defaultdict(set)
    for atom in topology.topology_atoms:
        names[atom.residue_number].add(atom.atom_name)
    assert {"O1", "HO1"} <= names[1] and "HO4" not in names[1]
    assert "HO4" in names[_DP] and not {"O1", "HO1"} & names[_DP]


def test_finite_chitin_renders_two_impropers_per_residue(force_field: Path) -> None:
    topology = native_finite_charmm36_topology(
        force_field, degree_polymerization=_DP, residue="BGLCNA"
    )
    impropers = [
        interaction
        for interaction in topology.interactions
        if interaction.section == "dihedrals" and interaction.parameters == ("2",)
    ]
    assert len(impropers) == 2 * _DP
    atoms = topology.topology_atoms
    assert {
        tuple(atoms[index].atom_name for index in interaction.atom_indices)
        for interaction in impropers
    } == {("C", "CT", "N", "O"), ("N", "C", "C2", "HN")}

    itp = render_finite_charmm36_itp(topology)
    improper_rows = [
        line.split()
        for line in itp.splitlines()
        if len(line.split()) == 5 and line.split()[-1] == "2"
    ]
    assert len(improper_rows) == 2 * _DP
    for atom_type in ("NC2D1", "CC2O1", "OC2D1", "CC331", "CC3161", "CC3162", "OC301"):
        assert f" {atom_type} " in itp


def test_bglcna_names_map_onto_chitin_crystal_labels(force_field: Path) -> None:
    topology = native_finite_charmm36_topology(
        force_field, degree_polymerization=_DP, residue="BGLCNA"
    )
    assert {
        name: _crystal_name(name) for name in ("N", "HN", "C", "O", "CT", "HT1", "HT2", "HT3")
    } == {
        "N": "N2",
        "HN": "HN2",
        "C": "C7",
        "O": "O7",
        "CT": "C8",
        "HT1": "H81",
        "HT2": "H82",
        "HT3": "H83",
    }
    internal = {
        _crystal_name(atom.atom_name)
        for atom in topology.topology_atoms
        if atom.residue_number == 2
    }
    for allomorph in ("alpha", "beta"):
        labels = set(load_allomorph(allomorph, polymer="chitin").asymmetric_unit.atom_labels)
        assert internal == labels


def test_default_residue_is_unchanged_bglc(force_field: Path) -> None:
    default = native_finite_charmm36_topology(force_field, degree_polymerization=_DP)
    explicit = native_finite_charmm36_topology(
        force_field, degree_polymerization=_DP, residue="BGLC"
    )
    assert render_finite_charmm36_itp(default) == render_finite_charmm36_itp(explicit)
    assert {atom.residue_name for atom in default.topology_atoms} == {"BGLC"}
    assert not any(interaction.parameters == ("2",) for interaction in default.interactions)
    assert _crystal_name("O2") == "O2"
