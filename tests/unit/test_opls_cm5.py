"""Tests for the experimental OPLS-CM5 reference adapter."""

from pathlib import Path

import numpy as np
import pytest

from cellulose_forge.crystallography import load_allomorph
from cellulose_forge.geometry import build_bulk_crystal
from cellulose_forge.opls_cm5 import (
    PeriodicOplsTopology,
    _InteractionTemplate,
    _read_first_reference_chain,
    _read_topology_atoms,
    _TopologyAtom,
    derive_periodic_opls_cm5_topology,
    map_periodic_opls_cm5,
    periodic_interaction_counts,
    render_periodic_opls_cm5_itp,
)


def test_opls_reference_topology_must_be_neutral(tmp_path: Path) -> None:
    topology = tmp_path / "CNC.itp"
    rows = [
        f"{index} opls_140 1 CNC HC {index} {0.1 if index == 1 else 0.0} 1.008"
        for index in range(1, 424)
    ]
    topology.write_text("[ atoms ]\n" + "\n".join(rows) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="not neutral"):
        _read_topology_atoms(topology)


def test_opls_reference_gro_requires_expected_dp20_layout(tmp_path: Path) -> None:
    reference = tmp_path / "reference.gro"
    residue_counts = [22, *([21] * 18), 23]
    lines = ["reference", "423"]
    atom_index = 1
    for residue_number, count in enumerate(residue_counts, start=1):
        for _ in range(count):
            lines.append(
                f"{residue_number:5d}{'BGLC':<5}{'C1':>5}{atom_index:5d}"
                f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}"
            )
            atom_index += 1
    lines.append("1.0 1.0 1.0")
    reference.write_text("\n".join(lines) + "\n", encoding="utf-8")

    atoms = _read_first_reference_chain(reference)

    assert len(atoms) == 423
    assert atoms[0].residue_number == 1
    assert atoms[-1].residue_number == 20


def test_periodic_opls_requires_even_degree_polymerization() -> None:
    with pytest.raises(ValueError, match="even degree of polymerization"):
        derive_periodic_opls_cm5_topology("missing.itp", "missing.gro", degree_polymerization=3)


def test_periodic_opls_render_wraps_interactions_and_neutralizes_charge() -> None:
    topology = PeriodicOplsTopology(
        degree_polymerization=4,
        coordinate_atom_names=("A", "B"),
        reference_coordinates=np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        topology_atoms=(
            _TopologyAtom("opls_1", "A", 0.1, 12.0),
            _TopologyAtom("opls_2", "B", 0.1, 1.0),
        ),
        interactions=tuple(
            _InteractionTemplate(section, atoms, parameters)
            for section, atoms, parameters in (
                ("bonds", ((0, "A"), (1, "B")), ("1", "0.1", "100")),
                ("angles", ((0, "A"), (0, "B"), (1, "A")), ("1", "100", "10")),
                (
                    "dihedrals",
                    ((0, "A"), (0, "B"), (1, "A"), (1, "B")),
                    ("3", "1", "2", "3", "4", "5", "6"),
                ),
                ("pairs", ((0, "A"), (1, "B")), ("1",)),
            )
        ),
        source_repeat_charge=0.2,
    )

    rendered = render_periodic_opls_cm5_itp(topology)

    assert topology.corrected_repeat_charge == pytest.approx(0.0)
    assert periodic_interaction_counts(topology) == {
        "bonds": 4,
        "angles": 4,
        "dihedrals": 4,
        "pairs": 4,
    }
    assert "7 2 1 0.1 100" in rendered
    atom_rows = [line.split() for line in rendered.splitlines() if "opls_" in line]
    assert sum(float(row[6]) for row in atom_rows) == pytest.approx(0.0)


def test_periodic_opls_maps_cellulose_ii_antiparallel_chains() -> None:
    source = build_bulk_crystal(load_allomorph("i_beta"), (1, 1, 1))
    source_indices = [
        index
        for index, (chain, residue) in enumerate(
            zip(source.chain_numbers, source.residue_numbers, strict=True)
        )
        if chain == 1 and residue == 1
    ]
    coordinate_names = tuple(source.atom_labels[index] for index in source_indices)
    reference_coordinates = source.cartesian_coordinates[source_indices]
    topology = PeriodicOplsTopology(
        degree_polymerization=12,
        coordinate_atom_names=coordinate_names,
        reference_coordinates=reference_coordinates,
        topology_atoms=tuple(
            _TopologyAtom("opls_1", name, 0.0, 1.0) for name in coordinate_names
        ),
        interactions=(),
        source_repeat_charge=0.0,
    )
    crystal = build_bulk_crystal(load_allomorph("ii"), (1, 1, 6))

    mapped = map_periodic_opls_cm5(crystal, topology)

    assert mapped.chain_count == 2
    for chain_number in (1, 2):
        for residue_number in range(1, 13):
            global_residue = (chain_number - 1) * 12 + residue_number
            residue = {
                name: coordinate
                for name, coordinate, chain, mapped_residue in zip(
                    mapped.atom_names,
                    mapped.cartesian_coordinates,
                    mapped.chain_numbers,
                    mapped.residue_numbers,
                    strict=True,
                )
                if chain == chain_number and mapped_residue == global_residue
            }
            for oxygen, hydrogen in (("O2", "HO2"), ("O3", "HO3"), ("O6", "HO6")):
                assert 0.8 < np.linalg.norm(residue[hydrogen] - residue[oxygen]) < 1.1
        for residue_number in range(1, 12):
            global_residue = (chain_number - 1) * 12 + residue_number
            c1 = next(
                coordinate
                for name, coordinate, chain, mapped_residue in zip(
                    mapped.atom_names,
                    mapped.cartesian_coordinates,
                    mapped.chain_numbers,
                    mapped.residue_numbers,
                    strict=True,
                )
                if chain == chain_number
                and mapped_residue == global_residue
                and name == "C1"
            )
            o4 = next(
                coordinate
                for name, coordinate, chain, mapped_residue in zip(
                    mapped.atom_names,
                    mapped.cartesian_coordinates,
                    mapped.chain_numbers,
                    mapped.residue_numbers,
                    strict=True,
                )
                if chain == chain_number
                and mapped_residue == global_residue + 1
                and name == "O4"
            )
            displacement = o4 - c1
            fractional = np.linalg.solve(mapped.basis_vectors.T, displacement)
            fractional -= np.round(fractional)
            assert np.linalg.norm(fractional @ mapped.basis_vectors) < 1.7
