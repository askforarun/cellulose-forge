"""Tests for chain-level position-restraint selection and ordering."""

import numpy as np
import pytest

from cellulose_forge.crystallography import MillerIndex, load_allomorph
from cellulose_forge.geometry import build_36_chain_fibril, build_chain_parallel_slab
from cellulose_forge.glycam_mapping import GlycamMappedSystem
from cellulose_forge.restraints import (
    canonical_atom_name,
    group_mapped_chains,
    position_restraints_section,
    resolve_chain_restraints,
    resolve_mapped_restraint_atoms,
    transverse_chain_layers,
)


def test_position_restraints_include_only_heavy_atoms() -> None:
    section = position_restraints_section(
        ("C1", "H1", "O5", "HO5"),
        (1000.0, 500.0, 0.0),
    )

    assert "     1      1   1000.000    500.000      0.000" in section
    assert "     3      1   1000.000    500.000      0.000" in section
    assert "     2      1" not in section
    assert "     4      1" not in section


def test_position_restraints_filter_canonical_atom_names() -> None:
    section = position_restraints_section(
        ("C1", "H1", "C6", "O6", "O5"),
        (1000.0, 1000.0, 1000.0),
        selected_atom_names=("C1", "O5"),
    )

    assert "     1      1" in section
    assert "     5      1" in section
    assert "     3      1" not in section
    assert "     4      1" not in section


def test_shared_selector_uses_canonical_chitin_names_and_checks_each_chain() -> None:
    system = GlycamMappedSystem(
        atom_names=("C", "CT", "HN", "C", "CT", "HN"),
        residue_names=("BGLCNA",) * 6,
        residue_numbers=(1,) * 6,
        chain_numbers=(1, 1, 1, 2, 2, 2),
        cartesian_coordinates=np.zeros((6, 3)),
        basis_vectors=np.eye(3),
    )

    selected = resolve_mapped_restraint_atoms(
        system,
        (1, 2),
        requested_atom_names=("C7",),
        polymer="chitin",
        naming="charmm",
    )

    assert selected.atom_indices == (0, 3)
    assert selected.atoms_per_chain == ((1, 1), (2, 1))
    assert canonical_atom_name("CT", polymer="chitin", naming="charmm") == "C8"


def test_shared_selector_rejects_a_selected_chain_with_no_matching_atom() -> None:
    system = GlycamMappedSystem(
        atom_names=("O1", "C1"),
        residue_names=("BGC", "BGC"),
        residue_numbers=(1, 1),
        chain_numbers=(1, 2),
        cartesian_coordinates=np.zeros((2, 3)),
        basis_vectors=np.eye(3),
    )

    with pytest.raises(ValueError, match="chain IDs: 2"):
        resolve_mapped_restraint_atoms(
            system,
            (1, 2),
            requested_atom_names=("O1",),
        )


def test_resolves_universal_chain_ids() -> None:
    crystal = build_36_chain_fibril(
        load_allomorph("i_beta"),
        degree_polymerization=4,
        vacuum=0.0,
        chain_boundary="periodic",
    )

    restraints = resolve_chain_restraints(
        crystal,
        chain_ids=[31, 6, 11, 6],
        force_constants=[1000.0, 500.0, 0.0],
    )

    assert restraints is not None
    assert restraints.chain_ids == (6, 11, 31)
    assert restraints.force_constants == (1000.0, 500.0, 0.0)
    assert restraints.selection_mode == "all_heavy_atoms"


def test_rejects_unknown_chain_id() -> None:
    crystal = build_36_chain_fibril(load_allomorph("i_beta"), degree_polymerization=4, vacuum=0.0)

    with pytest.raises(ValueError, match="unknown chain IDs: 37"):
        resolve_chain_restraints(
            crystal,
            chain_ids=[37],
            force_constants=[1000.0, 1000.0, 1000.0],
        )


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_rejects_nonfinite_force_constants(value: float) -> None:
    crystal = build_36_chain_fibril(load_allomorph("i_beta"), degree_polymerization=4, vacuum=0.0)
    with pytest.raises(ValueError, match="must be finite"):
        resolve_chain_restraints(
            crystal,
            chain_ids=[1],
            force_constants=[1000.0, value, 0.0],
        )


def test_resolves_layers() -> None:
    crystal = build_chain_parallel_slab(
        load_allomorph("i_beta"),
        face=MillerIndex(1, 1, 0),
        layers=3,
        chains_per_layer=4,
        degree_polymerization=4,
        vacuum=0.0,
        chain_boundary="finite",
    )

    restraints = resolve_chain_restraints(
        crystal,
        chain_ids=None,
        layer_ids=[2],
        force_constants=[1000.0, 1000.0, 1000.0],
    )

    assert restraints is not None
    assert restraints.layer_ids == (2,)
    assert restraints.layer_mapping == (
        (1, (1, 2, 3, 4)),
        (2, (5, 6, 7, 8)),
        (3, (9, 10, 11, 12)),
    )
    assert restraints.chain_ids == (5, 6, 7, 8)


def test_resolves_named_ring_atom_filter() -> None:
    crystal = build_36_chain_fibril(load_allomorph("i_beta"), degree_polymerization=4, vacuum=0.0)
    restraints = resolve_chain_restraints(
        crystal,
        chain_ids=[1],
        force_constants=[1000.0, 1000.0, 1000.0],
        atom_names=["c1", "C2", "C1", "C3", "C4", "C5", "O5"],
    )

    assert restraints is not None
    assert restraints.atom_names == ("C1", "C2", "C3", "C4", "C5", "O5")
    assert restraints.selection_mode == "named_heavy_atoms"


def test_rejects_atom_names_for_wrong_polymer() -> None:
    crystal = build_36_chain_fibril(load_allomorph("i_beta"), degree_polymerization=4, vacuum=0.0)
    with pytest.raises(ValueError, match=r"unknown cellulose.*N2"):
        resolve_chain_restraints(
            crystal,
            chain_ids=[1],
            force_constants=[1000.0, 1000.0, 1000.0],
            atom_names=["N2"],
            polymer="cellulose",
        )


def test_rejects_atom_names_without_chain_or_layer_scope() -> None:
    crystal = build_36_chain_fibril(load_allomorph("i_beta"), degree_polymerization=4, vacuum=0.0)
    with pytest.raises(ValueError, match="require a chain-restraint selector"):
        resolve_chain_restraints(
            crystal,
            chain_ids=None,
            atom_names=["C1"],
        )


def test_rejects_combined_chain_and_layer_selectors() -> None:
    crystal = build_36_chain_fibril(load_allomorph("i_beta"), degree_polymerization=4, vacuum=0.0)

    with pytest.raises(ValueError, match="mutually exclusive"):
        resolve_chain_restraints(
            crystal,
            chain_ids=[1],
            layer_ids=[1],
            force_constants=[1000.0, 1000.0, 1000.0],
        )


def test_fibril_layers_progress_from_low_to_high_z_after_rotation() -> None:
    crystal = build_36_chain_fibril(
        load_allomorph("i_beta"),
        degree_polymerization=4,
        vacuum=0.0,
        chain_boundary="periodic",
    )

    assert transverse_chain_layers(crystal) == (
        (1, (1, 2, 3)),
        (2, (4, 5, 6, 7)),
        (3, (8, 9, 10, 11, 12)),
        (4, (13, 14, 15, 16, 17, 18)),
        (5, (19, 20, 21, 22, 23, 24)),
        (6, (25, 26, 27, 28, 29)),
        (7, (30, 31, 32, 33)),
        (8, (34, 35, 36)),
    )

    coordinates = np.asarray(crystal.cartesian_coordinates)
    chain_numbers = np.asarray(crystal.chain_numbers)
    layer_centers = [
        np.mean([coordinates[chain_numbers == chain_id, 2].mean() for chain_id in chain_ids])
        for _, chain_ids in transverse_chain_layers(crystal)
    ]
    assert layer_centers == sorted(layer_centers)
    for _, chain_ids in transverse_chain_layers(crystal):
        chain_centers_x = [
            coordinates[chain_numbers == chain_id, 0].mean() for chain_id in chain_ids
        ]
        assert chain_centers_x == sorted(chain_centers_x)


def test_rejects_unknown_layer_id() -> None:
    crystal = build_36_chain_fibril(load_allomorph("i_beta"), degree_polymerization=4, vacuum=0.0)

    with pytest.raises(ValueError, match="unknown layer IDs: 99"):
        resolve_chain_restraints(
            crystal,
            chain_ids=None,
            layer_ids=[99],
            force_constants=[1000.0, 1000.0, 1000.0],
        )


def test_groups_free_chains_before_restrained_chains() -> None:
    system = GlycamMappedSystem(
        atom_names=("A", "B", "A", "B", "A", "B"),
        residue_names=("GLC",) * 6,
        residue_numbers=(1,) * 6,
        chain_numbers=(1, 1, 2, 2, 3, 3),
        cartesian_coordinates=np.arange(18).reshape(6, 3),
        basis_vectors=np.eye(3),
    )

    reordered, mapping = group_mapped_chains(system, (2,))

    assert reordered.chain_numbers == (1, 1, 3, 3, 2, 2)
    assert [entry["original_chain_id"] for entry in mapping] == [1, 3, 2]
    assert mapping[-1] == {
        "original_chain_id": 2,
        "gromacs_molecule_index": 3,
        "molecule_type": "CELLULOSE_CHAIN_RESTRAINED",
        "global_atom_start": 5,
        "global_atom_end": 6,
        "restraint_status": "restrained",
    }
