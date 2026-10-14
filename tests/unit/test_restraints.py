"""Tests for chain-level position-restraint selection and ordering."""

import numpy as np
import pytest

from cellulose_forge.crystallography import MillerIndex, load_allomorph
from cellulose_forge.geometry import build_36_chain_fibril, build_chain_parallel_slab
from cellulose_forge.glycam_mapping import GlycamMappedSystem
from cellulose_forge.restraints import (
    group_mapped_chains,
    resolve_chain_restraints,
    transverse_chain_layers,
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


def test_rejects_unknown_chain_id() -> None:
    crystal = build_36_chain_fibril(load_allomorph("i_beta"), degree_polymerization=4, vacuum=0.0)

    with pytest.raises(ValueError, match="unknown chain IDs: 37"):
        resolve_chain_restraints(
            crystal,
            chain_ids=[37],
            force_constants=[1000.0, 1000.0, 1000.0],
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


def test_rejects_combined_chain_and_layer_selectors() -> None:
    crystal = build_36_chain_fibril(load_allomorph("i_beta"), degree_polymerization=4, vacuum=0.0)

    with pytest.raises(ValueError, match="mutually exclusive"):
        resolve_chain_restraints(
            crystal,
            chain_ids=[1],
            layer_ids=[1],
            force_constants=[1000.0, 1000.0, 1000.0],
        )


def test_fibril_layers_progress_from_low_to_high_y_after_rotation() -> None:
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
        np.mean([coordinates[chain_numbers == chain_id, 1].mean() for chain_id in chain_ids])
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
