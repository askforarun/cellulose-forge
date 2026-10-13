"""Tests for chain-level position-restraint selection and ordering."""

import numpy as np
import pytest

from cellulose_forge.crystallography import load_allomorph
from cellulose_forge.geometry import build_36_chain_fibril
from cellulose_forge.glycam_mapping import GlycamMappedSystem
from cellulose_forge.restraints import group_mapped_chains, resolve_chain_restraints


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
    crystal = build_36_chain_fibril(
        load_allomorph("i_beta"), degree_polymerization=4, vacuum=0.0
    )

    with pytest.raises(ValueError, match="unknown chain IDs: 37"):
        resolve_chain_restraints(
            crystal,
            chain_ids=[37],
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
