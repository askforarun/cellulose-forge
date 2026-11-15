"""Tests for native NAMD/X-PLOR PSF output."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from cellulose_forge.charmm36 import (
    Charmm36TemplateAtom,
    FiniteCharmm36Topology,
    _IndexedInteraction,
)
from cellulose_forge.glycam_mapping import GlycamMappedSystem
from cellulose_forge.namd import (
    render_namd_configuration,
    render_pdb,
    render_psf,
    validate_namd_one_four_pairs,
    write_namd_system,
)


def _topology(*, include_pair: bool = True) -> FiniteCharmm36Topology:
    atoms = tuple(
        Charmm36TemplateAtom(1, "BGLC", name, charge, atom_type, mass)
        for name, charge, atom_type, mass in (
            ("C1", 0.2, "CC3162", 12.011),
            ("O5", -0.3, "OC3C61", 15.999),
            ("C5", 0.2, "CC3163", 12.011),
            ("O4", -0.1, "OC301", 15.999),
        )
    )
    interactions = [
        _IndexedInteraction("bonds", (0, 1), ("1",)),
        _IndexedInteraction("bonds", (1, 2), ("1",)),
        _IndexedInteraction("bonds", (2, 3), ("1",)),
        _IndexedInteraction("angles", (0, 1, 2), ("5",)),
        _IndexedInteraction("angles", (1, 2, 3), ("5",)),
        _IndexedInteraction("dihedrals", (0, 1, 2, 3), ("9",)),
        _IndexedInteraction("dihedrals", (1, 0, 2, 3), ("2",)),
    ]
    if include_pair:
        interactions.append(_IndexedInteraction("pairs", (0, 3), ("1",)))
    return FiniteCharmm36Topology(1, "TEST", atoms, tuple(interactions))


def _system() -> GlycamMappedSystem:
    coordinates = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.4, 0.0, 0.0],
            [2.8, 0.0, 0.0],
            [4.2, 0.0, 0.0],
            [0.0, 5.0, 0.0],
            [1.4, 5.0, 0.0],
            [2.8, 5.0, 0.0],
            [4.2, 5.0, 0.0],
        ]
    )
    return GlycamMappedSystem(
        atom_names=("C1", "O5", "C5", "O4") * 2,
        residue_names=("BGLC",) * 8,
        residue_numbers=(1,) * 8,
        chain_numbers=(1,) * 4 + (2,) * 4,
        cartesian_coordinates=coordinates,
        basis_vectors=np.asarray([[20.0, 0.0, 0.0], [2.0, 21.0, 0.0], [1.0, 3.0, 22.0]]),
    )


def test_psf_expands_chains_and_writes_zero_exclusion_array() -> None:
    text = render_psf(_system(), _topology())

    assert "         8 !NATOM" in text
    assert "         6 !NBOND: bonds" in text
    assert "         4 !NTHETA: angles" in text
    assert "         2 !NPHI: dihedrals" in text
    assert "         2 !NIMPHI: impropers" in text
    assert "C001" in text
    assert "C002" in text
    exclusion_line = text.split("!NNB\n", maxsplit=1)[1].splitlines()[0]
    assert [int(value) for value in exclusion_line.split()] == [0] * 8


def test_pair_guard_rejects_missing_graph_distance_three_pair() -> None:
    with pytest.raises(ValueError, match="missing 1, extra 0"):
        validate_namd_one_four_pairs(_topology(include_pair=False))


def test_pdb_order_and_triclinic_configuration() -> None:
    system = _system()
    pdb = render_pdb(system, _topology())
    configuration = render_namd_configuration(system, periodic_chains=True)

    atom_lines = [line for line in pdb.splitlines() if line.startswith("ATOM")]
    assert len(atom_lines) == 8
    assert atom_lines[0][12:16].strip() == "C1"
    assert atom_lines[0][17:21] == "BGLC"
    assert atom_lines[0][21] == "A"
    assert float(atom_lines[0][30:38]) == 0.0
    assert atom_lines[3][12:16].strip() == "O4"
    assert "C001" in atom_lines[0]
    assert "C002" in atom_lines[4]
    assert "cellBasisVector1 20.00000000 0.00000000 0.00000000" in configuration
    assert "cellBasisVector2 2.00000000 21.00000000 0.00000000" in configuration
    assert "cellBasisVector3 1.00000000 3.00000000 22.00000000" in configuration
    assert "cellOrigin 11.50000000 12.00000000 11.00000000" in configuration
    assert "coordinates structure_namd.pdb" in configuration
    assert "exclude scaled1-4" in configuration
    assert "oneFourScaling 1.0" in configuration
    assert "vdwForceSwitching on" in configuration
    assert "wrapAll off" in configuration


def test_pdb_uses_fixed_width_alias_for_six_character_chitin_residue() -> None:
    system = replace(_system(), residue_names=("BGLCNA",) * 8)

    atom_line = render_pdb(system, _topology()).splitlines()[0]

    assert atom_line[17:21].strip() == "NAG"
    assert float(atom_line[30:38]) == 0.0


def test_writer_copies_parameter_file_and_writes_complete_bundle(tmp_path: Path) -> None:
    parameters = tmp_path / "source.prm"
    parameters.write_text("* test parameters\n", encoding="utf-8")

    files = write_namd_system(
        _system(), _topology(), parameters, tmp_path / "namd", periodic_chains=False
    )

    assert files.psf.is_file()
    assert files.pdb.is_file()
    assert files.configuration.is_file()
    assert files.parameters.read_text(encoding="utf-8") == "* test parameters\n"
