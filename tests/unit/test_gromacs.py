"""Tests for mapped GROMACS output."""

from pathlib import Path

import numpy as np
import pytest

from cellulose_forge.gromacs import _gro_box_line, write_repeated_chain_topology


def test_gro_box_line_converts_orthogonal_angstrom_vectors_to_nanometers() -> None:
    line = _gro_box_line(np.diag([10.0, 20.0, 30.0]))

    assert [float(value) for value in line.split()] == [1.0, 2.0, 3.0]


def test_gro_box_line_writes_gromacs_triclinic_component_order() -> None:
    vectors = np.array([[10.0, 1.0, 2.0], [3.0, 20.0, 4.0], [5.0, 6.0, 30.0]])

    line = _gro_box_line(vectors)

    assert [float(value) for value in line.split()] == pytest.approx(
        [1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    )


def test_repeated_chain_topology_sets_exact_molecule_count(tmp_path: Path) -> None:
    template = tmp_path / "template.top"
    template.write_text(
        "[ moleculetype ]\nCELLULOSE_CHAIN 3\n\n"
        "[ system ]\nCellulose\n\n"
        "[ molecules ]\n; Compound count\nCELLULOSE_CHAIN 1\n",
        encoding="utf-8",
    )
    output = tmp_path / "system.top"

    write_repeated_chain_topology(template, output, chain_count=80)

    assert output.read_text(encoding="utf-8").endswith("CELLULOSE_CHAIN  80\n")
