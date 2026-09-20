"""Tests for the experimental OPLS-CM5 reference adapter."""

from pathlib import Path

import pytest

from cellulose_forge.opls_cm5 import _read_first_reference_chain, _read_topology_atoms


def test_opls_reference_topology_must_be_neutral(tmp_path: Path) -> None:
    topology = tmp_path / "CNC.itp"
    rows = [
        f"{index} opls_140 1 CNC HC {index} {0.1 if index == 1 else 0.0} 1.008"
        for index in range(1, 424)
    ]
    topology.write_text("[ atoms ]\n" + "\n".join(rows) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="not neutral"):
        _read_topology_atoms(topology)


def test_opls_reference_gro_requires_published_dp20_layout(tmp_path: Path) -> None:
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
