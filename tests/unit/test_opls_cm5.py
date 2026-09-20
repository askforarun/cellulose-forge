"""Tests for the experimental OPLS-CM5 reference adapter."""

from pathlib import Path

import pytest

from cellulose_forge.opls_cm5 import (
    PeriodicOplsTopology,
    _InteractionTemplate,
    _read_first_reference_chain,
    _read_topology_atoms,
    _TopologyAtom,
    derive_periodic_opls_cm5_topology,
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


def test_periodic_opls_requires_even_degree_polymerization() -> None:
    with pytest.raises(ValueError, match="even degree of polymerization"):
        derive_periodic_opls_cm5_topology("missing.itp", "missing.gro", degree_polymerization=3)


def test_periodic_opls_render_wraps_interactions_and_neutralizes_charge() -> None:
    topology = PeriodicOplsTopology(
        degree_polymerization=4,
        coordinate_atom_names=("A", "B"),
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
