"""Tests for the experimental CHARMM36 GROMACS adapter."""

from pathlib import Path

import numpy as np
import pytest

import cellulose_forge.charmm36 as charmm36_module
from cellulose_forge.cellulose_ii import cellulose_ii_chain_roles
from cellulose_forge.charmm36 import (
    derive_periodic_charmm36_topology,
    finite_charmm36_interaction_counts,
    map_finite_charmm36,
    map_periodic_charmm36,
    native_finite_charmm36_topology,
    native_periodic_charmm36_topology,
    periodic_charmm36_interaction_counts,
    read_charmm36_chain_template,
    render_finite_charmm36_itp,
    render_periodic_charmm36_itp,
)
from cellulose_forge.clashes import validate_periodic_glycosidic_bonds
from cellulose_forge.crystallography import MillerIndex, load_allomorph
from cellulose_forge.geometry import (
    build_36_chain_fibril,
    build_bulk_crystal,
    build_chain_parallel_slab,
)


def _write_template(path: Path, *, charge: float = 0.0) -> None:
    crystal = build_bulk_crystal(load_allomorph("i_beta"), (1, 1, 1))
    labels = tuple(
        label
        for label, chain, residue in zip(
            crystal.atom_labels,
            crystal.chain_numbers,
            crystal.residue_numbers,
            strict=True,
        )
        if chain == 1 and residue == 1
    )
    lines = ["[ moleculetype ]", "CARB_6 3", "", "[ atoms ]"]
    atom_number = 0
    first_atoms: list[int] = []
    for residue in range(1, 7):
        first_atoms.append(atom_number + 1)
        for label in labels:
            atom_number += 1
            atom_charge = charge if atom_number == 1 else 0.0
            lines.append(
                f"{atom_number} CT {residue} BGLC {label} {atom_number} {atom_charge:.6f} 12.0"
            )
    for section, offsets, function in (
        ("bonds", (0, 1), 1),
        ("pairs", (0, 1), 1),
        ("angles", (0, 1, 2), 5),
        ("dihedrals", (0, 1, 2, 3), 9),
    ):
        lines.extend(("", f"[ {section} ]"))
        lines.extend(
            " ".join(str(first + offset) for offset in offsets) + f" {function}"
            for first in first_atoms
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_reads_and_maps_periodic_charmm36_template(tmp_path: Path) -> None:
    template_path = tmp_path / "chain.itp"
    _write_template(template_path)
    template = read_charmm36_chain_template(template_path)
    topology = derive_periodic_charmm36_topology(template_path, degree_polymerization=4)
    crystal = build_bulk_crystal(load_allomorph("i_beta"), (1, 1, 2))

    mapped = map_periodic_charmm36(crystal, topology)

    assert template.molecule_name == "CARB_6"
    assert template.degree_polymerization == 6
    assert template.total_charge == pytest.approx(0.0)
    assert topology.molecule_name == "CARB_4"
    assert periodic_charmm36_interaction_counts(topology) == {
        "bonds": 4,
        "pairs": 4,
        "angles": 4,
        "dihedrals": 4,
    }
    assert "CARB_4 3" in render_periodic_charmm36_itp(topology)
    assert mapped.chain_count == 2
    assert mapped.atom_count == crystal.atom_count
    assert set(mapped.residue_names) == {"BGLC"}


def test_builds_native_periodic_topology_from_bglc(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    atom_names = (
        "C1",
        "H1",
        "O1",
        "HO1",
        "C5",
        "H5",
        "O5",
        "C2",
        "H2",
        "O2",
        "HO2",
        "C3",
        "H3",
        "O3",
        "HO3",
        "C4",
        "H4",
        "O4",
        "HO4",
        "C6",
        "H61",
        "H62",
        "O6",
        "HO6",
    )
    atom_types = {
        name: (
            "CC3162"
            if name == "C1"
            else "CC3163"
            if name == "C5"
            else "CC3161"
            if name in {"C2", "C3", "C4"}
            else "CC321"
            if name == "C6"
            else "OC3C61"
            if name == "O5"
            else "OC311"
            if name.startswith("O")
            else "HCP1"
            if name.startswith("HO")
            else "HCA2"
            if name in {"H61", "H62"}
            else "HCA1"
        )
        for name in atom_names
    }
    charges = {
        "C1": 0.34,
        "H1": 0.09,
        "O1": -0.65,
        "HO1": 0.42,
        "C5": 0.11,
        "H5": 0.09,
        "O5": -0.40,
        "C2": 0.14,
        "H2": 0.09,
        "O2": -0.65,
        "HO2": 0.42,
        "C3": 0.14,
        "H3": 0.09,
        "O3": -0.65,
        "HO3": 0.42,
        "C4": 0.14,
        "H4": 0.09,
        "O4": -0.65,
        "HO4": 0.42,
        "C6": 0.05,
        "H61": 0.09,
        "H62": 0.09,
        "O6": -0.65,
        "HO6": 0.42,
    }
    bonds = [
        ("C1", "O1"),
        ("C1", "H1"),
        ("O1", "HO1"),
        ("C1", "O5"),
        ("C1", "C2"),
        ("C2", "H2"),
        ("C2", "O2"),
        ("O2", "HO2"),
        ("C2", "C3"),
        ("C3", "H3"),
        ("C3", "O3"),
        ("O3", "HO3"),
        ("C3", "C4"),
        ("C4", "H4"),
        ("C4", "O4"),
        ("O4", "HO4"),
        ("C4", "C5"),
        ("C5", "H5"),
        ("C5", "C6"),
        ("C6", "H61"),
        ("C6", "H62"),
        ("C6", "O6"),
        ("O6", "HO6"),
        ("C5", "O5"),
    ]
    monkeypatch.setattr(
        charmm36_module,
        "_read_bglc_rtp",
        lambda _: ([(name, atom_types[name], charges[name]) for name in atom_names], bonds),
    )
    monkeypatch.setattr(
        charmm36_module,
        "_read_atom_type_masses",
        lambda _: dict.fromkeys(set(atom_types.values()) | {"OC301"}, 12.0),
    )

    topology = native_periodic_charmm36_topology(tmp_path, degree_polymerization=4)

    assert topology.atom_count == 84
    assert topology.total_charge == pytest.approx(0.0)
    assert {atom.atom_name for atom in topology.topology_atoms}.isdisjoint({"O1", "HO1", "HO4"})
    assert periodic_charmm36_interaction_counts(topology) == {
        "bonds": 88,
        "pairs": 252,
        "angles": 164,
        "dihedrals": 264,
    }

    finite = native_finite_charmm36_topology(tmp_path, degree_polymerization=4)
    finite_crystal = build_chain_parallel_slab(
        load_allomorph("i_beta"),
        MillerIndex(1, 1, 0),
        layers=1,
        chains_per_layer=2,
        degree_polymerization=4,
        vacuum=20.0,
        padding_z=12.0,
        chain_boundary="finite",
    )
    finite_mapped = map_finite_charmm36(finite_crystal, finite)
    names_by_residue = {
        residue: {
            atom.atom_name for atom in finite.topology_atoms if atom.residue_number == residue
        }
        for residue in range(1, 5)
    }
    assert finite.atom_count == 87
    assert finite.total_charge == pytest.approx(0.0)
    assert {"O1", "HO1"}.issubset(names_by_residue[1])
    assert "HO4" not in names_by_residue[1]
    assert {"O1", "HO1"}.isdisjoint(names_by_residue[4])
    assert "HO4" in names_by_residue[4]
    assert finite_charmm36_interaction_counts(finite) == {
        "bonds": 90,
        "pairs": 252,
        "angles": 165,
        "dihedrals": 264,
    }
    assert finite_mapped.atom_count == finite_mapped.chain_count * finite.atom_count
    finite_path = tmp_path / "finite.itp"
    finite_path.write_text(
        render_finite_charmm36_itp(finite) + "#ifdef REST_ON\n#endif\n",
        encoding="utf-8",
    )
    reread_finite = read_charmm36_chain_template(
        finite_path, require_uniform_residues=False
    )
    assert len(reread_finite.atoms) == finite.atom_count


def test_rejects_charged_charmm36_template(tmp_path: Path) -> None:
    template_path = tmp_path / "charged.itp"
    _write_template(template_path, charge=0.1)

    with pytest.raises(ValueError, match="not neutral"):
        read_charmm36_chain_template(template_path)


def test_rejects_charmm36_dp_mismatch(tmp_path: Path) -> None:
    template_path = tmp_path / "chain.itp"
    _write_template(template_path)
    topology = derive_periodic_charmm36_topology(template_path, degree_polymerization=4)
    crystal = build_bulk_crystal(load_allomorph("i_beta"), (1, 1, 1))

    with pytest.raises(ValueError, match="topology DP does not match"):
        map_periodic_charmm36(crystal, topology)


@pytest.mark.parametrize(
    ("allomorph", "cells"),
    [
        ("i_alpha", (2, 2, 2)),
        ("i_beta", (2, 2, 2)),
        ("ii", (2, 2, 2)),
        ("iii_i", (2, 2, 2)),
    ],
)
def test_maps_periodic_charmm36_allomorphs(
    tmp_path: Path,
    allomorph: str,
    cells: tuple[int, int, int],
) -> None:
    template_path = tmp_path / "chain.itp"
    _write_template(template_path)
    topology = derive_periodic_charmm36_topology(template_path, degree_polymerization=4)
    crystal = build_bulk_crystal(load_allomorph(allomorph), cells)

    mapped = map_periodic_charmm36(crystal, topology)

    assert mapped.atom_count == mapped.chain_count * topology.atom_count
    assert set(mapped.atom_names) == {atom.atom_name for atom in topology.topology_atoms}
    validation = validate_periodic_glycosidic_bonds(mapped)
    assert validation["bond_count"] == mapped.chain_count * topology.degree_polymerization
    assert validation["maximum_distance_angstrom"] < 2.0


def test_periodic_charmm36_cellulose_ii_hydroxyls_follow_langan_network(
    tmp_path: Path,
) -> None:
    template_path = tmp_path / "chain.itp"
    _write_template(template_path)
    topology = derive_periodic_charmm36_topology(template_path, degree_polymerization=4)
    crystal = build_bulk_crystal(load_allomorph("ii"), (2, 2, 2))
    mapped = map_periodic_charmm36(crystal, topology)
    roles = cellulose_ii_chain_roles(crystal)
    expected = {
        ("origin", "HO2"): ("origin", "O6"),
        ("origin", "HO6"): ("center", "O6"),
        ("center", "HO2"): ("origin", "O2"),
        ("center", "HO6"): ("center", "O2"),
    }
    names = np.asarray(mapped.atom_names)
    chains = np.asarray(mapped.chain_numbers)
    residues = np.asarray(mapped.residue_numbers)
    xyz = np.asarray(mapped.cartesian_coordinates)
    inverse = np.linalg.inv(mapped.basis_vectors)
    for chain, role in roles.items():
        residue = int(residues[chains == chain].min())
        for hydrogen in ("HO2", "HO6"):
            target_role, target_name = expected[role, hydrogen]
            coordinate = xyz[(chains == chain) & (residues == residue) & (names == hydrogen)][0]
            candidates = np.flatnonzero(
                (names == target_name)
                & np.asarray([roles[int(item)] == target_role for item in chains])
            )
            fractional = (xyz[candidates] - coordinate) @ inverse
            fractional -= np.round(fractional)
            distances = np.linalg.norm(fractional @ mapped.basis_vectors, axis=1)
            assert float(distances.min()) < 2.3


def test_maps_cellulose_ii_fibril_in_glycosidic_order(tmp_path: Path) -> None:
    template_path = tmp_path / "chain.itp"
    _write_template(template_path)
    topology = derive_periodic_charmm36_topology(template_path, degree_polymerization=4)
    crystal = build_36_chain_fibril(
        load_allomorph("ii"),
        degree_polymerization=4,
        vacuum=30.0,
        chain_boundary="periodic",
    )

    mapped = map_periodic_charmm36(crystal, topology)
    inverse_box = np.linalg.inv(mapped.basis_vectors)
    linkage_lengths: list[float] = []
    for chain_number in range(1, mapped.chain_count + 1):
        residue_atoms = {
            (residue, atom): coordinate
            for atom, residue, chain, coordinate in zip(
                mapped.atom_names,
                mapped.residue_numbers,
                mapped.chain_numbers,
                mapped.cartesian_coordinates,
                strict=True,
            )
            if chain == chain_number
        }
        first_residue = (chain_number - 1) * topology.degree_polymerization + 1
        for offset in range(topology.degree_polymerization):
            residue = first_residue + offset
            following = first_residue + (offset + 1) % topology.degree_polymerization
            delta = residue_atoms[(following, "C1")] - residue_atoms[(residue, "O4")]
            fractional = delta @ inverse_box
            fractional -= np.round(fractional)
            linkage_lengths.append(float(np.linalg.norm(fractional @ mapped.basis_vectors)))

    assert max(linkage_lengths) < 2.0
