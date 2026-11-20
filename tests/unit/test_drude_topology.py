from __future__ import annotations

import os
from pathlib import Path

import pytest
from drude_fixtures import MAIN, PATCH, RESIDUE, carbohydrate

from cellulose_forge.drude_topology import (
    DRUDE_CARBOHYDRATE_FILE,
    DRUDE_MAIN_FILE,
    DRUDE_MASS,
    NAMD_COULOMB_CONSTANT,
    DrudeChainTopology,
    anisotropy_constants,
    build_drude_cellulose_chain,
    drude_charge,
    parse_drude_cellulose_definitions,
    read_drude_cellulose_definitions,
    thole_screened_pairs,
)


def _chain(dp: int, *, periodic: bool, **kwargs: str) -> DrudeChainTopology:
    definitions = parse_drude_cellulose_definitions(MAIN, carbohydrate(**kwargs))
    return build_drude_cellulose_chain(definitions, dp, periodic=periodic)


def _names(topology: DrudeChainTopology, indices: tuple[int, ...]) -> tuple[tuple[int, str], ...]:
    return tuple((topology.sites[index].residue, topology.sites[index].name) for index in indices)


def test_parser_reads_every_supported_directive() -> None:
    definitions = parse_drude_cellulose_definitions(MAIN, carbohydrate())
    assert [atom.name for atom in definitions.residue.atoms][:4] == ["C1", "H1", "O1", "HO1"]
    assert len(definitions.residue.bonds) == 11
    assert len(definitions.residue.lone_pairs) == 4
    assert len(definitions.residue.anisotropies) == 2
    assert definitions.patch.deleted_atoms[5] == "2DO1"
    assert definitions.patch.lone_pairs[0].mode == "bisector"
    assert definitions.parameters.drude_force_constant == 500.0


@pytest.mark.parametrize(
    "line",
    [
        "IMPR C1 H1 C4 O1",
        "DOUBLE C1 C4",
        "TRIPLE C1 C4",
        "CMAP C1 C4 O4 HO4 C1 C4 O4 HO4",
        "ANGL O1 C1 C4",
        "DIHE O1 C1 C4 O4",
        "dele bond C1 C4",
        "dele ANISOTROPY O1 C1 LP1A LP1B",
        "BOND C1 -",
        "BOND C1",
        "ATOM X1 TC 0.0 TYPE 3",
        "ATOM X1 TC 0.0 ALPHA -1.0",
        "LONEPAIR colinear LPZ O1 C1 HO1 distance 0.35 angle 110.0 dihe 91.0",
        "PATCH FIRST NTER LAST NONE",
        "ATOM C1 TC 0.100",
    ],
)
def test_parser_refuses_unsupported_or_malformed_directives(line: str) -> None:
    with pytest.raises(ValueError, match=r"unsupported|malformed|continuation|more than once|must"):
        parse_drude_cellulose_definitions(MAIN, carbohydrate(residue=RESIDUE + line + "\n"))


def test_parser_refuses_duplicate_definitions_and_missing_masses() -> None:
    with pytest.raises(ValueError, match="more than once"):
        parse_drude_cellulose_definitions(MAIN, carbohydrate() + RESIDUE)
    with pytest.raises(ValueError, match="lacks masses"):
        parse_drude_cellulose_definitions(MAIN.replace("MASS -1 TE", "! MASS"), carbohydrate())
    with pytest.raises(ValueError, match="exactly one 'X DRUD'"):
        parse_drude_cellulose_definitions(MAIN.replace("X DRUD 500.0 0.0\n", ""), carbohydrate())


def test_patch_may_delete_only_the_drude_of_a_deleted_polarizable_atom() -> None:
    with pytest.raises(ValueError, match="deletes absent atom 2:DUMMY"):
        _chain(2, periodic=False, patch=PATCH.replace("dele atom 2DO1", "dele atom 2DUMMY"))
    # DC4 is the Drude of a polarizable atom that the patch keeps.
    with pytest.raises(ValueError, match="deletes absent atom 2:DC4"):
        _chain(2, periodic=False, patch=PATCH.replace("dele atom 2DO1", "dele atom 2DC4"))
    # DH1 names no Drude: H1 is not polarizable.
    with pytest.raises(ValueError, match="deletes absent atom 2:DH1"):
        _chain(2, periodic=False, patch=PATCH.replace("dele atom 2DO1", "dele atom 2DH1"))


def test_finite_chain_patches_donor_and_acceptor_sides_of_internal_residues() -> None:
    topology = _chain(3, periodic=False)
    by_name = {(site.residue, site.name): site for site in topology.sites}
    # Residue 2 is both the O4 donor of link 2-3 and the C1 acceptor of link 1-2.
    assert by_name[(2, "O4")].atom_type == "TE"
    assert by_name[(2, "C4")].alpha == -1.050
    assert by_name[(2, "C1")].alpha == -1.020
    assert (2, "O1") not in by_name and (2, "HO4") not in by_name
    # Chain ends keep the unpatched caps.
    assert by_name[(1, "O1")].atom_type == "TO" and (1, "HO1") in by_name
    assert by_name[(3, "O4")].atom_type == "TO" and (3, "HO4") in by_name
    assert (1, "DO1") in by_name and (2, "DO1") not in by_name
    assert abs(topology.total_charge) < 1.0e-10


def test_sites_keep_stream_order_with_each_drude_after_its_parent() -> None:
    topology = _chain(2, periodic=False)
    assert [site.name for site in topology.sites if site.residue == 1] == [
        "C1", "DC1", "H1", "O1", "DO1", "HO1", "LP1A", "LP1B",
        "C4", "DC4", "H4", "O4", "DO4", "LPDA4", "LPDB4", "LPXD4",
    ]  # fmt: skip
    for drude, parent in topology.drude_parents:
        assert drude == parent + 1
        assert topology.sites[drude].kind == "drude"
        assert topology.sites[drude].mass == DRUDE_MASS
        assert topology.sites[drude].charge < 0.0


def test_parent_drude_bonds_follow_every_other_bond() -> None:
    # NAMD 3.0.3 rejects a Drude PSF whose parent--Drude bonds are mixed in.
    for periodic in (False, True):
        topology = _chain(3, periodic=periodic)
        kinds = [
            "drude" in {topology.sites[first].kind, topology.sites[second].kind}
            for first, second in topology.bonds
        ]
        count = len(topology.drude_parents)
        assert kinds == [False] * (len(kinds) - count) + [True] * count
        assert all(
            topology.sites[second].kind == "drude" for _first, second in topology.bonds[-count:]
        )
        assert len(set(topology.bonds)) == len(topology.bonds)


def test_drude_particles_conserve_charge_and_mass() -> None:
    definitions = parse_drude_cellulose_definitions(MAIN, carbohydrate())
    topology = build_drude_cellulose_chain(definitions, 2, periodic=False)
    parent = topology.sites[topology.index(1, "C1")]
    drude = topology.sites[topology.index(1, "DC1")]
    assert parent.charge + drude.charge == pytest.approx(0.100)
    assert parent.mass + drude.mass == pytest.approx(12.011)
    assert parent.alpha == -1.0 and parent.thole == 1.3
    assert drude.charge == pytest.approx(drude_charge(-1.0, 500.0, topology.coulomb_constant))


def test_angles_and_dihedrals_cover_the_linkage_and_skip_special_sites() -> None:
    topology = _chain(3, periodic=False)
    kinds = [site.kind for site in topology.sites]
    for term in (*topology.angles, *topology.dihedrals):
        assert all(kinds[index] == "atom" for index in term)
    angles = {_names(topology, angle) for angle in topology.angles}
    angles |= {angle[::-1] for angle in angles}
    assert ((1, "C4"), (1, "O4"), (2, "C1")) in angles
    assert ((1, "O4"), (2, "C1"), (2, "H1")) in angles
    dihedrals = {_names(topology, dihedral) for dihedral in topology.dihedrals}
    dihedrals |= {dihedral[::-1] for dihedral in dihedrals}
    assert ((1, "C1"), (1, "C4"), (1, "O4"), (2, "C1")) in dihedrals
    assert ((1, "C4"), (1, "O4"), (2, "C1"), (2, "C4")) in dihedrals
    assert len(topology.angles) == len(set(topology.angles))
    assert len(topology.dihedrals) == len(set(topology.dihedrals))


def test_cyclic_chain_closes_last_residue_onto_first() -> None:
    topology = _chain(3, periodic=True)
    bonds = {_names(topology, bond) for bond in topology.bonds}
    assert ((1, "C1"), (3, "O4")) in bonds
    names = {(site.residue, site.name) for site in topology.sites}
    assert not {name for _residue, name in names} & {"O1", "HO1", "HO4"}
    per_residue = [sum(site.residue == residue for site in topology.sites) for residue in (1, 2, 3)]
    assert len(set(per_residue)) == 1
    closing = next(
        item for item in topology.lone_pairs if _names(topology, (item.site,)) == ((3, "LPDA4"),)
    )
    assert _names(topology, closing.hosts) == ((3, "O4"), (3, "C4"), (1, "C1"))
    assert closing.mode == "bisector"
    assert abs(topology.total_charge) < 1.0e-10


def test_anisotropy_follows_the_patched_linkage_oxygen() -> None:
    topology = _chain(2, periodic=False)
    records = {_names(topology, item.atoms): item for item in topology.anisotropies}
    assert ((1, "O4"), (1, "LPXD4"), (1, "LPDA4"), (1, "LPDB4")) in records
    assert ((1, "O4"), (1, "C4"), (1, "LP4A"), (1, "LP4B")) not in records
    assert ((2, "O1"), (2, "C1"), (2, "LP1A"), (2, "LP1B")) not in records
    assert ((2, "O4"), (2, "C4"), (2, "LP4A"), (2, "LP4B")) in records
    assert len(topology.anisotropies) == 3


def test_missing_bonded_parameter_is_refused() -> None:
    definitions = parse_drude_cellulose_definitions(
        MAIN.replace("TC TE TC 50.0 115.0\n", ""), carbohydrate()
    )
    with pytest.raises(ValueError, match="no parameters for: angle TC-TE-TC"):
        build_drude_cellulose_chain(definitions, 2, periodic=False)


def test_charged_chain_and_short_chain_are_refused() -> None:
    with pytest.raises(ValueError, match="not neutral"):
        _chain(2, periodic=False, patch=PATCH.replace("ATOM 2C1 TC 0.150", "ATOM 2C1 TC 0.250"))
    with pytest.raises(ValueError, match="DP of at least 2"):
        _chain(1, periodic=False)


def test_drude_charge_and_anisotropy_constants_match_psfgen_values() -> None:
    # Values printed by psfgen 2.0 for BGLC C1 and the two BGLC anisotropy kinds.
    assert drude_charge(-1.041, 500.0, NAMD_COULOMB_CONSTANT) == pytest.approx(-1.770577, abs=5e-7)
    assert anisotropy_constants(0.8108, 1.2162, 500.0) == pytest.approx(
        (102.8003, -102.7580, 13.8746), abs=5e-5
    )
    assert anisotropy_constants(0.8889, 1.2222, 500.0) == pytest.approx(
        (0.0, -153.3946, 62.4930), abs=5e-5
    )
    with pytest.raises(ValueError, match="negative CHARMM ALPHA"):
        drude_charge(1.0, 500.0, NAMD_COULOMB_CONSTANT)
    with pytest.raises(ValueError, match="must be positive"):
        anisotropy_constants(2.0, 1.5, 500.0)


def test_thole_pairs_are_one_two_and_one_three_polarizable_pairs() -> None:
    topology = _chain(2, periodic=False)
    pairs = {
        (_names(topology, (first,))[0], _names(topology, (second,))[0]): factor
        for first, second, factor in thole_screened_pairs(topology)
    }
    assert pairs[((1, "C1"), (1, "O1"))] == pytest.approx(1.3 + 0.5)
    assert pairs[((1, "O4"), (2, "C1"))] == pytest.approx(1.0 + 0.46)
    assert pairs[((1, "C4"), (2, "C1"))] == pytest.approx(0.2 + 0.46)  # 1-3 across the link
    assert ((1, "O1"), (1, "O4")) not in pairs  # 1-4
    assert all(first < second for first, second, _factor in thole_screened_pairs(topology))


def test_reader_refuses_streams_that_were_not_audited(tmp_path: Path) -> None:
    main = tmp_path / DRUDE_MAIN_FILE
    stream = tmp_path / DRUDE_CARBOHYDRATE_FILE
    main.write_text(MAIN, encoding="utf-8")
    stream.write_text(carbohydrate(), encoding="utf-8")
    with pytest.raises(ValueError, match="is not the audited CHARMM Drude 2023a"):
        read_drude_cellulose_definitions(main, stream)
    with pytest.raises(ValueError, match="stream not found"):
        read_drude_cellulose_definitions(tmp_path / "absent.str", stream)


_TOPPAR = os.environ.get("CELLULOSE_FORGE_DRUDE_TOPPAR")


@pytest.mark.skipif(
    _TOPPAR is None, reason="set CELLULOSE_FORGE_DRUDE_TOPPAR to the official streams"
)
@pytest.mark.parametrize(
    ("dp", "periodic", "expected"),
    [
        # explicit sites, explicit bonds, lone pairs, anisotropy, Drudes, angles, dihedrals
        (2, False, (70, 71, 25, 11, 23, 83, 132)),
        (8, False, (268, 275, 97, 41, 89, 329, 528)),
        (8, True, (264, 272, 96, 40, 88, 328, 528)),
    ],
)
def test_official_stream_audit_counts(dp: int, periodic: bool, expected: tuple[int, ...]) -> None:
    assert _TOPPAR is not None
    root = Path(_TOPPAR)
    definitions = read_drude_cellulose_definitions(
        root / DRUDE_MAIN_FILE, root / DRUDE_CARBOHYDRATE_FILE
    )
    topology = build_drude_cellulose_chain(definitions, dp, periodic=periodic)
    assert (
        topology.explicit_site_count,
        topology.explicit_bond_count,
        len(topology.lone_pairs),
        len(topology.anisotropies),
        len(topology.drude_parents),
        len(topology.angles),
        len(topology.dihedrals),
    ) == expected
    assert abs(topology.total_charge) < 1.0e-8
