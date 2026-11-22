"""Regression tests for the optional experimental OpenFF Sage backend."""

from __future__ import annotations

import importlib
import importlib.util
import sys
import time
from collections import Counter
from dataclasses import replace
from importlib.util import find_spec
from pathlib import Path

import numpy as np
import pytest
from hydroxyl_contacts import modelled_hydroxyl_contacts

from cellulose_forge import cli
from cellulose_forge.cli import main
from cellulose_forge.crystallography import MillerIndex, load_allomorph
from cellulose_forge.geometry import build_bulk_crystal, build_chain_parallel_slab
from cellulose_forge.glycam_mapping import map_crystal_to_periodic_glycam
from cellulose_forge.openff_sage import (
    _REFERENCE_PARITY,
    _centre_parity,
    build_sage_molecule,
    build_sage_system,
    dispersion_correction_for_shape,
    reduced_periodic_box,
    require_openff,
    validate_declared_stereochemistry,
    validate_sage_stereochemistry,
)


def test_dispersion_policy_is_geometry_specific() -> None:
    assert dispersion_correction_for_shape("bulk") is True
    assert dispersion_correction_for_shape("crystal") is True
    assert dispersion_correction_for_shape("slab") is False
    assert dispersion_correction_for_shape("fibril") is False
    with pytest.raises(ValueError, match="undefined"):
        dispersion_correction_for_shape("cylinder")


def test_missing_optional_dependency_has_installation_message(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real_import = importlib.import_module

    def missing(name: str) -> object:
        if name == "openff.toolkit":
            raise ImportError("not installed")
        return real_import(name)

    monkeypatch.setattr(importlib, "import_module", missing)
    with pytest.raises(RuntimeError, match=r"environment-openff\.yml"):
        require_openff()


def test_library_api_refuses_mannan_as_cellulose_sage_graph() -> None:
    crystal = build_bulk_crystal(load_allomorph("i", polymer="mannan"), (2, 2, 2))
    with pytest.raises(ValueError, match=r"requires crystal residue BGC, not MAN"):
        build_sage_system(crystal, "periodic", "cellulose")


def test_triclinic_box_reduction_preserves_lattice_volume() -> None:
    box = np.asarray(((10.0, 0.0, 0.0), (-7.0, 9.0, 0.0), (6.0, 5.0, 8.0)))
    reduced = reduced_periodic_box(box)
    assert np.linalg.det(reduced) == pytest.approx(np.linalg.det(box))
    assert abs(reduced[1, 0]) <= reduced[0, 0] / 2
    assert abs(reduced[2, 0]) <= reduced[0, 0] / 2
    assert abs(reduced[2, 1]) <= reduced[1, 1] / 2


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
def test_periodic_mirror_reflects_coordinates_and_triclinic_box() -> None:
    root = Path(__file__).resolve().parents[2]
    validation = root / "validation"
    sys.path.insert(0, str(validation))
    spec = importlib.util.spec_from_file_location(
        "_test_validate_openff_periodic_slab",
        validation / "validate_openff_periodic_slab.py",
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    box = np.asarray(((4.0, 0.0, 0.0), (-1.0, 3.0, 0.0), (0.0, 0.0, 5.0)))
    positions = np.asarray(((0.2, 0.3, 0.4), (3.9, 2.8, 0.5)))
    displacement = positions[1] - positions[0]
    mirrored_positions, mirrored_box = module.mirror_x_periodic(positions, box)

    expected = module.minimum_image(displacement, box).copy()
    expected[0] *= -1.0
    observed = module.minimum_image(mirrored_positions[1] - mirrored_positions[0], mirrored_box)
    assert observed == pytest.approx(expected)


def _angle(left: np.ndarray, centre: np.ndarray, right: np.ndarray) -> float:
    first = left - centre
    second = right - centre
    cosine = np.dot(first, second) / (np.linalg.norm(first) * np.linalg.norm(second))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
def test_direct_graph_counts_and_periodic_scaling() -> None:
    toolkit = importlib.import_module("openff.toolkit")
    force_field_type, topology_type = toolkit.ForceField, toolkit.Topology
    force_field = force_field_type("openff-2.3.0.offxml")
    normalized: dict[int, dict[str, Counter[str]]] = {}
    for degree in (4, 6):
        molecule = build_sage_molecule(degree, "periodic")
        assert molecule.n_atoms == 21 * degree
        assert molecule.n_bonds == 22 * degree
        labels = force_field.label_molecules(topology_type.from_molecules([molecule]))[0]
        normalized[degree] = {}
        for handler in ("Bonds", "Angles", "ProperTorsions"):
            counts = Counter(parameter.id for parameter in labels[handler].values())
            assert all(count % degree == 0 for count in counts.values())
            normalized[degree][handler] = Counter(
                {parameter: count // degree for parameter, count in counts.items()}
            )
    assert normalized[4] == normalized[6]


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
@pytest.mark.parametrize("boundary", ["finite", "periodic"])
def test_chitin_direct_graph_has_glcnac_chemistry(boundary: str) -> None:
    molecule = build_sage_molecule(6, boundary, "chitin")
    assert molecule.n_atoms == (165 if boundary == "finite" else 162)
    assert molecule.n_bonds == (170 if boundary == "finite" else 168)
    carbonyls = [
        bond for bond in molecule.bonds if {bond.atom1.name, bond.atom2.name} == {"C7", "O7"}
    ]
    assert len(carbonyls) == 6
    assert all(int(bond.bond_order) == 2 for bond in carbonyls)
    assert {atom.metadata["residue_name"] for atom in molecule.atoms} == {"NAG"}
    terminal = {
        (int(atom.metadata["residue_number"]), atom.name): atom.stereochemistry
        for atom in molecule.atoms
        if atom.name in {"C1", "C3"}
    }
    if boundary == "finite":
        assert terminal[1, "C1"] == "R"
        assert terminal[6, "C3"] == "R"


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
@pytest.mark.parametrize("allomorph", ["alpha", "beta"])
def test_chitin_periodic_nagl_charges_and_stereochemistry(allomorph: str) -> None:
    crystal = build_bulk_crystal(load_allomorph(allomorph, polymer="chitin"), (4, 2, 2))
    build = build_sage_system(crystal, "periodic", "chitin")
    assert build.charge_provenance["net_charge_e"] == pytest.approx(0.0, abs=1.0e-6)
    assert build.charge_provenance["periodic_equivalent_atom_charges"] is True
    assert build.charge_provenance["polymer"] == "chitin"
    assert build.stereochemistry["passed"] is True
    assert set(build.mapped.residue_names) == {"NAG"}


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
def test_chitin_stereochemistry_rejects_mirror_and_local_inversion() -> None:
    crystal = build_bulk_crystal(load_allomorph("beta", polymer="chitin"), (4, 2, 2))
    build = build_sage_system(crystal, "periodic", "chitin")
    mirrored = build.mapped.cartesian_coordinates.copy()
    mirrored[:, 0] *= -1
    with pytest.raises(ValueError, match="OpenFF chitin stereochemistry check failed"):
        validate_sage_stereochemistry(
            build.molecule,
            replace(build.mapped, cartesian_coordinates=mirrored),
            "chitin",
        )

    inverted = build.mapped.cartesian_coordinates.copy()
    centre = next(
        atom
        for atom in build.molecule.atoms
        if atom.name == "C2" and int(atom.metadata["residue_number"]) == 1
    )
    first, second = sorted(centre.bonded_atoms, key=lambda atom: (atom.symbol, atom.name))[:2]
    inverted[[first.molecule_atom_index, second.molecule_atom_index]] = inverted[
        [second.molecule_atom_index, first.molecule_atom_index]
    ]
    with pytest.raises(ValueError, match="OpenFF chitin stereochemistry check failed"):
        validate_sage_stereochemistry(
            build.molecule,
            replace(build.mapped, cartesian_coordinates=inverted),
            "chitin",
        )


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
@pytest.mark.parametrize("polymer", ["cellulose", "chitin"])
@pytest.mark.parametrize("boundary", ["finite", "periodic"])
def test_declared_stereochemistry_matches_mapped_3d_coordinates(
    polymer: str, boundary: str
) -> None:
    """Every graph label must agree with stereochemistry perceived from its coordinates."""
    chem = importlib.import_module("rdkit.Chem")
    point_type = importlib.import_module("rdkit.Geometry").Point3D
    if polymer == "cellulose":
        definition = load_allomorph("i_beta")
        face = MillerIndex(1, 1, 0)
    else:
        definition = load_allomorph("beta", polymer="chitin")
        face = MillerIndex(1, 0, 0)
    crystal = build_chain_parallel_slab(
        definition,
        face,
        layers=2,
        chains_per_layer=2,
        degree_polymerization=4,
        vacuum=30.0,
        padding_z=12.0 if boundary == "finite" else 0.0,
        chain_boundary=boundary,
    )
    build = build_sage_system(crystal, boundary, polymer)
    rdkit_molecule = build.molecule.to_rdkit()
    conformer = chem.Conformer(rdkit_molecule.GetNumAtoms())
    for index, coordinate in enumerate(
        build.mapped.cartesian_coordinates[: rdkit_molecule.GetNumAtoms()]
    ):
        conformer.SetAtomPosition(index, point_type(*map(float, coordinate)))
    rdkit_molecule.RemoveAllConformers()
    rdkit_molecule.AddConformer(conformer)
    for atom in rdkit_molecule.GetAtoms():
        atom.SetChiralTag(chem.ChiralType.CHI_UNSPECIFIED)
        if atom.HasProp("_CIPCode"):
            atom.ClearProp("_CIPCode")
    chem.AssignStereochemistryFrom3D(rdkit_molecule, confId=0, replaceExistingTags=True)

    mismatches = []
    for index, atom in enumerate(build.molecule.atoms):
        if atom.stereochemistry is None:
            continue
        perceived = rdkit_molecule.GetAtomWithIdx(index)
        actual = perceived.GetProp("_CIPCode") if perceived.HasProp("_CIPCode") else None
        if actual != atom.stereochemistry:
            mismatches.append(
                (int(atom.metadata["residue_number"]), atom.name, atom.stereochemistry, actual)
            )
    assert mismatches == []
    declared = build.stereochemistry["declared_labels"]
    assert declared["passed"] is True
    assert declared["distinct_chain_types_checked"] >= 1


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
def test_build_time_stereochemistry_guard_rejects_wrong_graph_label() -> None:
    crystal = build_chain_parallel_slab(
        load_allomorph("beta", polymer="chitin"),
        MillerIndex(1, 0, 0),
        layers=2,
        chains_per_layer=2,
        degree_polymerization=4,
        vacuum=30.0,
        padding_z=12.0,
        chain_boundary="finite",
    )
    build = build_sage_system(crystal, "finite", "chitin")
    terminal_c3 = next(
        atom
        for atom in build.molecule.atoms
        if atom.name == "C3" and int(atom.metadata["residue_number"]) == 4
    )
    terminal_c3.stereochemistry = "S"
    with pytest.raises(ValueError, match="declared stereochemistry disagrees"):
        validate_declared_stereochemistry(build.molecule, build.mapped, "finite", "chitin")


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
@pytest.mark.parametrize("polymer", ["cellulose", "chitin"])
def test_declared_stereo_guard_handles_wrapped_periodic_chains(polymer: str) -> None:
    if polymer == "cellulose":
        crystal = build_bulk_crystal(load_allomorph("i_beta"), (3, 3, 3))
    else:
        crystal = build_bulk_crystal(load_allomorph("beta", polymer="chitin"), (4, 2, 2))
    build = build_sage_system(crystal, "periodic", polymer)
    inverse = np.linalg.inv(build.mapped.basis_vectors)
    fractional = build.mapped.cartesian_coordinates @ inverse
    wrapped = (fractional - np.floor(fractional)) @ build.mapped.basis_vectors
    result = validate_declared_stereochemistry(
        build.molecule,
        replace(build.mapped, cartesian_coordinates=wrapped),
        "periodic",
        polymer,
    )
    assert result["passed"] is True


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
@pytest.mark.parametrize(
    ("degree", "pubchem_cid", "isomeric_smiles"),
    (
        (
            2,
            10712,
            "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O[C@@H]2[C@H](O[C@H]"
            "([C@@H]([C@H]2O)O)O)CO)O)O)O)O",
        ),
        (
            6,
            440948,
            "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O[C@@H]2[C@H](O[C@H]"
            "([C@@H]([C@H]2O)O)O[C@@H]3[C@H](O[C@H]([C@@H]([C@H]3O)O)O"
            "[C@@H]4[C@H](O[C@H]([C@@H]([C@H]4O)O)O[C@@H]5[C@H](O[C@H]"
            "([C@@H]([C@H]5O)O)O[C@@H]6[C@H](O[C@H]([C@@H]([C@H]6O)O)O)"
            "CO)CO)CO)CO)CO)O)O)O)O",
        ),
    ),
)
def test_finite_graph_matches_pubchem_stereochemistry(
    degree: int, pubchem_cid: int, isomeric_smiles: str
) -> None:
    molecule_type = importlib.import_module("openff.toolkit").Molecule
    generated = build_sage_molecule(degree, "finite")
    reference = molecule_type.from_smiles(isomeric_smiles, allow_undefined_stereo=False)
    isomorphic, _ = molecule_type.are_isomorphic(
        generated,
        reference,
        atom_stereochemistry_matching=True,
        bond_stereochemistry_matching=True,
    )
    assert isomorphic, f"finite DP{degree} graph differs from PubChem CID {pubchem_cid}"


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
def test_cellulose_ii_modelled_hydrogens_are_tetrahedral() -> None:
    crystal = build_bulk_crystal(load_allomorph("ii"), (1, 1, 2))
    build = build_sage_system(crystal, "periodic")
    for chain in set(build.mapped.chain_numbers):
        for residue in set(build.mapped.residue_numbers):
            atoms = {
                name: coordinate
                for name, coordinate, atom_chain, atom_residue in zip(
                    build.mapped.atom_names,
                    build.mapped.cartesian_coordinates,
                    build.mapped.chain_numbers,
                    build.mapped.residue_numbers,
                    strict=True,
                )
                if atom_chain == chain and atom_residue == residue
            }
            for suffix in ("2", "3", "6"):
                assert np.linalg.norm(atoms[f"HO{suffix}"] - atoms[f"O{suffix}"]) == pytest.approx(
                    0.96, abs=1.0e-6
                )
                assert _angle(
                    atoms[f"C{suffix}"], atoms[f"O{suffix}"], atoms[f"HO{suffix}"]
                ) == pytest.approx(108.5, abs=1.0e-6)
    contacts = modelled_hydroxyl_contacts(build.mapped, {"HO2": "O2", "HO3": "O3", "HO6": "O6"})
    assert contacts["hydrogen"] >= 1.9
    # Hydrogen-bond acceptors may approach to ~1.7 A, as in neutron structures.
    assert contacts["oxygen"] >= 1.6
    assert contacts["other_heavy"] >= 2.0
    assert contacts["hydrogen_bonded_fraction"] >= 0.95
    # The intrachain O3-H...O5 bond is present in every cellulose allomorph.
    assert contacts["o3_to_intrachain_o5_fraction"] == 1.0


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
def test_cellulose_ii_hydrogen_placement_scales_to_large_bulk() -> None:
    crystal = build_bulk_crystal(load_allomorph("ii"), (6, 5, 4))
    start = time.perf_counter()
    build = build_sage_system(crystal, "periodic")
    assert time.perf_counter() - start < 60.0
    contacts = modelled_hydroxyl_contacts(build.mapped, {"HO2": "O2", "HO3": "O3", "HO6": "O6"})
    assert contacts["hydrogen_bonded_fraction"] >= 0.95
    # The intrachain O3-H...O5 bond is present in every cellulose allomorph.
    assert contacts["o3_to_intrachain_o5_fraction"] == 1.0


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
def test_stereochemical_negative_controls_reject_inversions_and_mirror() -> None:
    molecule = build_sage_molecule(4, "periodic")
    # An idealized local tetrahedron is sufficient to prove that each control
    # reverses the independent signed-parity observable used at build time.
    positions = np.zeros((molecule.n_atoms, 3), dtype=float)
    for atom in molecule.atoms:
        if atom.name not in {"C1", "C3"} or int(atom.metadata["residue_number"]) != 2:
            continue
        centre = atom.molecule_atom_index
        positions[centre] = (4.0, 4.0, 4.0)
        neighbours = sorted(atom.bonded_atoms, key=lambda item: (item.symbol, item.name))[:3]
        offsets = np.eye(3)
        if int(np.sign(np.linalg.det(offsets))) != _REFERENCE_PARITY[atom.name]:
            offsets[[0, 1]] = offsets[[1, 0]]
        for neighbour, offset in zip(neighbours, offsets, strict=True):
            positions[neighbour.molecule_atom_index] = positions[centre] + offset
        original = _centre_parity(positions, atom, np.eye(3) * 20.0)
        inverted = positions.copy()
        left, right = (item.molecule_atom_index for item in neighbours[:2])
        inverted[[left, right]] = inverted[[right, left]]
        assert original == _REFERENCE_PARITY[atom.name]
        assert _centre_parity(inverted, atom, np.eye(3) * 20.0) == -original
        mirrored = positions.copy()
        mirrored[:, 0] *= -1
        assert _centre_parity(mirrored, atom, np.eye(3) * 20.0) == -original


@pytest.mark.parametrize(
    "arguments",
    (
        ["--shape", "slab", "--face", "1", "1", "0", "--layers", "6", "--chains-per-layer", "4"],
        ["--shape", "fibril", "--fibril-model", "18-chain"],
    ),
)
def test_openff_vacuum_shapes_warn_below_cutoff_vacuum(
    tmp_path: Path, arguments: list[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    def stop(*_: object, **__: object) -> None:
        raise SystemExit(0)

    # Stop before any OpenFF work; only the CLI vacuum policy is under test.
    monkeypatch.setattr(cli, "build_sage_system", stop)
    with pytest.warns(UserWarning, match="less than 9 angstrom of transverse vacuum"):
        with pytest.raises(SystemExit):
            main(
                [
                    "build",
                    "--allomorph",
                    "I-beta",
                    *arguments,
                    "--degree-polymerization",
                    "6",
                    "--chain-boundary",
                    "periodic",
                    "--vacuum",
                    "5",
                    "--force-field",
                    "openff",
                    "--md-engine",
                    "gromacs",
                    "--output",
                    str(tmp_path / "out"),
                ]
            )


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
def test_glycam_and_openff_place_identical_cellulose_ii_hydrogens() -> None:
    crystal = build_bulk_crystal(load_allomorph("ii"), (2, 2, 6))
    glycam = map_crystal_to_periodic_glycam(
        crystal, Path("examples/glycam/dp12_periodic_reference/cellulose_periodic_chain.pdb")
    )
    openff = build_sage_system(crystal, "periodic").mapped

    def hydrogens_by_oxygen(
        mapped: object, names: dict[str, str]
    ) -> dict[tuple[float, float, float], np.ndarray]:
        atoms = list(
            zip(
                mapped.atom_names,  # type: ignore[attr-defined]
                mapped.chain_numbers,  # type: ignore[attr-defined]
                mapped.residue_numbers,  # type: ignore[attr-defined]
                np.asarray(mapped.cartesian_coordinates),  # type: ignore[attr-defined]
                strict=True,
            )
        )
        oxygen = {(name, chain, residue): xyz for name, chain, residue, xyz in atoms}
        return {
            tuple(np.round(oxygen[(names[name], chain, residue)], 4)): xyz
            for name, chain, residue, xyz in atoms
            if name in names
        }

    from_glycam = hydrogens_by_oxygen(glycam, {"H2O": "O2", "H3O": "O3", "H6O": "O6"})
    from_openff = hydrogens_by_oxygen(openff, {"HO2": "O2", "HO3": "O3", "HO6": "O6"})
    assert from_glycam.keys() == from_openff.keys()
    for oxygen_position, hydrogen in from_glycam.items():
        np.testing.assert_allclose(hydrogen, from_openff[oxygen_position], atol=1e-6)


def _guard_case(polymer: str, allomorph: str, boundary: str, degree: int) -> tuple[object, object]:
    from cellulose_forge.polymer_mapping import (  # noqa: PLC0415
        chain_template,
        map_crystal_to_chain_template,
    )

    if allomorph == "ii":
        crystal = build_chain_parallel_slab(
            load_allomorph("ii"),
            MillerIndex(1, 1, 0),
            layers=3,
            chains_per_layer=4,
            degree_polymerization=degree,
            vacuum=20.0,
            chain_boundary=boundary,
        )
    else:
        crystal = build_chain_parallel_slab(
            load_allomorph(allomorph, polymer),
            MillerIndex(1, 0, 0),
            layers=1,
            chains_per_layer=2,
            degree_polymerization=degree,
            vacuum=20.0,
            chain_boundary=boundary,
        )
    molecule = build_sage_molecule(degree, boundary, polymer)
    template = chain_template(
        (atom.name, int(atom.metadata["residue_number"])) for atom in molecule.atoms
    )
    residue = "NAG" if polymer == "chitin" else "BGC"
    return molecule, map_crystal_to_chain_template(crystal, template, boundary, residue)


def _full_rdkit_label(
    molecule: object, coordinates: np.ndarray, box: np.ndarray, centre: int
) -> str | None:
    """Reference: RDKit's full perception with the chain unwrapped around ``centre``."""
    chem = importlib.import_module("rdkit.Chem")
    point_type = importlib.import_module("rdkit.Geometry").Point3D
    atoms = molecule.n_atoms  # type: ignore[attr-defined]
    adjacency: list[list[int]] = [[] for _ in range(atoms)]
    for bond in molecule.bonds:  # type: ignore[attr-defined]
        left, right = bond.atom1.molecule_atom_index, bond.atom2.molecule_atom_index
        adjacency[left].append(right)
        adjacency[right].append(left)
    inverse = np.linalg.inv(box)
    unwrapped = np.empty((atoms, 3))
    unwrapped[centre] = coordinates[centre]
    seen = {centre}
    queue = [centre]
    while queue:
        parent = queue.pop(0)
        for child in adjacency[parent]:
            if child not in seen:
                fractional = (coordinates[child] - coordinates[parent]) @ inverse
                unwrapped[child] = unwrapped[parent] + (fractional - np.round(fractional)) @ box
                seen.add(child)
                queue.append(child)
    rdkit_molecule = molecule.to_rdkit()  # type: ignore[attr-defined]
    conformer = chem.Conformer(atoms)
    for index, coordinate in enumerate(unwrapped):
        conformer.SetAtomPosition(index, point_type(*map(float, coordinate)))
    rdkit_molecule.RemoveAllConformers()
    rdkit_molecule.AddConformer(conformer)
    for atom in rdkit_molecule.GetAtoms():
        atom.SetChiralTag(chem.ChiralType.CHI_UNSPECIFIED)
        if atom.HasProp("_CIPCode"):
            atom.ClearProp("_CIPCode")
    chem.AssignStereochemistryFrom3D(rdkit_molecule, confId=0, replaceExistingTags=True)
    perceived = rdkit_molecule.GetAtomWithIdx(centre)
    return perceived.GetProp("_CIPCode") if perceived.HasProp("_CIPCode") else None


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
@pytest.mark.parametrize("degree", [4, 6])
@pytest.mark.parametrize("boundary", ["finite", "periodic"])
@pytest.mark.parametrize(
    ("polymer", "allomorph"), [("cellulose", "i_beta"), ("cellulose", "ii"), ("chitin", "beta")]
)
def test_fast_stereo_guard_agrees_with_full_rdkit_perception(
    polymer: str, allomorph: str, boundary: str, degree: int
) -> None:
    """The linear-time guard reads the same R/S as RDKit at every centre, wrapped or not."""
    molecule, mapped = _guard_case(polymer, allomorph, boundary, degree)
    fractional = mapped.cartesian_coordinates @ np.linalg.inv(mapped.basis_vectors)  # type: ignore[attr-defined]
    wrapped = replace(
        mapped,  # type: ignore[type-var]
        cartesian_coordinates=(fractional - np.floor(fractional)) @ mapped.basis_vectors,  # type: ignore[attr-defined]
    )
    for coordinates in (mapped, wrapped):
        assert validate_declared_stereochemistry(molecule, coordinates, boundary, polymer)["passed"]
    chain = wrapped.cartesian_coordinates[: molecule.n_atoms]  # type: ignore[attr-defined]
    for atom in molecule.atoms:  # type: ignore[attr-defined]
        if atom.stereochemistry is not None:
            perceived = _full_rdkit_label(
                molecule,
                chain,
                wrapped.basis_vectors,
                atom.molecule_atom_index,  # type: ignore[attr-defined]
            )
            assert perceived == atom.stereochemistry, (atom.metadata["residue_number"], atom.name)


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
@pytest.mark.parametrize("boundary", ["finite", "periodic"])
@pytest.mark.parametrize(("polymer", "allomorph"), [("cellulose", "i_beta"), ("chitin", "beta")])
def test_fast_stereo_guard_catches_every_single_label_flip(
    polymer: str, allomorph: str, boundary: str
) -> None:
    molecule, mapped = _guard_case(polymer, allomorph, boundary, 4)
    labelled = [atom for atom in molecule.atoms if atom.stereochemistry is not None]  # type: ignore[attr-defined]
    for atom in labelled:
        original = atom.stereochemistry
        atom.stereochemistry = "R" if original == "S" else "S"
        with pytest.raises(ValueError, match="declared stereochemistry disagrees"):
            validate_declared_stereochemistry(molecule, mapped, boundary, polymer)
        atom.stereochemistry = original
    assert len(labelled) == 20


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
def test_fast_stereo_guard_stays_fast_for_long_periodic_chains() -> None:
    """RDKit's global perception took seconds per DP20 periodic chain; the guard must not."""
    molecule, mapped = _guard_case("cellulose", "i_beta", "periodic", 40)
    validate_declared_stereochemistry(molecule, mapped, "periodic", "cellulose")
    started = time.perf_counter()
    result = validate_declared_stereochemistry(molecule, mapped, "periodic", "cellulose")
    assert time.perf_counter() - started < 1.0
    assert result["labelled_centres_per_chain"] == 200


@pytest.mark.skipif(find_spec("openff") is None, reason="optional OpenFF environment")
def test_declared_stereo_guard_requires_explicit_valid_boundary() -> None:
    molecule, mapped = _guard_case("cellulose", "i_beta", "finite", 4)
    with pytest.raises(ValueError, match=r"chain_boundary.*'finite' or 'periodic'"):
        validate_declared_stereochemistry(molecule, mapped, "inferred", "cellulose")

    # The supplied boundary controls the reference graph; it is not inferred
    # from terminal atom names.
    with pytest.raises(ValueError, match="molecular graph"):
        validate_declared_stereochemistry(molecule, mapped, "periodic", "cellulose")
