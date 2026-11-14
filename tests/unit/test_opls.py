"""Experimental OPLS-AA and OPLS-AA-SEI cellulose topologies and engine outputs."""

from __future__ import annotations

import json
import math
import shutil
import subprocess
from importlib.util import find_spec
from pathlib import Path

import numpy as np
import pytest

from cellulose_forge.cellulose_ii import cellulose_ii_preferred_acceptor
from cellulose_forge.cli import main
from cellulose_forge.crystallography import MillerIndex, load_allomorph
from cellulose_forge.geometry import (
    build_36_chain_fibril,
    build_bulk_crystal,
    build_chain_parallel_slab,
)
from cellulose_forge.opls import (
    _TORSIONS,
    KCAL_TO_KJ,
    SEI_ONE_FIVE_SCALE,
    SEI_ONE_SIX_SCALE,
    angle_parameters,
    bond_parameters,
    build_opls_system,
    build_opls_topology,
    nonbonded_exceptions,
    rb_coefficients,
)
from cellulose_forge.opls_engines import (
    SEI_TABLE_RANGE_ANGSTROM,
    _sei_table,
    write_opls_lammps,
)

DP = 6
# Per-residue SEI 1-5 pairs in a periodic chain (no chain ends).
SEI_ONE_FIVE_PER_RESIDUE = 47


@pytest.mark.parametrize("variant", ["opls-aa", "opls-aa-sei"])
@pytest.mark.parametrize("boundary", ["finite", "periodic"])
def test_chains_are_neutral_and_fully_parameterized(variant: str, boundary: str) -> None:
    topology = build_opls_topology(variant, DP, boundary)
    assert abs(sum(topology.charges)) < 1e-9
    for first, second in topology.bonds:
        bond_parameters(topology, first, second)
    for first, centre, last in topology.angles:
        angle_parameters(topology, first, centre, last)
    assert {name for *_, name in topology.dihedrals} <= set(_TORSIONS[variant])


def test_terminal_residue_charges_balance() -> None:
    topology = build_opls_topology("opls-aa", DP, "finite")
    by_residue: dict[int, float] = {}
    for atom, charge in zip(topology.atoms, topology.charges, strict=True):
        by_residue[atom.residue] = by_residue.get(atom.residue, 0.0) + charge
    assert by_residue[1] == pytest.approx(-0.2)  # reducing-end hemiacetal
    assert by_residue[DP] == pytest.approx(0.2)  # non-reducing-end C4 hydroxyl
    for residue in range(2, DP):
        assert by_residue[residue] == pytest.approx(0.0, abs=1e-12)


def test_periodic_chain_has_identical_residues() -> None:
    topology = build_opls_topology("opls-aa-sei", DP, "periodic")
    counts = topology.class_counts()
    assert all(count % DP == 0 for count in counts.values())
    one_five = [pair for pair in topology.scaled_pairs if pair[3] == "1-5"]
    assert len(one_five) == SEI_ONE_FIVE_PER_RESIDUE * DP
    assert not [pair for pair in topology.scaled_pairs if pair[3] == "1-6"]


def test_sei_scaling_is_intraresidue_and_one_six_only_at_hydroxyl_c4() -> None:
    topology = build_opls_topology("opls-aa-sei", DP, "finite")
    assert not build_opls_topology("opls-aa", DP, "finite").scaled_pairs
    for first, second, scale, kind in topology.scaled_pairs:
        assert topology.atoms[first].residue == topology.atoms[second].residue
        divisor = SEI_ONE_FIVE_SCALE if kind == "1-5" else SEI_ONE_SIX_SCALE
        assert scale == pytest.approx(1.0 / divisor)
    one_six = {
        frozenset((topology.atoms[a].name, topology.atoms[b].name))
        for a, b, _, kind in topology.scaled_pairs
        if kind == "1-6"
    }
    assert one_six == {frozenset({"HO6", "O4"}), frozenset({"O6", "HO4"})}
    residues = {
        topology.atoms[a].residue for a, _, _, kind in topology.scaled_pairs if kind == "1-6"
    }
    assert residues == {DP}


def test_every_pair_has_one_nonbonded_treatment() -> None:
    topology = build_opls_topology("opls-aa-sei", DP, "finite")
    entries = nonbonded_exceptions(topology)
    kinds = {kind for *_, kind in entries}
    assert kinds == {"excluded", "1-4", "1-5", "1-6"}
    assert sum(1 for *_, kind in entries if kind == "1-4") == len(topology.pairs_14)


@pytest.mark.parametrize("phi", np.linspace(-math.pi, math.pi, 13))
def test_rb_conversion_reproduces_the_fourier_series(phi: float) -> None:
    v1, v2, v3 = (1.501, -0.729, 0.870)
    fourier = KCAL_TO_KJ * (
        0.5 * v1 * (1 + math.cos(phi))
        + 0.5 * v2 * (1 - math.cos(2 * phi))
        + 0.5 * v3 * (1 + math.cos(3 * phi))
    )
    psi = phi - math.pi
    rb = sum(c * math.cos(psi) ** n for n, c in enumerate(rb_coefficients((v1, v2, v3))))
    assert rb == pytest.approx(fourier, abs=1e-9)


def test_kony_table_values_are_transcribed() -> None:
    assert _TORSIONS["opls-aa"]["OH-C-C-OH"] == (9.066, 0.0, 0.0)
    assert _TORSIONS["opls-aa-sei"]["OH-C-C-OH"] == (8.519, -6.829, 4.793)
    assert _TORSIONS["opls-aa-sei"]["C-C-C-OS_ace"] == (1.501, 0.0, 0.0)
    assert _TORSIONS["opls-aa-sei"]["H-C-C-H"] == _TORSIONS["opls-aa"]["H-C-C-H"]


def test_sei_table_is_full_lennard_jones_plus_scaled_coulomb() -> None:
    sigma, epsilon, charge_product, scale = 0.3, 0.5, -0.28, 1 / 1.26
    rows = [line.split() for line in _sei_table(sigma, epsilon, charge_product, scale)]
    assert float(rows[0][1]) == SEI_TABLE_RANGE_ANGSTROM[0]
    def energy_at(r: float) -> float:
        s6 = (10 * sigma / r) ** 6
        lj = 4 * epsilon / KCAL_TO_KJ * (s6 * s6 - s6)
        return lj + scale * 332.06371 * charge_product / r

    r, energy, force = (float(value) for value in rows[2000][1:])
    assert energy == pytest.approx(energy_at(r), rel=1e-10)
    step = 1e-6
    derivative = (energy_at(r + step) - energy_at(r - step)) / (2 * step)
    assert force == pytest.approx(-derivative, rel=1e-6)


def _bulk_build(variant: str):  # type: ignore[no-untyped-def]
    crystal = build_bulk_crystal(load_allomorph("ii"), (3, 3, 3))
    return build_opls_system(crystal, variant, "periodic")


@pytest.mark.parametrize("shape", ["slab", "fibril"])
def test_cellulose_ii_surface_table_3_acceptors_take_priority(shape: str) -> None:
    """A surface fallback must not occupy an available Table-3 acceptor."""
    definition = load_allomorph("ii")
    if shape == "slab":
        crystal = build_chain_parallel_slab(
            definition,
            MillerIndex(1, 1, 0),
            layers=3,
            chains_per_layer=4,
            degree_polymerization=4,
            vacuum=30.0,
            chain_boundary="periodic",
        )
    else:
        crystal = build_36_chain_fibril(
            definition,
            degree_polymerization=4,
            vacuum=30.0,
            chain_boundary="periodic",
        )
    mapped = build_opls_system(crystal, "opls-aa", "periodic").mapped
    labels = np.asarray(crystal.atom_labels)
    crystal_chains = np.asarray(crystal.chain_numbers)
    crystal_residues = np.asarray(crystal.residue_numbers)
    crystal_xyz = np.asarray(crystal.cartesian_coordinates)
    names = np.asarray(mapped.atom_names)
    chains = np.asarray(mapped.chain_numbers)
    residues = np.asarray(mapped.residue_numbers)
    xyz = np.asarray(mapped.cartesian_coordinates)
    inverse = np.linalg.inv(mapped.basis_vectors)

    def minimum_image(delta: np.ndarray) -> np.ndarray:
        fractional = delta @ inverse
        fractional -= np.round(fractional)
        return fractional @ mapped.basis_vectors

    checked = 0
    for chain in sorted(set(crystal.chain_numbers)):
        for residue in sorted(set(crystal_residues[crystal_chains == chain])):
            for donor, hydrogen in (("O2", "HO2"), ("O6", "HO6")):
                donor_coordinate = crystal_xyz[
                    (crystal_chains == chain)
                    & (crystal_residues == residue)
                    & (labels == donor)
                ][0]
                acceptor = cellulose_ii_preferred_acceptor(
                    crystal,
                    chain_number=chain,
                    residue_number=residue,
                    donor=donor,
                )
                if np.linalg.norm(minimum_image(acceptor - donor_coordinate)) > 3.2:
                    continue
                # Origin-chain residue numbering differs between the crystal
                # and mapped force-field system. Locate the mapped donor by
                # position rather than assuming those residue IDs coincide.
                donor_candidates = np.flatnonzero((chains == chain) & (names == donor))
                offsets = np.linalg.norm(
                    np.vstack(
                        [
                            minimum_image(xyz[index] - donor_coordinate)
                            for index in donor_candidates
                        ]
                    ),
                    axis=1,
                )
                mapped_donor = donor_candidates[int(np.argmin(offsets))]
                assert offsets.min() < 1.0e-3
                mapped_residue = residues[mapped_donor]
                hydrogen_coordinate = xyz[
                    (chains == chain) & (residues == mapped_residue) & (names == hydrogen)
                ][0]
                delta = minimum_image(acceptor - hydrogen_coordinate)
                assert np.linalg.norm(delta) < 2.6
                to_donor = minimum_image(donor_coordinate - hydrogen_coordinate)
                angle = np.degrees(
                    np.arccos(
                        np.clip(
                            np.dot(to_donor, delta)
                            / (np.linalg.norm(to_donor) * np.linalg.norm(delta)),
                            -1.0,
                            1.0,
                        )
                    )
                )
                assert angle > 120.0
                checked += 1
    assert checked > 0


@pytest.mark.parametrize("variant", ["opls-aa", "opls-aa-sei"])
def test_lammps_files_use_opls_styles(tmp_path: Path, variant: str) -> None:
    build = _bulk_build(variant)
    files, metadata = write_opls_lammps(build, tmp_path, dispersion_correction=True)
    text = files.input.read_text()
    assert "pair_modify mix geometric tail yes" in text
    assert "special_bonds lj/coul 0.0 0.0 0.5 angle yes dihedral yes" in text
    assert "dihedral_style opls" in text
    data = files.data.read_text()
    chains = build.mapped.chain_count
    assert f"{chains * len(build.topology.angles)} angles" in data
    if variant == "opls-aa-sei":
        assert files.table is not None and files.table.is_file()
        assert "bond_style hybrid harmonic table" in text
        bonds = len(build.topology.bonds) + len(build.topology.scaled_pairs)
        assert f"{chains * bonds} bonds" in data
        assert metadata["sei_pseudo_bond_count"] == chains * len(build.topology.scaled_pairs)
    else:
        assert files.table is None
        assert "bond_style harmonic" in text


@pytest.mark.skipif(shutil.which("lmp_serial") is None, reason="LAMMPS not installed")
def test_lammps_sei_energy_equals_plain_plus_scaled_pairs(tmp_path: Path) -> None:
    """The SEI pseudo-bonds change only the scaled Coulomb terms and torsions."""
    energies = {}
    for variant in ("opls-aa", "opls-aa-sei"):
        directory = tmp_path / variant
        write_opls_lammps(_bulk_build(variant), directory, dispersion_correction=False)
        output = subprocess.run(
            ["lmp_serial", "-in", "in.cellulose", "-log", "none"],
            cwd=directory, capture_output=True, text=True, check=True,
        ).stdout.splitlines()
        header = next(i for i, line in enumerate(output) if line.split()[:2] == ["Step", "Temp"])
        energies[variant] = dict(
            zip(output[header].split(), map(float, output[header + 1].split()), strict=True)
        )
    for term in ("E_angle", "E_long"):
        assert energies["opls-aa"][term] == pytest.approx(energies["opls-aa-sei"][term])


def _cli(tmp_path: Path, *extra: str) -> list[str]:
    return [
        "build", "--allomorph", "II", "--shape", "bulk", "--cells", "3", "3", "3",
        "--chain-boundary", "periodic", "--output", str(tmp_path / "out"), *extra,
    ]


def test_cli_writes_gromacs_and_lammps(tmp_path: Path) -> None:
    with pytest.warns(UserWarning, match="experimental"):
        main(_cli(tmp_path, "--force-field", "opls-aa-sei", "--md-engine", "gromacs",
                  "--md-engine", "lammps"))
    report = json.loads((tmp_path / "out" / "build_report.json").read_text())
    assert report["outputs"]["lammps_sei_table"] == "sei_pairs.table"
    assert report["outputs"]["gromacs_topology"] == "system.top"
    assert report["dispersion_correction"] is True
    assert "lammps" in report["engine_validation"]


@pytest.mark.skipif(find_spec("openmm") is None, reason="optional OpenMM environment")
def test_cli_writes_openmm(tmp_path: Path) -> None:
    with pytest.warns(UserWarning, match="experimental"):
        main(_cli(tmp_path, "--force-field", "opls-aa", "--md-engine", "openmm"))
    report = json.loads((tmp_path / "out" / "build_report.json").read_text())
    assert report["engine_validation"]["openmm"]["status"] == "passed"
    assert not (tmp_path / "out" / "system.top").exists()


@pytest.mark.parametrize(
    ("extra", "message"),
    [
        (["--md-engine", "amber"], "support --md-engine gromacs, openmm, and lammps"),
        (["--restrain-chains", "1"], "position restraints are not implemented"),
    ],
)
def test_cli_rejects_unsupported_requests(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], extra: list[str], message: str
) -> None:
    with pytest.raises(SystemExit):
        main(_cli(tmp_path, "--force-field", "opls-aa", *extra))
    assert message in capsys.readouterr().err


def test_cli_rejects_unit_cells(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit):
        main(["build", "--allomorph", "I-beta", "--shape", "unit-cell", "--force-field",
              "opls-aa", "--output", str(tmp_path / "out")])
    # Unit cells fail the box rule before the OPLS shape check; either refusal is
    # correct, and neither writes output.
    error = capsys.readouterr().err
    assert "periodic box must exceed" in error or "periodic bulk crystals" in error
    assert not (tmp_path / "out").exists()
