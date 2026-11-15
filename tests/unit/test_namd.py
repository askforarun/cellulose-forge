"""Tests for native NAMD/X-PLOR PSF output."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from cellulose_forge.charmm36 import (
    Charmm36TemplateAtom,
    FiniteCharmm36Topology,
    _IndexedInteraction,
)
from cellulose_forge.cli import main
from cellulose_forge.glycam_mapping import GlycamMappedSystem
from cellulose_forge.namd import (
    render_namd_configuration,
    render_pdb,
    render_psf,
    render_t0_comparison_script,
    validate_charmm_parameter_file,
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


def _cli_arguments(output: Path, *extra: str) -> list[str]:
    return [
        "build",
        "--polymer",
        "cellulose",
        "--allomorph",
        "I-beta",
        "--shape",
        "slab",
        "--face",
        "1",
        "1",
        "0",
        "--layers",
        "2",
        "--chains-per-layer",
        "4",
        "--degree-polymerization",
        "6",
        "--chain-boundary",
        "periodic",
        "--vacuum",
        "30",
        "--force-field",
        "charmm36",
        *extra,
        "--output",
        str(output),
    ]


def test_psf_expands_chains_and_writes_zero_exclusion_array() -> None:
    text = render_psf(_system(), _topology())

    assert "         8 !NATOM" in text
    assert "         6 !NBOND: bonds" in text
    assert "         4 !NTHETA: angles" in text
    assert "         2 !NPHI: dihedrals" in text
    assert "         2 !NIMPHI: impropers" in text
    assert "C001" in text
    assert "C002" in text
    assert "!NUMLP NUMLPH" in text
    assert "!NCRTERM: cross-terms" in text
    angle_rows = text.split("!NTHETA: angles\n", maxsplit=1)[1].split("\n\n", maxsplit=1)[0]
    assert all(len(line.split()) % 3 == 0 for line in angle_rows.splitlines() if line.strip())
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
    assert pdb.startswith("CRYST1")
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
    assert configuration.index("minimize 1000") < configuration.index("reinitvels 300")
    assert configuration.index("run 0") < configuration.index("minimize 1000")


def test_pdb_uses_fixed_width_alias_for_six_character_chitin_residue() -> None:
    system = replace(_system(), residue_names=("BGLCNA",) * 8)

    atom_line = next(
        line for line in render_pdb(system, _topology()).splitlines() if line.startswith("ATOM")
    )

    assert atom_line[17:21].strip() == "NAG"
    assert float(atom_line[30:38]) == 0.0


def test_parameter_file_checksum_is_enforced(tmp_path: Path) -> None:
    parameters = tmp_path / "wrong.prm"
    parameters.write_text("* wrong release\n", encoding="utf-8")

    with pytest.raises(ValueError, match="expected SHA-256"):
        validate_charmm_parameter_file(parameters)


def test_cli_rejects_namd_inputs_before_creating_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    output = tmp_path / "out"
    missing = tmp_path / "missing.prm"
    with pytest.raises(SystemExit, match="2"):
        main(
            _cli_arguments(
                output,
                "--md-engine",
                "namd",
                "--charmm-parameter-file",
                str(missing),
            )
        )
    assert "requires an existing CHARMM parameter file" in capsys.readouterr().err
    assert not output.exists()

    parameters = tmp_path / "parameters.prm"
    parameters.write_text("* test parameters\n", encoding="utf-8")
    checksum = hashlib.sha256(parameters.read_bytes()).hexdigest()
    monkeypatch.setattr("cellulose_forge.namd.CHARMM36_CARB_PARAMETER_SHA256", checksum)
    with pytest.raises(SystemExit, match="2"):
        main(
            _cli_arguments(
                output,
                "--md-engine",
                "namd",
                "--charmm-parameter-file",
                str(parameters),
                "--restrain-layers",
                "1",
            )
        )
    assert "NAMD position restraints are not implemented" in capsys.readouterr().err
    assert not output.exists()


def test_cli_rejects_unused_parameter_file(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    parameters = tmp_path / "parameters.prm"
    parameters.write_text("* test parameters\n", encoding="utf-8")

    with pytest.raises(SystemExit, match="2"):
        main(_cli_arguments(tmp_path / "out", "--charmm-parameter-file", str(parameters)))

    assert "requires --md-engine namd" in capsys.readouterr().err


def test_cli_applies_namd_twelve_angstrom_cutoff(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    parameters = tmp_path / "parameters.prm"
    parameters.write_text("* test parameters\n", encoding="utf-8")
    checksum = hashlib.sha256(parameters.read_bytes()).hexdigest()
    monkeypatch.setattr("cellulose_forge.namd.CHARMM36_CARB_PARAMETER_SHA256", checksum)

    with pytest.raises(SystemExit, match="2"):
        main(
            _cli_arguments(
                tmp_path / "out",
                "--md-engine",
                "namd",
                "--charmm-parameter-file",
                str(parameters),
            )
        )

    assert "must exceed 24 angstrom" in capsys.readouterr().err
    assert not (tmp_path / "out").exists()


def test_writer_copies_parameter_file_and_writes_complete_bundle(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    parameters = tmp_path / "source.prm"
    parameters.write_text("* test parameters\n", encoding="utf-8")
    checksum = hashlib.sha256(parameters.read_bytes()).hexdigest()
    monkeypatch.setattr("cellulose_forge.namd.CHARMM36_CARB_PARAMETER_SHA256", checksum)

    files = write_namd_system(
        _system(), _topology(), parameters, tmp_path / "namd", periodic_chains=False
    )

    assert files.psf.is_file()
    assert files.pdb.is_file()
    assert files.configuration.is_file()
    assert files.parameters.read_text(encoding="utf-8") == "* test parameters\n"
    assert "vdw-modifier = Force-switch" in files.gromacs_reference.read_text()
    assert "periodic-molecules = no" in files.gromacs_reference.read_text()
    assert "gromacs_t0_energy.xvg" in files.slurm.read_text()
    comparison = files.energy_comparison.read_text(encoding="utf-8")
    compile(comparison, str(files.energy_comparison), "exec")


def test_t0_comparison_script_groups_matching_terms(tmp_path: Path) -> None:
    script = tmp_path / "compare.py"
    script.write_text(render_t0_comparison_script(), encoding="utf-8")
    namd = tmp_path / "namd.log"
    namd.write_text(
        "ENERGY: 0 1 5 3 4 26 19 0 0 0 58 300 58 58 300 0 0 1000 0 0\n",
        encoding="utf-8",
    )
    xvg = tmp_path / "energy.xvg"
    names = (
        "Bond",
        "U-B",
        "Angle",
        "Proper Dih.",
        "Improper Dih.",
        "LJ-14",
        "Coulomb-14",
        "LJ (SR)",
        "Coulomb (SR)",
        "Coul. recip.",
        "Potential",
    )
    values_kcal = (1, 2, 3, 3, 4, 7, 8, 12, 9, 9, 58)
    xvg.write_text(
        "\n".join(
            [
                *(f'@ s{index} legend "{name}"' for index, name in enumerate(names)),
                "0 " + " ".join(str(value * 4.184) for value in values_kcal),
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    completed = subprocess.run(
        [sys.executable, str(script), str(namd), str(xvg)],
        check=True,
        capture_output=True,
        text=True,
    )
    result = json.loads(completed.stdout)

    for term in ("BOND", "ANGLE", "DIHED", "IMPRP", "ELECT", "VDW", "POTENTIAL"):
        assert result["comparison"][term]["difference_kcal_mol"] == pytest.approx(0.0)
