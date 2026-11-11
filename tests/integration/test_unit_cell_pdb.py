"""Integration tests for unit-cell PDB generation through the CLI."""

import json
import shutil
from pathlib import Path

import pytest

from cellulose_forge.cli import main
from cellulose_forge.glycam import GlycamChainFiles


def _use_committed_glycam_reference(
    monkeypatch: pytest.MonkeyPatch,
    *,
    chain_boundary: str,
) -> None:
    reference_name = "dp12_periodic_reference" if chain_boundary == "periodic" else "dp12_reference"
    source_prefix = (
        "cellulose_periodic_chain" if chain_boundary == "periodic" else "cellulose_chain"
    )

    def prepare(**kwargs: object) -> tuple[GlycamChainFiles, dict[str, object]]:
        directory = Path(str(kwargs["directory"]))
        files = GlycamChainFiles(directory, source_prefix)
        source = Path("examples/glycam") / reference_name
        shutil.copyfile(source / f"{source_prefix}.pdb", files.pdb)
        shutil.copyfile(source / f"{source_prefix}.top", files.gromacs_topology)
        return files, {"test_fixture": reference_name}

    monkeypatch.setattr("cellulose_forge.cli.generate_glycam_gromacs_reference", prepare)


def test_cli_writes_annotated_i_beta_unit_cell(tmp_path: Path) -> None:
    output = tmp_path / "i_beta_unit_cell.pdb"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "unit-cell",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("HETATM")]
    assert len(atom_lines) == 84
    assert any(line.startswith("CRYST1") for line in lines)
    assert {line[21] for line in atom_lines} == {"A", "B"}
    assert {(line[21], int(line[22:26])) for line in atom_lines} == {
        ("A", 1),
        ("A", 2),
        ("B", 1),
        ("B", 2),
    }
    report = json.loads(
        (tmp_path / "i_beta_unit_cell.build_report.json").read_text(encoding="utf-8")
    )
    assert report["force_field"] == "none"
    assert report["md_engines"] == []
    assert report["degree_polymerization"] == 2
    assert report["chain_count"] == 2
    assert report["outputs"]["coordinates"] == output.name


def test_cli_writes_i_beta_bulk_supercell(tmp_path: Path) -> None:
    output = tmp_path / "i_beta_2x2x2.pdb"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "bulk",
            "--cells",
            "2",
            "2",
            "2",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("HETATM")]
    cryst1 = next(line for line in lines if line.startswith("CRYST1"))
    assert len(atom_lines) == 84 * 2 * 2 * 2
    assert (float(cryst1[6:15]), float(cryst1[15:24]), float(cryst1[24:33])) == (
        15.568,
        16.402,
        20.760,
    )
    assert {line[21] for line in atom_lines} == set("ABCDEFGH")
    assert {int(line[22:26]) for line in atom_lines} == {1, 2, 3, 4}


def test_cli_writes_i_alpha_unit_cell(tmp_path: Path) -> None:
    output = tmp_path / "i_alpha_unit_cell.pdb"

    result = main(
        [
            "build",
            "--allomorph",
            "I-alpha",
            "--shape",
            "unit-cell",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("HETATM")]
    cryst1 = next(line for line in lines if line.startswith("CRYST1"))
    assert len(atom_lines) == 42
    assert {line[21] for line in atom_lines} == {"A"}
    assert {int(line[22:26]) for line in atom_lines} == {1, 2}
    assert (float(cryst1[6:15]), float(cryst1[15:24]), float(cryst1[24:33])) == (
        10.400,
        6.717,
        5.962,
    )


def test_cli_writes_iii_i_unit_cell(tmp_path: Path) -> None:
    output = tmp_path / "iii_i_unit_cell.pdb"

    assert (
        main(
            [
                "build",
                "--allomorph",
                "III-I",
                "--shape",
                "unit-cell",
                "--output",
                str(output),
            ]
        )
        == 0
    )
    lines = output.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("HETATM")]
    assert len(atom_lines) == 42
    assert {line[21] for line in atom_lines} == {"A"}
    assert {int(line[22:26]) for line in atom_lines} == {1, 2}


def test_cli_writes_ii_unit_cell_without_inventing_hydrogens(tmp_path: Path) -> None:
    output = tmp_path / "ii_unit_cell.pdb"

    assert (
        main(
            [
                "build",
                "--allomorph",
                "II",
                "--shape",
                "unit-cell",
                "--output",
                str(output),
            ]
        )
        == 0
    )
    lines = output.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("HETATM")]
    assert len(atom_lines) == 72
    assert {line[21] for line in atom_lines} == {"A", "B"}
    assert not any(line[12:16].strip() in {"HO2", "HO3", "HO6"} for line in atom_lines)


def test_cli_writes_i_beta_110_slab(tmp_path: Path) -> None:
    output = tmp_path / "i_beta_110_slab.pdb"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "slab",
            "--face",
            "1",
            "1",
            "0",
            "--layers",
            "4",
            "--chains-per-layer",
            "20",
            "--degree-polymerization",
            "12",
            "--chain-boundary",
            "finite",
            "--vacuum",
            "30",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("HETATM")]
    assert len(atom_lines) == 20_160
    assert any("shape=slab face=(1,1,0) layers=4" in line for line in lines)
    assert any("chains_per_layer=20 total_chains=80" in line for line in lines)
    assert any("degree_polymerization=12" in line for line in lines)
    assert any("vacuum=30.000000 angstrom" in line for line in lines)
    assert atom_lines[0][72:76] == "C001"
    assert {line[72:76] for line in atom_lines} == {f"C{chain:03d}" for chain in range(1, 81)}


def test_cli_writes_slab_from_physical_dimensions(tmp_path: Path) -> None:
    output = tmp_path / "i_beta_110_sized_slab.pdb"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "slab",
            "--face",
            "1",
            "1",
            "0",
            "--slab-width",
            "40",
            "--slab-thickness",
            "18",
            "--slab-length",
            "25",
            "--chain-boundary",
            "periodic",
            "--vacuum",
            "10",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    assert any("size_mode=physical_dimensions" in line for line in lines)
    assert any("requested_width=40.000000" in line for line in lines)
    assert any("requested_thickness=18.000000" in line for line in lines)
    assert any("requested_length=25.000000" in line for line in lines)
    assert any("chain_boundary=periodic" in line for line in lines)


def test_cli_writes_slab_with_mixed_axis_sizing(tmp_path: Path) -> None:
    output = tmp_path / "i_beta_110_mixed_slab.pdb"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "slab",
            "--face",
            "1",
            "1",
            "0",
            "--layers",
            "4",
            "--slab-width",
            "40",
            "--slab-length",
            "40",
            "--chain-boundary",
            "finite",
            "--vacuum",
            "10",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    assert any("size_mode=mixed" in line for line in lines)
    assert any("shape=slab face=(1,1,0) layers=4" in line for line in lines)
    assert any("chains_per_layer=8 total_chains=32" in line for line in lines)
    assert any("degree_polymerization=8" in line for line in lines)
    report = json.loads(
        (tmp_path / "i_beta_110_mixed_slab.build_report.json").read_text(encoding="utf-8")
    )
    assert report["degree_polymerization"] == 8
    assert report["chain_count"] == 32


def test_cli_rejects_two_sizing_options_for_same_slab_axis(tmp_path: Path) -> None:
    with pytest.raises(SystemExit, match="2"):
        main(
            [
                "build",
                "--allomorph",
                "I-beta",
                "--shape",
                "slab",
                "--face",
                "1",
                "1",
                "0",
                "--layers",
                "4",
                "--slab-thickness",
                "20",
                "--slab-width",
                "40",
                "--slab-length",
                "40",
                "--chain-boundary",
                "finite",
                "--output",
                str(tmp_path / "invalid.pdb"),
            ]
        )


def test_cli_writes_literature_i_beta_36_chain_fibril(tmp_path: Path) -> None:
    output = tmp_path / "i_beta_36_chain_fibril.pdb"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "fibril",
            "--fibril-model",
            "36-chain",
            "--degree-polymerization",
            "6",
            "--chain-boundary",
            "periodic",
            "--vacuum",
            "20",
            "--padding-x",
            "24",
            "--padding-y",
            "32",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("HETATM")]
    assert len(atom_lines) == 36 * 6 * 21
    assert {line[72:76] for line in atom_lines} == {f"C{chain:03d}" for chain in range(1, 37)}
    assert any("shape=fibril fibril_model=36-chain total_chains=36" in line for line in lines)
    assert any("padding_x=24.000000" in line for line in lines)
    assert any("padding_y=32.000000" in line for line in lines)


@pytest.mark.parametrize(
    ("model", "chain_count"),
    [("18-chain", 18), ("24-chain", 24)],
)
def test_cli_writes_named_i_beta_fibril_profiles(
    tmp_path: Path,
    model: str,
    chain_count: int,
) -> None:
    output = tmp_path / f"i_beta_{model}.pdb"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "fibril",
            "--fibril-model",
            model,
            "--degree-polymerization",
            "4",
            "--chain-boundary",
            "periodic",
            "--vacuum",
            "20",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    assert len([line for line in lines if line.startswith("HETATM")]) == (
        chain_count * 4 * 21
    )
    assert any(f"fibril_model={model}" in line for line in lines)


def test_cli_writes_custom_i_beta_fibril_profile(tmp_path: Path) -> None:
    output = tmp_path / "i_beta_custom_fibril.pdb"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "fibril",
            "--fibril-profile",
            "1,2,3,2,1",
            "--degree-polymerization",
            "4",
            "--chain-boundary",
            "periodic",
            "--vacuum",
            "20",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    assert len([line for line in lines if line.startswith("HETATM")]) == 9 * 4 * 21
    assert any("fibril_profile_bottom_to_top=1,2,3,2,1" in line for line in lines)


def test_cli_rejects_named_and_custom_fibril_profiles_together(tmp_path: Path) -> None:
    with pytest.raises(SystemExit, match="2"):
        main(
            [
                "build",
                "--allomorph",
                "I-beta",
                "--shape",
                "fibril",
                "--fibril-model",
                "18-chain",
                "--fibril-profile",
                "2,3,4,4,3,2",
                "--degree-polymerization",
                "4",
                "--chain-boundary",
                "periodic",
                "--output",
                str(tmp_path / "invalid.pdb"),
            ]
        )


def test_cli_rejects_native_amber_for_periodic_chains(tmp_path: Path) -> None:
    with pytest.raises(SystemExit, match="2"):
        main(
            [
                "build",
                "--allomorph",
                "I-beta",
                "--shape",
                "fibril",
                "--fibril-model",
                "18-chain",
                "--degree-polymerization",
                "4",
                "--chain-boundary",
                "periodic",
                "--force-field",
                "glycam",
                "--md-engine",
                "amber",
                "--output",
                str(tmp_path / "unsafe_periodic_amber"),
            ]
        )


def test_cli_rejects_md_engine_without_force_field(tmp_path: Path) -> None:
    with pytest.raises(SystemExit, match="2"):
        main(
            [
                "build",
                "--allomorph",
                "I-beta",
                "--shape",
                "unit-cell",
                "--md-engine",
                "gromacs",
                "--output",
                str(tmp_path / "invalid.pdb"),
            ]
        )


def test_cli_rejects_periodic_native_openmm_engine(tmp_path: Path) -> None:
    with pytest.raises(SystemExit, match="2"):
        main(
            [
                "build",
                "--allomorph",
                "I-beta",
                "--shape",
                "unit-cell",
                "--chain-boundary",
                "periodic",
                "--force-field",
                "glycam",
                "--md-engine",
                "openmm",
                "--output",
                str(tmp_path / "invalid"),
            ]
        )


def test_cli_rejects_md_engine_and_legacy_topology_format(tmp_path: Path) -> None:
    with pytest.raises(SystemExit, match="2"):
        main(
            [
                "build",
                "--allomorph",
                "I-beta",
                "--shape",
                "unit-cell",
                "--chain-boundary",
                "periodic",
                "--force-field",
                "glycam",
                "--md-engine",
                "gromacs",
                "--topology-format",
                "gromacs",
                "--output",
                str(tmp_path / "invalid"),
            ]
        )


def test_cli_writes_finite_fibril_with_longitudinal_padding(tmp_path: Path) -> None:
    output = tmp_path / "i_beta_finite_padded_fibril.pdb"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "fibril",
            "--fibril-model",
            "36-chain",
            "--degree-polymerization",
            "6",
            "--chain-boundary",
            "finite",
            "--padding-x",
            "20",
            "--padding-y",
            "30",
            "--padding-longitudinal",
            "40",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    lines = output.read_text(encoding="utf-8").splitlines()
    assert any("padding_z=40.000000" in line for line in lines)


def test_cli_rejects_longitudinal_padding_for_periodic_fibril(tmp_path: Path) -> None:
    with pytest.raises(SystemExit, match="2"):
        main(
            [
                "build",
                "--allomorph",
                "I-beta",
                "--shape",
                "fibril",
                "--fibril-model",
                "36-chain",
                "--degree-polymerization",
                "6",
                "--chain-boundary",
                "periodic",
                "--padding-longitudinal",
                "20",
                "--output",
                str(tmp_path / "invalid_periodic_fibril.pdb"),
            ]
        )


def test_cli_writes_periodic_glycam_fibril_system(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _use_committed_glycam_reference(monkeypatch, chain_boundary="periodic")
    output = tmp_path / "ibeta_fibril"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "fibril",
            "--fibril-model",
            "36-chain",
            "--degree-polymerization",
            "12",
            "--chain-boundary",
            "periodic",
            "--vacuum",
            "30",
            "--force-field",
            "glycam",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    assert (output / "structure.pdb").is_file()
    assert (output / "structure.gro").is_file()
    topology = (output / "system.top").read_text(encoding="utf-8")
    assert "CELLULOSE_CHAIN  36" in topology
    report = json.loads((output / "build_report.json").read_text(encoding="utf-8"))
    assert report["force_field"] == "GLYCAM"
    assert report["md_engines"] == ["gromacs"]
    assert report["chain_boundary"] == "periodic"
    assert report["degree_polymerization"] == 12
    assert report["chain_count"] == 36
    assert report["force_field_atom_count"] == 9072
    chirality = report["chirality_validation"]
    assert chirality["status"] == "passed"
    assert chirality["configuration"] == "D-cellulose"
    assert chirality["residue_count"] == 432


def test_cli_writes_selectively_restrained_fibril_chains(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _use_committed_glycam_reference(monkeypatch, chain_boundary="periodic")
    output = tmp_path / "restrained_fibril"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "fibril",
            "--fibril-model",
            "36-chain",
            "--degree-polymerization",
            "12",
            "--chain-boundary",
            "periodic",
            "--force-field",
            "glycam",
            "--restrain-chains",
            "6",
            "11",
            "16",
            "21",
            "26",
            "31",
            "--restraint-force",
            "1000",
            "1000",
            "1000",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    topology = (output / "system.top").read_text(encoding="utf-8")
    assert "CELLULOSE_CHAIN  30" in topology
    assert "CELLULOSE_CHAIN_RESTRAINED  6" in topology
    assert "#ifdef POSRES_CELLULOSE" in topology
    report = json.loads((output / "build_report.json").read_text(encoding="utf-8"))
    restraints = report["position_restraints"]
    assert restraints["selected_chain_ids"] == [6, 11, 16, 21, 26, 31]
    assert restraints["force_constants_kj_mol_nm2"] == [1000.0, 1000.0, 1000.0]
    assert restraints["atom_selection"] == "heavy atoms in each selected chain"
    assert [
        entry["original_chain_id"]
        for entry in restraints["chain_mapping"]
        if entry["restraint_status"] == "restrained"
    ] == [6, 11, 16, 21, 26, 31]


def test_cli_rejects_combined_restraint_selectors(tmp_path: Path) -> None:
    with pytest.raises(SystemExit, match="2"):
        main(
            [
                "build",
                "--allomorph",
                "I-beta",
                "--shape",
                "fibril",
                "--fibril-model",
                "36-chain",
                "--degree-polymerization",
                "12",
                "--restrain-chains",
                "1",
                "--restrain-layers",
                "1",
                "--restraint-force",
                "1000",
                "1000",
                "1000",
                "--output",
                str(tmp_path / "invalid_restraints"),
            ]
        )


def test_cli_writes_finite_glycam_slab_system(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _use_committed_glycam_reference(monkeypatch, chain_boundary="finite")
    output = tmp_path / "ibeta_slab"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "slab",
            "--face",
            "1",
            "1",
            "0",
            "--layers",
            "1",
            "--chains-per-layer",
            "2",
            "--degree-polymerization",
            "12",
            "--chain-boundary",
            "finite",
            "--vacuum",
            "30",
            "--padding-longitudinal",
            "12",
            "--force-field",
            "glycam",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    report = json.loads((output / "build_report.json").read_text(encoding="utf-8"))
    assert report["chain_boundary"] == "finite"
    assert report["chain_count"] == 2
    assert report["force_field_atom_count"] == 510
    assert report["closest_longitudinal_periodic_image_contact"][
        "distance_angstrom"
    ] >= 1.0
    assert "CELLULOSE_CHAIN  2" in (output / "system.top").read_text(encoding="utf-8")


def test_cli_writes_and_validates_finite_openmm_system(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _use_committed_glycam_reference(monkeypatch, chain_boundary="finite")

    def write_amber(
        system: object,
        template_prmtop: object,
        template_inpcrd: object,
        output_directory: str | Path,
        **kwargs: object,
    ) -> tuple[Path, Path, dict[str, object]]:
        del system, template_prmtop, template_inpcrd, kwargs
        directory = Path(output_directory)
        topology = directory / "structure.prmtop"
        coordinates = directory / "structure.inpcrd"
        topology.touch()
        coordinates.touch()
        return topology, coordinates, {"status": "passed"}

    monkeypatch.setattr("cellulose_forge.cli.write_glycam_amber_system", write_amber)
    monkeypatch.setattr(
        "cellulose_forge.cli.validate_openmm_system",
        lambda *args, **kwargs: {
            "status": "passed",
            "openmm_version": "8.5",
            "platform": "CPU",
            "initial_potential_energy_kj_mol": 100.0,
            "minimized_potential_energy_kj_mol": 80.0,
        },
    )
    output = tmp_path / "ibeta_openmm"

    result = main(
        [
            "build",
            "--allomorph",
            "I-beta",
            "--shape",
            "slab",
            "--face",
            "1",
            "1",
            "0",
            "--layers",
            "1",
            "--chains-per-layer",
            "2",
            "--degree-polymerization",
            "12",
            "--chain-boundary",
            "finite",
            "--vacuum",
            "30",
            "--padding-longitudinal",
            "12",
            "--force-field",
            "glycam",
            "--md-engine",
            "gromacs",
            "--md-engine",
            "openmm",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    report = json.loads((output / "build_report.json").read_text(encoding="utf-8"))
    assert report["md_engines"] == ["gromacs", "openmm"]
    assert report["openmm_validation"]["status"] == "passed"
    assert report["outputs"]["amber_topology"] == "structure.prmtop"
    assert report["outputs"]["openmm_runner"] == "run_openmm.py"
    assert report["outputs"]["gromacs_topology"] == "system.top"


def test_cli_rejects_unpadded_finite_force_field_slab(tmp_path: Path) -> None:
    with pytest.raises(SystemExit, match="2"):
        main(
            [
                "build",
                "--allomorph",
                "I-beta",
                "--shape",
                "slab",
                "--face",
                "1",
                "1",
                "0",
                "--layers",
                "1",
                "--chains-per-layer",
                "2",
                "--degree-polymerization",
                "12",
                "--chain-boundary",
                "finite",
                "--force-field",
                "glycam",
                "--output",
                str(tmp_path / "unsafe"),
            ]
        )


def test_cli_accepts_deprecated_longitudinal_padding_alias(tmp_path: Path) -> None:
    output = tmp_path / "finite_coordinates.pdb"
    with pytest.warns(FutureWarning, match="padding-z is deprecated"):
        assert (
            main(
                [
                    "build",
                    "--allomorph",
                    "I-beta",
                    "--shape",
                    "slab",
                    "--face",
                    "1",
                    "1",
                    "0",
                    "--layers",
                    "1",
                    "--chains-per-layer",
                    "2",
                    "--degree-polymerization",
                    "12",
                    "--chain-boundary",
                    "finite",
                    "--padding-z",
                    "12",
                    "--output",
                    str(output),
                ]
            )
            == 0
        )
