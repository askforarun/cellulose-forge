"""Tests for finite-chain GLYCAM reference generation."""

from pathlib import Path

import pytest

from cellulose_forge.glycam import (
    GlycamChainFiles,
    GlycamChainSpec,
    GlycamPeriodicChainSpec,
    _gromacs_section_rows,
    _replace_gromacs_molecule_name,
    _tleap_output_failed,
    write_tleap_chain_input,
    write_tleap_periodic_chain_input,
)


@pytest.mark.parametrize(
    "output",
    [
        "Errors = 1",
        "Errors=2",
        "Errors   =   17",
        "FATAL: missing atom type",
        "Parameter file was not saved",
    ],
)
def test_tleap_failure_detection_catches_any_reported_error(output: str) -> None:
    assert _tleap_output_failed(output)


@pytest.mark.parametrize("output", ["Errors = 0", "Checking Unit.", ""])
def test_tleap_failure_detection_accepts_success_output(output: str) -> None:
    assert not _tleap_output_failed(output)


def test_dp12_glycam_sequence_has_correct_terminals_and_internal_residues() -> None:
    spec = GlycamChainSpec(12)

    assert spec.tleap_residues == ("ROH", *("4GB",) * 11, "0GB")
    assert spec.chemical_residues == ("0GB", *("4GB",) * 11, "ROH")


def test_dp6_chitin_glycam_sequence_uses_glcnac_residues() -> None:
    spec = GlycamChainSpec(6, residue_code="Y")

    assert spec.tleap_residues == ("ROH", *("4YB",) * 5, "0YB")
    assert spec.chemical_residues == ("0YB", *("4YB",) * 5, "ROH")


def test_periodic_dp6_chitin_uses_neutral_internal_glcnac() -> None:
    spec = GlycamPeriodicChainSpec(6, residue_code="Y")

    assert spec.tleap_residues == ("4YB",) * 6
    assert spec.internal_atom_count == 27


@pytest.mark.parametrize("degree_polymerization", [0, -1, 1.5, True])
def test_glycam_sequence_rejects_invalid_dp(degree_polymerization: object) -> None:
    with pytest.raises(ValueError, match="positive integer"):
        GlycamChainSpec(degree_polymerization)  # type: ignore[arg-type]


def test_periodic_dp12_uses_only_internal_residues_and_closure_bond(tmp_path: Path) -> None:
    spec = GlycamPeriodicChainSpec(12)

    files = write_tleap_periodic_chain_input(spec, tmp_path)
    text = files.leap.read_text(encoding="utf-8")

    assert spec.tleap_residues == ("4GB",) * 12
    assert "chain = sequence { " + " ".join(("4GB",) * 12) + " }" in text
    assert "bond chain.12.O4 chain.1.C1" in text
    assert "ROH" not in text
    assert "0GB" not in text


def test_periodic_dp12_topology_has_complete_closure_interactions() -> None:
    topology = Path(
        "examples/glycam/dp12_periodic_reference/cellulose_periodic_chain.top"
    ).read_text(encoding="utf-8")
    atom_residues = {
        int(fields[0]): int(fields[2])
        for row in _gromacs_section_rows(topology, "atoms")
        if (fields := row.split())
    }

    def boundary_count(section: str, atom_columns: int) -> int:
        count = 0
        for row in _gromacs_section_rows(topology, section):
            atom_indices = [int(value) for value in row.split()[:atom_columns]]
            interaction_residues = {atom_residues[index] for index in atom_indices}
            if {1, 12}.issubset(interaction_residues):
                count += 1
        return count

    assert boundary_count("bonds", 2) == 1
    assert boundary_count("angles", 3) == 4
    assert boundary_count("dihedrals", 4) == 16
    assert boundary_count("pairs", 2) == 10


def test_tleap_input_uses_pinned_glycam_release(tmp_path: Path) -> None:
    spec = GlycamChainSpec(4)

    files = write_tleap_chain_input(spec, tmp_path)

    assert files.leap.read_text(encoding="utf-8") == "\n".join(
        (
            "source leaprc.GLYCAM_06j-1",
            "chain = sequence { ROH 4GB 4GB 4GB 0GB }",
            "check chain",
            "savepdb chain cellulose_chain.pdb",
            "saveamberparm chain cellulose_chain.prmtop cellulose_chain.inpcrd",
            "quit",
            "",
        )
    )


def test_glycam_chain_files_include_gromacs_outputs(tmp_path: Path) -> None:
    files = GlycamChainFiles(tmp_path, "dp12")

    assert files.gromacs_topology == tmp_path / "dp12.top"
    assert files.gromacs_coordinates == tmp_path / "dp12.gro"


def test_gromacs_section_rows_ignores_comments_and_other_sections() -> None:
    topology = """
[ defaults ]
; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ
1 2 yes 1 1

[ molecules ]
; molecule count
CELLULOSE_CHAIN 8
"""

    assert _gromacs_section_rows(topology, "defaults") == ["1 2 yes 1 1"]
    assert _gromacs_section_rows(topology, "molecules") == ["CELLULOSE_CHAIN 8"]


def test_replace_gromacs_molecule_name_updates_definition_and_count_entry() -> None:
    topology = """
[ moleculetype ]
; Name nrexcl
system1 3

[ molecules ]
system1 1
"""

    renamed = _replace_gromacs_molecule_name(topology, molecule_name="CELLULOSE_CHAIN")

    assert _gromacs_section_rows(renamed, "moleculetype") == ["CELLULOSE_CHAIN  3"]
    assert _gromacs_section_rows(renamed, "molecules") == ["CELLULOSE_CHAIN  1"]
