from collections.abc import Mapping
from pathlib import Path

import pytest

from cellulose_forge.energy_drivers.amber import (
    AMBER_COULOMB_KCAL_ANGSTROM,
    COULOMB_SCALE_TO_GROMACS,
    GROMACS_COULOMB_KCAL_ANGSTROM,
    KJ_PER_KCAL,
    parse_amber_energy,
)
from cellulose_forge.energy_drivers.amber import (
    normalize_glycam_amber_energies as _normalize_glycam_amber_energies,
)
from cellulose_forge.energy_drivers.report import EnergyReport, EnergyValue, validate_potential_sum

_SANDER_BLOCK = """
imin=1, maxcyc=1000, ncyc=200,
 NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
 Etot   =      81.0000  EKtot   =      0.0000  EPtot      =      81.0000
 BOND   =      10.0000  ANGLE   =       5.0000  DIHED      =       9.0000
 1-4 NB =       7.0000  1-4 EEL =      11.0000  VDWAALS    =       8.0000
 EELEC  =      19.0000  EHBOND  =       0.0000  RESTRAINT  =      12.0000
 EAMBER (non-restraint)  =      69.0000
"""


def normalize_glycam_amber_energies(
    raw_terms: Mapping[str, EnergyValue],
    *,
    coordinate_sha256: str,
    periodic_improper: float = 0.0,
) -> EnergyReport:
    return _normalize_glycam_amber_energies(
        raw_terms,
        coordinate_sha256=coordinate_sha256,
        periodic_improper=periodic_improper,
    )


def _kcal_terms() -> dict[str, float]:
    return {
        "NSTEP": 0,
        "Etot": 81.0,
        "EKtot": 0.0,
        "EPtot": 81.0,
        "BOND": 10.0,
        "ANGLE": 5.0,
        "DIHED": 9.0,
        "1-4 NB": 7.0,
        "1-4 EEL": 11.0,
        "VDWAALS": 8.0,
        "EELEC": 19.0,
        "EHBOND": 0.0,
        "RESTRAINT": 12.0,
        "EAMBER (non-restraint)": 69.0,
    }


def test_normalize_glycam_amber_energies_uses_canonical_complete_terms() -> None:
    report = normalize_glycam_amber_energies(_kcal_terms(), coordinate_sha256="abc")

    assert report.engine == "amber"
    assert report.extra["scope"] == "finite-chain GLYCAM"
    assert report.terms == {
        "bond": 10.0 * KJ_PER_KCAL,
        "angle": 5.0 * KJ_PER_KCAL,
        "dihedral": 9.0 * KJ_PER_KCAL,
        "improper": 0.0,
        "vdw": 15.0 * KJ_PER_KCAL,
        "electrostatics": 30.0 * KJ_PER_KCAL * COULOMB_SCALE_TO_GROMACS,
        "restraint": 12.0 * KJ_PER_KCAL,
        "dispersion_correction": 0.0,
        "potential": 81.0 * KJ_PER_KCAL + 30.0 * KJ_PER_KCAL * (COULOMB_SCALE_TO_GROMACS - 1.0),
    }
    validate_potential_sum(report, absolute_kj_mol=1e-9, relative=0.0)


def test_coulomb_constant_is_explicit_and_canonical_terms_are_reference_comparable() -> None:
    report = normalize_glycam_amber_energies(_kcal_terms(), coordinate_sha256="abc")
    assert report.corrections is not None
    correction = report.corrections[0]
    raw_elec_kcal = 30.0
    raw_elec_kj = raw_elec_kcal * KJ_PER_KCAL

    assert correction["kind"] == "scale"
    assert correction["factor"] == COULOMB_SCALE_TO_GROMACS
    assert correction["engine_constant_kcal_angstrom_per_mol_e2"] == AMBER_COULOMB_KCAL_ANGSTROM
    assert correction["reference_constant_kcal_angstrom_per_mol_e2"] == (
        GROMACS_COULOMB_KCAL_ANGSTROM
    )
    assert correction["raw_kcal_mol"] == raw_elec_kcal
    assert correction["raw_kj_mol"] == raw_elec_kj
    assert correction["name"] == "amber_coulomb_constant"
    assert correction["correction_kj_mol"] == pytest.approx(
        raw_elec_kj * (COULOMB_SCALE_TO_GROMACS - 1.0)
    )
    assert correction["comparable_kj_mol"] == pytest.approx(raw_elec_kj * COULOMB_SCALE_TO_GROMACS)
    assert correction["comparable_kj_mol"] != correction["raw_kj_mol"]
    assert report.terms["electrostatics"] == pytest.approx(raw_elec_kj * COULOMB_SCALE_TO_GROMACS)
    gromacs_convention_back_to_amber = (
        correction["comparable_kj_mol"]
        * AMBER_COULOMB_KCAL_ANGSTROM
        / GROMACS_COULOMB_KCAL_ANGSTROM
    )
    assert gromacs_convention_back_to_amber == pytest.approx(raw_elec_kj)


def test_periodic_improper_is_split_from_amber_dihed() -> None:
    report = normalize_glycam_amber_energies(
        _kcal_terms(), coordinate_sha256="abc", periodic_improper=2.5
    )

    assert report.terms["dihedral"] == 6.5 * KJ_PER_KCAL
    assert report.terms["improper"] == 2.5 * KJ_PER_KCAL
    validate_potential_sum(report, absolute_kj_mol=1e-9, relative=0.0)


def test_parse_amber_energy_ignores_namelist_and_later_blocks() -> None:
    text = (
        _SANDER_BLOCK
        + """
 NSTEP =      100   TIME(PS) =       0.100  TEMP(K) =   300.00  PRESS =     0.0
 Etot   =     999.0000  EKtot   =    100.0000  EPtot      =     899.0000
 BOND   =     999.0000  ANGLE   =     999.0000  DIHED      =     999.0000
 1-4 NB =     999.0000  1-4 EEL =     999.0000  VDWAALS    =     999.0000
 EELEC  =     999.0000  EHBOND  =       0.0000  RESTRAINT  =       0.0000
"""
    )

    parsed = parse_amber_energy(text)

    assert parsed["BOND"] == 10.0
    assert parsed["EPtot"] == 81.0
    assert parsed["1-4 NB"] == 7.0
    assert "imin" not in parsed
    assert parsed["NSTEP"] == 0


def test_parse_amber_energy_path(tmp_path: Path) -> None:
    path = tmp_path / "amber_sp.out"
    path.write_text(_SANDER_BLOCK, encoding="utf-8")

    parsed = parse_amber_energy(path)
    report = normalize_glycam_amber_energies(parsed, coordinate_sha256="abc")

    assert parsed["EELEC"] == 19.0
    assert report.terms["bond"] == 10.0 * KJ_PER_KCAL


def test_parse_amber_minimization_single_point_without_eptot() -> None:
    text = """
   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
      0       81.0000       1.0000       2.0000          C1       1
 BOND   =      10.0000  ANGLE   =       5.0000  DIHED      =       9.0000
 1-4 NB =       7.0000  1-4 EEL =      11.0000  VDWAALS    =       8.0000
 EELEC  =      19.0000  EHBOND  =       0.0000  RESTRAINT  =      12.0000
"""

    parsed = parse_amber_energy(text)

    assert parsed["EPtot"] == 81.0
    assert parsed["BOND"] == 10.0


def test_parse_amber_minimization_accepts_correct_five_digit_rounding() -> None:
    text = """
   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
      0       1.2345E+04    1.0000       2.0000          C1       1
 BOND   =   12300.0000  ANGLE   =      20.0000  DIHED      =      10.0000
 1-4 NB =       5.0000  1-4 EEL =       4.0000  VDWAALS    =       3.0000
 EELEC  =       2.0000  EHBOND  =       0.0000  RESTRAINT  =       1.4377
"""

    parsed = parse_amber_energy(text)

    assert parsed["MINIMIZATION ENERGY"] == 12345.0
    assert parsed["EPtot"] == pytest.approx(12345.4377)


def test_parse_amber_minimization_rejects_difference_beyond_printed_precision() -> None:
    text = """
   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
      0       1.2344E+04    1.0000       2.0000          C1       1
 BOND   =   12300.0000  ANGLE   =      20.0000  DIHED      =      10.0000
 1-4 NB =       5.0000  1-4 EEL =       4.0000  VDWAALS    =       3.0000
 EELEC  =       2.0000  EHBOND  =       0.0000  RESTRAINT  =       1.4377
"""

    with pytest.raises(ValueError, match="differs from its printed component sum"):
        parse_amber_energy(text)


def test_parse_amber_minimization_rejects_multiple_energy_rows() -> None:
    text = """
   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
      0       81.0000       1.0000       2.0000          C1       1
      1       80.0000       1.0000       2.0000          C1       1
 BOND   =      10.0000  ANGLE   =       5.0000  DIHED      =       9.0000
 1-4 NB =       7.0000  1-4 EEL =      11.0000  VDWAALS    =       8.0000
 EELEC  =      19.0000  EHBOND  =       0.0000  RESTRAINT  =      12.0000
"""

    with pytest.raises(ValueError, match="must contain one ENERGY row, found 2"):
        parse_amber_energy(text)


def test_periodic_improper_classification_is_required() -> None:
    with pytest.raises(TypeError, match="periodic_improper"):
        _normalize_glycam_amber_energies(_kcal_terms(), coordinate_sha256="abc")  # type: ignore[call-arg]


def test_one_four_vdw_alias_is_accepted() -> None:
    raw = _kcal_terms()
    raw["1-4 VDW"] = raw.pop("1-4 NB")

    report = normalize_glycam_amber_energies(raw, coordinate_sha256="abc")

    assert report.terms["vdw"] == 15.0 * KJ_PER_KCAL


def test_conflicting_one_four_labels_fail() -> None:
    raw = _kcal_terms()
    raw["1-4 VDW"] = 1.0

    with pytest.raises(ValueError, match="conflicting energy labels"):
        normalize_glycam_amber_energies(raw, coordinate_sha256="abc")


def test_charmm_terms_fail_closed_for_glycam_adapter() -> None:
    raw = _kcal_terms()
    raw["UB"] = 1.5

    with pytest.raises(ValueError, match=r"finite-chain GLYCAM.*UB"):
        normalize_glycam_amber_energies(raw, coordinate_sha256="abc")


def test_ehbond_nonzero_is_not_glycam() -> None:
    raw = _kcal_terms()
    raw["EHBOND"] = 0.25

    with pytest.raises(ValueError, match="EHBOND"):
        normalize_glycam_amber_energies(raw, coordinate_sha256="abc")


def test_unknown_zero_diagnostic_is_allowed() -> None:
    raw = _kcal_terms()
    raw["CMAP"] = 0.0
    raw["Density"] = 1.2

    report = normalize_glycam_amber_energies(raw, coordinate_sha256="abc")

    assert report.terms["dispersion_correction"] == 0.0


def test_edisper_maps_to_dispersion_correction() -> None:
    raw = _kcal_terms()
    raw["EDISPER"] = -1.0
    raw["EPtot"] = 80.0

    report = normalize_glycam_amber_energies(raw, coordinate_sha256="abc")

    assert report.terms["dispersion_correction"] == -1.0 * KJ_PER_KCAL
    validate_potential_sum(report, absolute_kj_mol=1e-9, relative=0.0)


def test_normalize_preserves_raw_kcal_diagnostics() -> None:
    raw = _kcal_terms()
    report = normalize_glycam_amber_energies(raw, coordinate_sha256="abc")
    raw["BOND"] = 999.0

    assert report.extra["raw_terms"]["BOND"] == 10.0
    assert report.extra["raw_units"] == "kcal/mol"


@pytest.mark.parametrize("missing", ["BOND", "ANGLE", "DIHED", "VDWAALS", "EPtot", "1-4 EEL"])
def test_normalize_rejects_missing_glycam_term(missing: str) -> None:
    raw = _kcal_terms()
    del raw[missing]

    with pytest.raises(ValueError, match="missing"):
        normalize_glycam_amber_energies(raw, coordinate_sha256="abc")


def test_normalize_rejects_missing_one_four_vdw() -> None:
    raw = _kcal_terms()
    del raw["1-4 NB"]

    with pytest.raises(ValueError, match="1-4 van der Waals"):
        normalize_glycam_amber_energies(raw, coordinate_sha256="abc")


@pytest.mark.parametrize("value", [float("nan"), float("inf"), True, "1.0"])
def test_normalize_rejects_invalid_values(value: object) -> None:
    raw: dict[str, object] = dict(_kcal_terms())
    raw["BOND"] = value

    with pytest.raises((TypeError, ValueError), match="BOND"):
        normalize_glycam_amber_energies(raw, coordinate_sha256="abc")  # type: ignore[arg-type]


def test_parse_amber_energy_rejects_missing_or_nonfinite_block() -> None:
    with pytest.raises(ValueError, match="no EPtot"):
        parse_amber_energy("imin=1, ntb=1\n")
    with pytest.raises(ValueError, match="not finite"):
        parse_amber_energy(" EPtot      =       1e999\n")
