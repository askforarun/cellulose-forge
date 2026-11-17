import struct
from pathlib import Path
from types import SimpleNamespace

import pytest

from cellulose_forge.energy_drivers.coordinates import write_coordinate_reference
from cellulose_forge.energy_drivers.namd import (
    KJ_PER_KCAL,
    parse_namd_energy_log,
    run_namd_single_point,
)
from cellulose_forge.energy_drivers.namd import (
    normalize_namd_energies as _normalize_namd_energies,
)
from cellulose_forge.energy_drivers.report import EnergyReport, EnergyValue, validate_potential_sum

_STANDARD_TITLES = (
    "TS",
    "BOND",
    "ANGLE",
    "DIHED",
    "IMPRP",
    "ELECT",
    "VDW",
    "BOUNDARY",
    "MISC",
    "KINETIC",
    "TOTAL",
    "TEMP",
    "POTENTIAL",
    "VOLUME",
)


def normalize_namd_energies(
    raw_terms: dict[str, EnergyValue],
    *,
    coordinate_sha256: str,
    classified_terms: dict[str, str] | None = None,
) -> EnergyReport:
    return _normalize_namd_energies(
        raw_terms,
        coordinate_sha256=coordinate_sha256,
        merge_crossterms=False,
        lj_correction=False,
        classified_terms=classified_terms,
    )


def _kcal_terms() -> dict[str, float]:
    return {
        "TS": 0,
        "BOND": 1.0,
        "ANGLE": 5.0,
        "DIHED": 3.0,
        "IMPRP": 4.0,
        "ELECT": 26.0,
        "VDW": 19.0,
        "BOUNDARY": 0.0,
        "MISC": 0.0,
        "KINETIC": 12.0,
        "TOTAL": 70.0,
        "TEMP": 300.0,
        "POTENTIAL": 58.0,
        "VOLUME": 1000.0,
    }


def _energy_log(titles: tuple[str, ...], values: dict[str, float], extra_rows: str = "") -> str:
    row = " ".join(str(values[name]) for name in titles)
    return f"ETITLE: {' '.join(titles)}\nENERGY: {row}\n{extra_rows}"


def test_normalize_namd_energies_uses_charmm_semantics_and_complete_terms() -> None:
    report = normalize_namd_energies(_kcal_terms(), coordinate_sha256="abc")

    assert report.engine == "namd"
    assert report.terms == {
        "bond": 1.0 * KJ_PER_KCAL,
        "angle": 5.0 * KJ_PER_KCAL,
        "dihedral": 3.0 * KJ_PER_KCAL,
        "improper": 4.0 * KJ_PER_KCAL,
        "vdw": 19.0 * KJ_PER_KCAL,
        "electrostatics": 26.0 * KJ_PER_KCAL,
        "restraint": 0.0,
        "dispersion_correction": 0.0,
        "potential": 58.0 * KJ_PER_KCAL,
    }
    assert report.extra["raw_units"] == "kcal/mol"
    validate_potential_sum(report, absolute_kj_mol=1e-9, relative=0.0)


def test_namd_angle_is_not_split_from_urey_bradley() -> None:
    raw = _kcal_terms()
    report = normalize_namd_energies(raw, coordinate_sha256="abc")

    assert report.terms["angle"] == raw["ANGLE"] * KJ_PER_KCAL
    assert report.terms["bond"] == raw["BOND"] * KJ_PER_KCAL
    assert report.terms["dihedral"] == raw["DIHED"] * KJ_PER_KCAL
    assert report.terms["improper"] == raw["IMPRP"] * KJ_PER_KCAL


@pytest.mark.parametrize("term", ["BOUNDARY", "MISC"])
def test_unclassified_nonzero_boundary_or_misc_fails(term: str) -> None:
    raw = _kcal_terms()
    raw[term] = 1.5
    raw["POTENTIAL"] += 1.5

    with pytest.raises(ValueError, match=rf"unclassified nonzero BOUNDARY/MISC.*{term}"):
        normalize_namd_energies(raw, coordinate_sha256="abc")


def test_explicitly_classified_boundary_and_misc_restraints_pass() -> None:
    raw = _kcal_terms()
    raw["MISC"] = 2.0
    raw["BOUNDARY"] = 1.5
    raw["POTENTIAL"] = 61.5

    report = normalize_namd_energies(
        raw,
        coordinate_sha256="abc",
        classified_terms={"BOUNDARY": "restraint", "MISC": "restraint"},
    )

    assert report.terms["restraint"] == 3.5 * KJ_PER_KCAL
    assert report.extra["classified_terms"] == {"BOUNDARY": "restraint", "MISC": "restraint"}
    validate_potential_sum(report, absolute_kj_mol=1e-9, relative=0.0)


def test_zero_boundary_and_misc_do_not_require_classification() -> None:
    report = normalize_namd_energies(_kcal_terms(), coordinate_sha256="abc")

    assert report.terms["restraint"] == 0.0
    assert "classified_terms" not in report.extra
    validate_potential_sum(report, absolute_kj_mol=1e-9, relative=0.0)


@pytest.mark.parametrize(
    ("settings", "message"),
    [
        ({"merge_crossterms": True, "lj_correction": False}, "mergeCrossterms no"),
        ({"merge_crossterms": False, "lj_correction": True}, "LJcorrection off"),
    ],
)
def test_namd_single_point_settings_fail_closed(settings: dict[str, bool], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        _normalize_namd_energies(
            _kcal_terms(), coordinate_sha256="abc", classified_terms=None, **settings
        )


def test_parse_namd_energy_log_uses_etitle_names_not_fixed_columns() -> None:
    values = _kcal_terms()
    shuffled = (
        "POTENTIAL",
        "VDW",
        "TS",
        "ELECT",
        "IMPRP",
        "BOND",
        "ANGLE",
        "DIHED",
        "MISC",
        "BOUNDARY",
        "KINETIC",
        "TOTAL",
        "TEMP",
        "VOLUME",
    )

    parsed = parse_namd_energy_log(_energy_log(shuffled, values))
    report = normalize_namd_energies(parsed, coordinate_sha256="abc")

    assert parsed["BOND"] == 1.0
    assert parsed["ANGLE"] == 5.0
    assert report.terms["bond"] == 1.0 * KJ_PER_KCAL
    assert report.terms["angle"] == 5.0 * KJ_PER_KCAL
    assert report.terms["potential"] == 58.0 * KJ_PER_KCAL


def test_parse_namd_energy_log_refuses_more_than_one_energy_row(tmp_path: Path) -> None:
    path = tmp_path / "namd.log"
    extra = "ENERGY: 1 999 999 999 999 999 999 0 0 0 999 300 999 0\n"
    path.write_text(_energy_log(_STANDARD_TITLES, _kcal_terms(), extra), encoding="utf-8")

    with pytest.raises(ValueError, match="must contain one ENERGY row"):
        parse_namd_energy_log(path)


def test_run_namd_single_point_checks_shared_binary_coordinates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference = tmp_path / "coordinate_reference.json"
    checksum = write_coordinate_reference(
        reference, [[1.0, 2.0, 3.0]], [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    )
    binary = tmp_path / "single_point.coor"
    binary.write_bytes(struct.pack("=i3d", 1, 1.0, 2.0, 3.0))
    configuration = tmp_path / "single_point_namd.conf"
    configuration.write_text(
        "binCoordinates single_point.coor\ncellBasisVector1 10 0 0\n"
        "cellBasisVector2 0 10 0\ncellBasisVector3 0 0 10\n"
        "mergeCrossterms no\nLJcorrection off\nrigidBonds none\nrun 0\n",
        encoding="utf-8",
    )
    log = _energy_log(_STANDARD_TITLES, _kcal_terms())
    monkeypatch.setattr(
        "cellulose_forge.energy_drivers.namd.subprocess.run",
        lambda *args, **kwargs: SimpleNamespace(returncode=0, stdout=log, stderr=""),
    )

    report = run_namd_single_point(
        tmp_path,
        configuration_file=configuration,
        coordinate_reference_file=reference,
        binary_coordinate_file=binary,
        coordinate_sha256=checksum,
    )

    assert report.terms["potential"] == 58.0 * KJ_PER_KCAL
    assert (tmp_path / "single_point_namd.log").read_text() == log


def test_run_namd_single_point_rejects_unmatched_binary_coordinates(tmp_path: Path) -> None:
    reference = tmp_path / "coordinate_reference.json"
    checksum = write_coordinate_reference(
        reference, [[1.0, 2.0, 3.0]], [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    )
    binary = tmp_path / "single_point.coor"
    binary.write_bytes(struct.pack("=i3d", 1, 1.0, 2.0, 3.1))
    configuration = tmp_path / "single_point_namd.conf"
    configuration.write_text(
        "binCoordinates single_point.coor\ncellBasisVector1 10 0 0\n"
        "cellBasisVector2 0 10 0\ncellBasisVector3 0 0 10\n"
        "mergeCrossterms no\nLJcorrection off\nrigidBonds none\nrun 0\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="coordinate reference"):
        run_namd_single_point(
            tmp_path,
            configuration_file=configuration,
            coordinate_reference_file=reference,
            binary_coordinate_file=binary,
            coordinate_sha256=checksum,
        )


@pytest.mark.parametrize(
    ("replacement", "message"),
    [
        ("binCoordinates other.coor", "binCoordinates does not name"),
        ("cellBasisVector1 11 0 0", "box differs"),
    ],
)
def test_run_namd_single_point_checks_configured_native_files(
    tmp_path: Path, replacement: str, message: str
) -> None:
    reference = tmp_path / "coordinate_reference.json"
    checksum = write_coordinate_reference(
        reference, [[1.0, 2.0, 3.0]], [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    )
    binary = tmp_path / "single_point.coor"
    binary.write_bytes(struct.pack("=i3d", 1, 1.0, 2.0, 3.0))
    configuration = tmp_path / "single_point_namd.conf"
    lines = [
        "binCoordinates single_point.coor",
        "cellBasisVector1 10 0 0",
        "cellBasisVector2 0 10 0",
        "cellBasisVector3 0 0 10",
        "mergeCrossterms no",
        "LJcorrection off",
        "rigidBonds none",
        "run 0",
    ]
    if replacement.startswith("binCoordinates"):
        lines[0] = replacement
    else:
        lines[1] = replacement
    configuration.write_text("\n".join(lines) + "\n")

    with pytest.raises(ValueError, match=message):
        run_namd_single_point(
            tmp_path,
            configuration_file=configuration,
            coordinate_reference_file=reference,
            binary_coordinate_file=binary,
            coordinate_sha256=checksum,
        )


@pytest.mark.parametrize(
    "directive", ["cwd elsewhere", "extendedSystem state.xsc", "source override.tcl"]
)
def test_run_namd_single_point_refuses_configuration_indirection(
    tmp_path: Path, directive: str
) -> None:
    reference = tmp_path / "coordinate_reference.json"
    checksum = write_coordinate_reference(
        reference, [[1.0, 2.0, 3.0]], [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    )
    binary = tmp_path / "single_point.coor"
    binary.write_bytes(struct.pack("=i3d", 1, 1.0, 2.0, 3.0))
    configuration = tmp_path / "single_point_namd.conf"
    configuration.write_text(
        "binCoordinates single_point.coor\ncellBasisVector1 10 0 0\n"
        "cellBasisVector2 0 10 0\ncellBasisVector3 0 0 10\n"
        "mergeCrossterms no\nLJcorrection off\nrigidBonds none\n"
        f"{directive}\nrun 0\n"
    )

    with pytest.raises(ValueError, match="outside the verified subset"):
        run_namd_single_point(
            tmp_path,
            configuration_file=configuration,
            coordinate_reference_file=reference,
            binary_coordinate_file=binary,
            coordinate_sha256=checksum,
        )


@pytest.mark.parametrize(
    "injection",
    [
        "cutoff 12.0 ; cd /elsewhere",
        "temperature [cd /elsewhere; expr 0]",
        "exclude scaled1-4; source evil.tcl",
        "cutoff $unsafe",
        "cutoff {12.0}",
        "cutoff 12.0 \\",
    ],
)
def test_run_namd_single_point_refuses_tcl_command_injection(
    tmp_path: Path, injection: str
) -> None:
    reference = tmp_path / "coordinate_reference.json"
    checksum = write_coordinate_reference(
        reference, [[1.0, 2.0, 3.0]], [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    )
    binary = tmp_path / "single_point.coor"
    binary.write_bytes(struct.pack("=i3d", 1, 1.0, 2.0, 3.0))
    configuration = tmp_path / "single_point_namd.conf"
    configuration.write_text(
        "binCoordinates single_point.coor\ncellBasisVector1 10 0 0\n"
        "cellBasisVector2 0 10 0\ncellBasisVector3 0 0 10\n"
        "mergeCrossterms no\nLJcorrection off\nrigidBonds none\n"
        f"{injection}\nrun 0\n"
    )

    with pytest.raises(ValueError, match="Tcl metacharacters"):
        run_namd_single_point(
            tmp_path,
            configuration_file=configuration,
            coordinate_reference_file=reference,
            binary_coordinate_file=binary,
            coordinate_sha256=checksum,
        )


def test_run_namd_single_point_requires_configuration_in_engine_directory(tmp_path: Path) -> None:
    engine = tmp_path / "namd"
    engine.mkdir()
    reference = engine / "coordinate_reference.json"
    checksum = write_coordinate_reference(
        reference, [[1.0, 2.0, 3.0]], [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    )
    binary = engine / "single_point.coor"
    binary.write_bytes(struct.pack("=i3d", 1, 1.0, 2.0, 3.0))
    configuration = tmp_path / "single_point_namd.conf"
    configuration.write_text("binCoordinates single_point.coor\n")

    with pytest.raises(ValueError, match="inside its engine directory"):
        run_namd_single_point(
            engine,
            configuration_file=configuration,
            coordinate_reference_file=reference,
            binary_coordinate_file=binary,
            coordinate_sha256=checksum,
        )


def test_run_namd_single_point_refuses_repeated_checked_directive(tmp_path: Path) -> None:
    reference = tmp_path / "coordinate_reference.json"
    checksum = write_coordinate_reference(
        reference, [[1.0, 2.0, 3.0]], [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    )
    binary = tmp_path / "single_point.coor"
    binary.write_bytes(struct.pack("=i3d", 1, 1.0, 2.0, 3.0))
    configuration = tmp_path / "single_point_namd.conf"
    configuration.write_text(
        "binCoordinates single_point.coor\nbinCoordinates other.coor\n"
        "cellBasisVector1 10 0 0\ncellBasisVector2 0 10 0\n"
        "cellBasisVector3 0 0 10\nmergeCrossterms no\nLJcorrection off\n"
        "rigidBonds none\nrun 0\n"
    )

    with pytest.raises(ValueError, match="exactly once"):
        run_namd_single_point(
            tmp_path,
            configuration_file=configuration,
            coordinate_reference_file=reference,
            binary_coordinate_file=binary,
            coordinate_sha256=checksum,
        )


def test_namd_attestations_are_required() -> None:
    with pytest.raises(TypeError, match="merge_crossterms"):
        _normalize_namd_energies(_kcal_terms(), coordinate_sha256="abc")  # type: ignore[call-arg]


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        ({"merge_crossterms": 0, "lj_correction": False}, "merge_crossterms"),
        ({"merge_crossterms": False, "lj_correction": 0}, "lj_correction"),
    ],
)
def test_namd_attestations_require_actual_booleans(
    arguments: dict[str, object], message: str
) -> None:
    with pytest.raises(TypeError, match=message):
        _normalize_namd_energies(  # type: ignore[arg-type]
            _kcal_terms(),
            coordinate_sha256="abc",
            **arguments,
        )


def test_normalize_preserves_raw_kcal_diagnostics() -> None:
    raw = _kcal_terms()
    report = normalize_namd_energies(raw, coordinate_sha256="abc")
    raw["BOND"] = 999.0

    assert report.extra["raw_terms"]["BOND"] == 1.0
    assert report.extra["raw_terms"]["KINETIC"] == 12.0


def test_crossterm_fails_closed_when_nonzero() -> None:
    raw = _kcal_terms()
    raw["CROSSTERM"] = 3.0

    with pytest.raises(ValueError, match=r"unsupported nonzero.*CROSSTERM"):
        normalize_namd_energies(raw, coordinate_sha256="abc")


def test_unknown_zero_and_nonzero_kinetic_diagnostics_are_allowed() -> None:
    raw = _kcal_terms()
    raw["CROSSTERM"] = 0.0
    raw["KINETIC"] = 250.0

    report = normalize_namd_energies(raw, coordinate_sha256="abc")

    assert report.terms["potential"] == 58.0 * KJ_PER_KCAL


@pytest.mark.parametrize(
    "missing", ["BOND", "ANGLE", "DIHED", "IMPRP", "ELECT", "VDW", "POTENTIAL"]
)
def test_normalize_rejects_missing_charmm_term(missing: str) -> None:
    raw = _kcal_terms()
    del raw[missing]

    with pytest.raises(ValueError, match="missing required CHARMM"):
        normalize_namd_energies(raw, coordinate_sha256="abc")


@pytest.mark.parametrize("value", [float("nan"), float("inf"), True, "1.0"])
def test_normalize_rejects_invalid_values(value: object) -> None:
    raw: dict[str, object] = dict(_kcal_terms())
    raw["BOND"] = value

    with pytest.raises((TypeError, ValueError), match="BOND"):
        normalize_namd_energies(raw, coordinate_sha256="abc")  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("log", "message"),
    [
        ("ENERGY: 0 1 2\n", "ENERGY before ETITLE"),
        ("ETITLE: TS BOND ANGLE\nENERGY: 0 1\n", "field count differs"),
        ("ETITLE: TS BOND BOND\nENERGY: 0 1 2\n", "duplicate field names"),
        ("ETITLE:\nENERGY: 0\n", "no fields"),
        ("ETITLE: TS BOND\nENERGY: 0 not-a-number\n", "not numeric"),
        ("ETITLE: TS BOND\nENERGY: 0 nan\n", "not finite"),
        ("# no energy\n", "no ENERGY row"),
    ],
)
def test_parse_namd_energy_log_rejects_malformed_logs(log: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        parse_namd_energy_log(log)
