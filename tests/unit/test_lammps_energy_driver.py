from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest

from cellulose_forge.energy_drivers import lammps
from cellulose_forge.energy_drivers.coordinates import write_coordinate_reference
from cellulose_forge.energy_drivers.lammps import (
    CHARMM_COULOMB_RATIO,
    KJ_PER_KCAL,
    CharmmfswConventions,
    charmmfsw_one_four_offset_kj_mol,
    charmmfsw_switch_shift_kj_mol,
    lammps_energy_report,
    normalize_lammps_energies,
    parse_lammps_input,
    parse_lammps_thermo,
    run_lammps_single_point,
)
from cellulose_forge.energy_drivers.report import validate_potential_sum

_VALIDATOR = Path(__file__).resolve().parents[2] / "validation" / "validate_charmm36_engines.py"

_CHARMM_INPUT = """units real
atom_style full
bond_style harmonic
angle_style charmm
dihedral_style charmmfsw
improper_style harmonic
pair_style lj/charmmfsw/coul/long 10.0 12.0
pair_modify mix arithmetic tail no
special_bonds charmm
read_data system.data
kspace_style pppm 1.0e-5
thermo 100
thermo_style custom step temp pe ebond eangle edihed evdwl ecoul elong etail eimp
run 0
"""

_PLAIN_INPUT = """units real
atom_style full
bond_style harmonic
angle_style harmonic
dihedral_style charmm
pair_style lj/cut/coul/long 10.0
pair_modify mix arithmetic tail yes
special_bonds lj/coul 0.0 0.0 1.0
read_data system.data
kspace_style pppm 1.0e-5
group cellulose_restraints id 1 2 3 &
  4 5
fix cellulose_position_restraints cellulose_restraints spring/self 2.39 xyz
fix_modify cellulose_position_restraints energy yes
thermo_style custom step temp pe ebond eangle edihed evdwl ecoul elong etail \
f_cellulose_position_restraints
run 0
"""
# The trailing backslash above is a Python line continuation, so the
# thermo_style command is a single LAMMPS line.

_CHARMM_LABELS = (
    "Step Temp PotEng E_bond E_angle E_dihed E_vdwl E_coul E_long E_tail E_impro"
).split()
_PLAIN_LABELS = (
    "Step Temp PotEng E_bond E_angle E_dihed E_vdwl E_coul E_long E_tail "
    "f_cellulose_position_restraints"
).split()

# kcal/mol, internally consistent: PotEng is the sum of the components (evdwl
# already includes etail, so etail is not added again).
_CHARMM_ROW = {
    "Step": 0.0,
    "Temp": 0.0,
    "E_bond": 100.0,
    "E_angle": 50.0,
    "E_dihed": 30.0,
    "E_vdwl": -20.0,
    "E_coul": 400.0,
    "E_long": -300.0,
    "E_tail": 0.0,
    "E_impro": 2.0,
}
_CHARMM_ROW["PotEng"] = sum(
    value for key, value in _CHARMM_ROW.items() if key not in {"Step", "Temp", "E_tail"}
)
_PLAIN_ROW = {
    "Step": 0.0,
    "Temp": 0.0,
    "E_bond": 100.0,
    "E_angle": 50.0,
    "E_dihed": 30.0,
    "E_vdwl": -20.0,
    "E_coul": 400.0,
    "E_long": -300.0,
    "E_tail": -5.0,
    "f_cellulose_position_restraints": 1.5,
}
_PLAIN_ROW["PotEng"] = sum(
    value for key, value in _PLAIN_ROW.items() if key not in {"Step", "Temp", "E_tail"}
)


def _log(labels: list[str], row: dict[str, float]) -> str:
    header = "   " + " ".join(f"{label:>14}" for label in labels)
    values = " ".join(f"{row[label]:>14.10g}" for label in labels)
    return (
        "LAMMPS (22 Jul 2025)\n"
        "thermo_style custom step temp pe\n"
        "Per MPI rank memory allocation (min/avg/max) = 12 | 12 | 12 Mbytes\n"
        f"{header}\n{values}\n"
        "Loop time of 0.01 on 1 procs for 0 steps with 7 atoms\n"
    )


def _data() -> str:
    """A small triclinic CHARMM-style LAMMPS data file with a periodic contact."""
    atoms = [
        (1, 1, 1, 0.1, 1.0, 1.0, 1.0),
        (2, 1, 2, -0.2, 2.5, 1.0, 1.0),
        (3, 1, 1, 0.1, 4.0, 1.0, 1.0),
        (4, 1, 2, -0.2, 5.5, 1.0, 1.0),
        (5, 1, 1, 0.1, 7.0, 1.0, 1.0),
        (6, 1, 2, 0.1, 8.5, 1.0, 1.0),
        # 9 angstrom from atom 1 through the z image; inside the 10-12 angstrom
        # switching region of atom 6.
        (7, 2, 1, 0.0, 1.0, 1.0, 22.0),
    ]
    atom_rows = "\n".join(
        f"{index} {molecule} {kind} {charge} {x} {y} {z}"
        for index, molecule, kind, charge, x, y, z in atoms
    )
    return (
        "LAMMPS data file\n\n"
        "7 atoms\n2 atom types\n5 bonds\n3 dihedrals\n\n"
        "0.0 30.0 xlo xhi\n0.0 30.0 ylo yhi\n0.0 30.0 zlo zhi\n1.5 0.0 0.0 xy xz yz\n\n"
        "Masses\n\n1 12.011\n2 15.999\n\n"
        f"Atoms # full\n\n{atom_rows}\n\n"
        "Bonds\n\n1 1 1 2\n2 1 2 3\n3 1 3 4\n4 1 4 5\n5 1 5 6\n\n"
        "Dihedral Coeffs\n\n1 0.2 1 180 1.0\n2 0.3 3 0 0.0\n\n"
        "Dihedrals\n\n1 1 1 2 3 4\n2 1 2 3 4 5\n3 2 3 4 5 6\n\n"
        "PairIJ Coeffs\n\n1 1 0.1 3.5 0.05 3.2\n1 2 0.12 3.3 0.06 3.1\n"
        "2 2 0.15 3.1 0.07 3.0\n"
    )


def _matching_gro() -> str:
    coordinates = (
        (0.100, 0.100, 0.100),
        (0.250, 0.100, 0.100),
        (0.400, 0.100, 0.100),
        (0.550, 0.100, 0.100),
        (0.700, 0.100, 0.100),
        (0.850, 0.100, 0.100),
        (0.100, 0.100, 2.200),
    )
    atoms = "\n".join(
        f"    1BGLC    C1{index:5d}{x:8.3f}{y:8.3f}{z:8.3f}"
        for index, (x, y, z) in enumerate(coordinates, start=1)
    )
    return f"shared\n7\n{atoms}\n3.0 3.0 3.0 0.0 0.0 0.15 0.0 0.0 0.0\n"


def test_run_lammps_single_point_checks_and_normalizes_shared_coordinates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_file = tmp_path / "in.cellulose"
    input_file.write_text(_PLAIN_INPUT, encoding="utf-8")
    data_file = tmp_path / "system.data"
    data_file.write_text(_data(), encoding="utf-8")
    reference = tmp_path / "coordinate_reference.json"
    checksum = write_coordinate_reference(
        reference,
        [[1, 1, 1], [2.5, 1, 1], [4, 1, 1], [5.5, 1, 1], [7, 1, 1], [8.5, 1, 1], [1, 1, 22]],
        [[30, 0, 0], [1.5, 30, 0], [0, 0, 30]],
    )
    log = _log(_PLAIN_LABELS, _PLAIN_ROW)
    monkeypatch.setattr(
        "cellulose_forge.energy_drivers.lammps.subprocess.run",
        lambda *args, **kwargs: SimpleNamespace(returncode=0, stdout=log, stderr=""),
    )

    report = run_lammps_single_point(
        tmp_path,
        input_file=input_file,
        data_file=data_file,
        coordinate_reference_file=reference,
        coordinate_sha256=checksum,
        periodic_improper_count=0,
    )

    assert report.terms["potential"] == pytest.approx(_PLAIN_ROW["PotEng"] * KJ_PER_KCAL)
    assert (tmp_path / "single_point_lammps.log").read_text() == log


def test_run_lammps_single_point_requires_checked_read_data_target(tmp_path: Path) -> None:
    input_file = tmp_path / "in.cellulose"
    input_file.write_text(_PLAIN_INPUT.replace("read_data system.data", "read_data other.data"))
    data_file = tmp_path / "system.data"
    data_file.write_text(_data())
    (tmp_path / "other.data").write_text(_data().replace("1.0 1.0 1.0", "1.3 1.0 1.0", 1))
    reference = tmp_path / "coordinate_reference.json"
    checksum = write_coordinate_reference(reference, [[0.0, 0.0, 0.0]], np.eye(3))

    with pytest.raises(ValueError, match="read_data does not name"):
        run_lammps_single_point(
            tmp_path,
            input_file=input_file,
            data_file=data_file,
            coordinate_reference_file=reference,
            coordinate_sha256=checksum,
            periodic_improper_count=0,
        )


@pytest.mark.parametrize(
    "mutation",
    [
        "displace_atoms all random 0.4 0.4 0.4 12345",
        '"displace_atoms" all random 0.4 0.4 0.4 4928',
        "'set' atom 1 x 2.0",
        "set atom 1 x 2.0",
        "change_box all x scale 1.1",
        "read_dump altered.dump 0 x y z box yes replace yes",
        "replicate 2 1 1",
        "create_atoms 1 single 2 2 2",
        "delete_atoms group all",
        "include override.in",
        "jump SELF altered",
    ],
)
def test_run_lammps_single_point_refuses_coordinate_mutation_commands(
    tmp_path: Path, mutation: str
) -> None:
    input_file = tmp_path / "in.cellulose"
    input_file.write_text(_PLAIN_INPUT.replace("run 0", f"{mutation}\nrun 0"))
    data_file = tmp_path / "system.data"
    data_file.write_text(_data())
    reference = tmp_path / "coordinate_reference.json"
    checksum = write_coordinate_reference(
        reference,
        [[1, 1, 1], [2.5, 1, 1], [4, 1, 1], [5.5, 1, 1], [7, 1, 1], [8.5, 1, 1], [1, 1, 22]],
        [[30, 0, 0], [1.5, 30, 0], [0, 0, 30]],
    )

    with pytest.raises(ValueError, match="outside the verified subset"):
        run_lammps_single_point(
            tmp_path,
            input_file=input_file,
            data_file=data_file,
            coordinate_reference_file=reference,
            coordinate_sha256=checksum,
            periodic_improper_count=0,
        )


@pytest.mark.parametrize(
    ("replacement", "message"),
    [
        ("read_data system.data shift 1 0 0", "may not use keywords"),
        ("boundary s p p", "requires boundary p p p"),
    ],
)
def test_run_lammps_single_point_refuses_identity_changing_loader_options(
    tmp_path: Path, replacement: str, message: str
) -> None:
    input_file = tmp_path / "in.cellulose"
    input_text = _PLAIN_INPUT.replace("read_data system.data", replacement)
    if replacement.startswith("boundary"):
        input_text = input_text.replace("atom_style full", f"atom_style full\n{replacement}")
    input_file.write_text(input_text)
    data_file = tmp_path / "system.data"
    data_file.write_text(_data())
    reference = tmp_path / "coordinate_reference.json"
    checksum = write_coordinate_reference(
        reference,
        [[1, 1, 1], [2.5, 1, 1], [4, 1, 1], [5.5, 1, 1], [7, 1, 1], [8.5, 1, 1], [1, 1, 22]],
        [[30, 0, 0], [1.5, 30, 0], [0, 0, 30]],
    )

    with pytest.raises(ValueError, match=message):
        run_lammps_single_point(
            tmp_path,
            input_file=input_file,
            data_file=data_file,
            coordinate_reference_file=reference,
            coordinate_sha256=checksum,
            periodic_improper_count=0,
        )


@pytest.mark.parametrize("mutation", ["displace_atoms all move 1 0 0", "set atom 1 x 2.0"])
def test_run_lammps_single_point_refuses_comment_continuation_bypass(
    tmp_path: Path, mutation: str
) -> None:
    input_file = tmp_path / "in.cellulose"
    input_file.write_text(
        _PLAIN_INPUT.replace(
            "kspace_style pppm 1.0e-5",
            f"kspace_style pppm 1.0e-5 & # accuracy\n{mutation}",
        )
    )
    data_file = tmp_path / "system.data"
    data_file.write_text(_data())
    reference = tmp_path / "coordinate_reference.json"
    checksum = write_coordinate_reference(
        reference,
        [[1, 1, 1], [2.5, 1, 1], [4, 1, 1], [5.5, 1, 1], [7, 1, 1], [8.5, 1, 1], [1, 1, 22]],
        [[30, 0, 0], [1.5, 30, 0], [0, 0, 30]],
    )

    with pytest.raises(ValueError, match="combine & continuation and comments"):
        run_lammps_single_point(
            tmp_path,
            input_file=input_file,
            data_file=data_file,
            coordinate_reference_file=reference,
            coordinate_sha256=checksum,
            periodic_improper_count=0,
        )


@pytest.mark.parametrize(
    "intervening",
    ["# accuracy", ""],
)
@pytest.mark.parametrize(
    "mutation",
    ["displace_atoms all random 0.4 0.4 0.4 4928", "set atom 1 x 2.0"],
)
def test_run_lammps_single_point_refuses_interrupted_continuation(
    tmp_path: Path, intervening: str, mutation: str
) -> None:
    input_file = tmp_path / "in.cellulose"
    input_file.write_text(
        _PLAIN_INPUT.replace(
            "kspace_style pppm 1.0e-5",
            f"kspace_style pppm 1.0e-5 &\n{intervening}\n{mutation}",
        )
    )
    data_file = tmp_path / "system.data"
    data_file.write_text(_data())
    reference = tmp_path / "coordinate_reference.json"
    checksum = write_coordinate_reference(
        reference,
        [[1, 1, 1], [2.5, 1, 1], [4, 1, 1], [5.5, 1, 1], [7, 1, 1], [8.5, 1, 1], [1, 1, 22]],
        [[30, 0, 0], [1.5, 30, 0], [0, 0, 30]],
    )

    with pytest.raises(ValueError, match="blank or comment line"):
        run_lammps_single_point(
            tmp_path,
            input_file=input_file,
            data_file=data_file,
            coordinate_reference_file=reference,
            coordinate_sha256=checksum,
            periodic_improper_count=0,
        )


@pytest.mark.parametrize(
    "extra",
    [
        "fix move_atoms all move linear 1.0 0.0 0.0 units box",
        "run 0 every 1 'set atom 1 x 2.0'",
    ],
)
def test_run_lammps_single_point_refuses_mutation_hidden_in_allowed_command(
    tmp_path: Path, extra: str
) -> None:
    input_file = tmp_path / "in.cellulose"
    input_file.write_text(_PLAIN_INPUT.replace("run 0", f"{extra}\nrun 0"))
    data_file = tmp_path / "system.data"
    data_file.write_text(_data())
    reference = tmp_path / "coordinate_reference.json"
    checksum = write_coordinate_reference(
        reference,
        [[1, 1, 1], [2.5, 1, 1], [4, 1, 1], [5.5, 1, 1], [7, 1, 1], [8.5, 1, 1], [1, 1, 22]],
        [[30, 0, 0], [1.5, 30, 0], [0, 0, 30]],
    )

    with pytest.raises(ValueError, match=r"run 0|spring/self"):
        run_lammps_single_point(
            tmp_path,
            input_file=input_file,
            data_file=data_file,
            coordinate_reference_file=reference,
            coordinate_sha256=checksum,
            periodic_improper_count=0,
        )


def _validator() -> ModuleType:
    sys.path.insert(0, str(_VALIDATOR.parent))
    spec = importlib.util.spec_from_file_location("validate_charmm36_engines", _VALIDATOR)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _charmm_report(**overrides: object):  # type: ignore[no-untyped-def]
    arguments: dict[str, object] = {
        "settings": parse_lammps_input(_CHARMM_INPUT),
        "coordinate_sha256": "abc",
        "periodic_improper_count": 0,
        "charmm_conventions": CharmmfswConventions(
            one_four_lj_offset_kj_mol=1.25, switch_shift_kj_mol=7.0
        ),
    }
    arguments.update(overrides)
    return normalize_lammps_energies(dict(_CHARMM_ROW), **arguments)  # type: ignore[arg-type]


def test_charmm_conventions_match_the_documented_validator() -> None:
    validator = _validator()
    data = _data()

    assert CHARMM_COULOMB_RATIO == validator._COULOMB_RATIO
    assert lammps.LAMMPS_CHARMM_COULOMB == validator._LAMMPS_CHARMM_COULOMB
    assert charmmfsw_one_four_offset_kj_mol(data, 10.0, 12.0) == pytest.approx(
        validator._charmmfsw_one_four_offset(data), rel=1e-12, abs=1e-15
    )
    shift = charmmfsw_switch_shift_kj_mol(data, 10.0, 12.0)
    assert shift == pytest.approx(validator._charmmfsw_short_range_offset(data), rel=1e-12)
    assert shift != 0.0


def test_charmm_corrections_remain_visible_with_raw_and_comparable_values() -> None:
    report = _charmm_report()
    raw = {label: value * KJ_PER_KCAL for label, value in _CHARMM_ROW.items()}
    corrections = {item["name"]: item for item in report.corrections or ()}

    assert list(corrections) == [
        "lammps_charmm_coulomb_constant",
        "charmmfsw_one_four_lj_offset",
        "charmmfsw_force_switch_potential_shift",
    ]
    coulomb = corrections["lammps_charmm_coulomb_constant"]
    raw_electrostatics = raw["E_coul"] + raw["E_long"]
    assert coulomb["kind"] == "scale"
    assert coulomb["factor"] == CHARMM_COULOMB_RATIO
    assert coulomb["raw_kj_mol"] == pytest.approx(raw_electrostatics)
    assert coulomb["comparable_kj_mol"] == pytest.approx(raw_electrostatics * CHARMM_COULOMB_RATIO)
    assert corrections["charmmfsw_one_four_lj_offset"]["correction_kj_mol"] == -1.25
    shift = corrections["charmmfsw_force_switch_potential_shift"]
    assert shift["raw_kj_mol"] == pytest.approx(raw["E_vdwl"] - 1.25)
    assert shift["comparable_kj_mol"] == pytest.approx(raw["E_vdwl"] - 1.25 - 7.0)

    terms = report.terms
    assert terms["electrostatics"] == pytest.approx(raw_electrostatics * CHARMM_COULOMB_RATIO)
    assert terms["vdw"] == pytest.approx(raw["E_vdwl"] - 8.25)
    total = sum(item["correction_kj_mol"] for item in corrections.values())
    assert terms["potential"] == pytest.approx(raw["PotEng"] + total)
    # Urey-Bradley stays in angle (angle_style charmm); impropers stay separate.
    assert terms["angle"] == pytest.approx(raw["E_angle"])
    assert terms["improper"] == pytest.approx(raw["E_impro"])
    assert terms["dihedral"] == pytest.approx(raw["E_dihed"])
    validate_potential_sum(report, absolute_kj_mol=1e-6, relative=0.0)

    extra = report.to_dict()
    assert extra["raw_thermo_kcal_mol"]["E_coul"] == _CHARMM_ROW["E_coul"]
    assert extra["raw_terms_kj_mol"]["vdw"] == pytest.approx(raw["E_vdwl"])
    assert extra["raw_terms_kj_mol"]["potential"] == pytest.approx(raw["PotEng"])


def test_full_report_from_log_input_and_data_files(tmp_path: Path) -> None:
    log = tmp_path / "in.cellulose.log"
    log.write_text(_log(_CHARMM_LABELS, _CHARMM_ROW), encoding="utf-8")
    data = _data()

    report = lammps_energy_report(
        log,
        input_text=_CHARMM_INPUT,
        coordinate_sha256="abc",
        periodic_improper_count=0,
        data_text=data,
    )

    corrections = {item["name"]: item for item in report.corrections or ()}
    assert corrections["charmmfsw_one_four_lj_offset"]["correction_kj_mol"] == pytest.approx(
        -charmmfsw_one_four_offset_kj_mol(data, 10.0, 12.0)
    )
    assert report.engine == "lammps"
    with pytest.raises(ValueError, match="require the LAMMPS data file"):
        lammps_energy_report(
            log, input_text=_CHARMM_INPUT, coordinate_sha256="abc", periodic_improper_count=0
        )


def test_plain_pair_style_splits_tail_and_restraint_without_corrections(tmp_path: Path) -> None:
    log = tmp_path / "in.cellulose.log"
    log.write_text(_log(_PLAIN_LABELS, _PLAIN_ROW), encoding="utf-8")

    report = lammps_energy_report(
        log, input_text=_PLAIN_INPUT, coordinate_sha256="abc", periodic_improper_count=0
    )

    terms = report.terms
    assert report.corrections == ()
    assert terms["dispersion_correction"] == pytest.approx(-5.0 * KJ_PER_KCAL)
    assert terms["vdw"] == pytest.approx((-20.0 + 5.0) * KJ_PER_KCAL)
    assert terms["restraint"] == pytest.approx(1.5 * KJ_PER_KCAL)
    assert terms["electrostatics"] == pytest.approx(100.0 * KJ_PER_KCAL)
    assert terms["potential"] == pytest.approx(_PLAIN_ROW["PotEng"] * KJ_PER_KCAL)
    validate_potential_sum(report, absolute_kj_mol=1e-6, relative=0.0)
    decomposition = report.to_dict()["decompositions"][0]
    assert decomposition["name"] == "lammps_evdwl_includes_etail"


def test_periodic_impropers_are_reclassified_from_the_topology() -> None:
    settings = parse_lammps_input(_PLAIN_INPUT)
    common = {"settings": settings, "coordinate_sha256": "abc"}

    with pytest.raises(ValueError, match="supply their separately evaluated energy"):
        normalize_lammps_energies(dict(_PLAIN_ROW), periodic_improper_count=4, **common)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="topology has none"):
        normalize_lammps_energies(
            dict(_PLAIN_ROW),
            periodic_improper_count=0,
            periodic_improper_kj_mol=1.0,
            **common,  # type: ignore[arg-type]
        )
    report = normalize_lammps_energies(
        dict(_PLAIN_ROW),
        periodic_improper_count=4,
        periodic_improper_kj_mol=12.0,
        **common,  # type: ignore[arg-type]
    )
    assert report.terms["improper"] == pytest.approx(12.0)
    assert report.terms["dihedral"] == pytest.approx(30.0 * KJ_PER_KCAL - 12.0)
    assert report.terms["potential"] == pytest.approx(_PLAIN_ROW["PotEng"] * KJ_PER_KCAL)
    validate_potential_sum(report, absolute_kj_mol=1e-6, relative=0.0)


def test_explicit_pair_table_must_be_split_out_of_bonds() -> None:
    text = _PLAIN_INPUT.replace(
        "bond_style harmonic", "bond_style hybrid harmonic table linear 9901"
    )
    settings = parse_lammps_input(text)
    common = {"settings": settings, "coordinate_sha256": "abc", "periodic_improper_count": 0}

    with pytest.raises(ValueError, match="tallied in ebond"):
        normalize_lammps_energies(dict(_PLAIN_ROW), **common)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="exactly vdw and electrostatics"):
        normalize_lammps_energies(dict(_PLAIN_ROW), explicit_pair_kj_mol={"vdw": 1.0}, **common)  # type: ignore[arg-type]
    report = normalize_lammps_energies(
        dict(_PLAIN_ROW),
        explicit_pair_kj_mol={"vdw": 3.0, "electrostatics": -40.0},
        **common,  # type: ignore[arg-type]
    )
    assert report.terms["bond"] == pytest.approx(100.0 * KJ_PER_KCAL + 37.0)
    assert report.terms["vdw"] == pytest.approx(-15.0 * KJ_PER_KCAL + 3.0)
    assert report.terms["electrostatics"] == pytest.approx(100.0 * KJ_PER_KCAL - 40.0)
    validate_potential_sum(report, absolute_kj_mol=1e-6, relative=0.0)
    with pytest.raises(ValueError, match="no pair table"):
        normalize_lammps_energies(
            dict(_PLAIN_ROW),
            settings=parse_lammps_input(_PLAIN_INPUT),
            coordinate_sha256="abc",
            periodic_improper_count=0,
            explicit_pair_kj_mol={"vdw": 0.0, "electrostatics": 0.0},
        )


@pytest.mark.parametrize(
    ("label", "value"),
    [("c_extra", 3.0), ("v_bias", -0.5), ("E_newterm", 1.0)],
)
def test_unknown_nonzero_columns_fail_closed(label: str, value: float) -> None:
    text = _CHARMM_INPUT.replace("etail eimp", f"etail eimp {label}")
    row = {**_CHARMM_ROW, label: value}

    with pytest.raises(ValueError, match="unsupported nonzero energy terms"):
        normalize_lammps_energies(
            row,
            settings=parse_lammps_input(text),
            coordinate_sha256="abc",
            periodic_improper_count=0,
            charmm_conventions=CharmmfswConventions(1.25, 7.0),
        )
    zero = normalize_lammps_energies(
        {**row, label: 0.0},
        settings=parse_lammps_input(text),
        coordinate_sha256="abc",
        periodic_improper_count=0,
        charmm_conventions=CharmmfswConventions(1.25, 7.0),
    )
    assert zero.terms["potential"] == _charmm_report().terms["potential"]


@pytest.mark.parametrize(
    ("edit", "error"),
    [
        (lambda text: text.replace("units real", "units metal"), "units 'metal'"),
        (lambda text: text + "thermo_modify norm yes\n", "normalized per atom"),
        (
            lambda text: text.replace("lj/charmmfsw/coul/long", "lj/charmm/coul/long"),
            "unsupported LAMMPS pair_style",
        ),
        (
            lambda text: text.replace("dihedral_style charmmfsw", "dihedral_style charmm"),
            "requires dihedral_style charmmfsw",
        ),
        (
            lambda text: text.replace("10.0 12.0", ""),
            "inner and outer switch distances",
        ),
    ],
)
def test_unsupported_settings_are_refused(edit, error: str) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(ValueError, match=error):
        _charmm_report(settings=parse_lammps_input(edit(_CHARMM_INPUT)))


def test_missing_required_columns_and_conventions_are_refused() -> None:
    no_impropers = {key: value for key, value in _CHARMM_ROW.items() if key != "E_impro"}
    text = _CHARMM_INPUT.replace(" etail eimp", " etail")
    with pytest.raises(ValueError, match="E_impro"):
        normalize_lammps_energies(
            no_impropers,
            settings=parse_lammps_input(text),
            coordinate_sha256="abc",
            periodic_improper_count=0,
            charmm_conventions=CharmmfswConventions(1.25, 7.0),
        )
    no_tail_column = _PLAIN_INPUT.replace("elong etail f_", "elong f_")
    with pytest.raises(ValueError, match="E_tail"):
        normalize_lammps_energies(
            {key: value for key, value in _PLAIN_ROW.items() if key != "E_tail"},
            settings=parse_lammps_input(no_tail_column),
            coordinate_sha256="abc",
            periodic_improper_count=0,
        )
    with pytest.raises(ValueError, match="thermo columns do not match"):
        _charmm_report(settings=parse_lammps_input(_PLAIN_INPUT), charmm_conventions=None)
    with pytest.raises(ValueError, match="requires the CHARMM convention energies"):
        _charmm_report(charmm_conventions=None)
    with pytest.raises(ValueError, match="apply only to lj/charmmfsw"):
        normalize_lammps_energies(
            dict(_PLAIN_ROW),
            settings=parse_lammps_input(_PLAIN_INPUT),
            coordinate_sha256="abc",
            periodic_improper_count=0,
            charmm_conventions=CharmmfswConventions(1.25, 7.0),
        )
    with pytest.raises(ValueError, match="not finite"):
        _charmm_report(charmm_conventions=CharmmfswConventions(math.nan, 7.0))


def test_thermo_parser_requires_one_single_point_row() -> None:
    text = _log(_CHARMM_LABELS, _CHARMM_ROW)
    parsed = parse_lammps_thermo(text)
    assert parsed["PotEng"] == pytest.approx(_CHARMM_ROW["PotEng"])

    row = text.splitlines()[4]
    with pytest.raises(ValueError, match="one thermo row"):
        parse_lammps_thermo(text.replace(row, f"{row}\n{row}"))
    with pytest.raises(ValueError, match="one thermo block"):
        parse_lammps_thermo(text + text)
    with pytest.raises(ValueError, match="does not match its header"):
        parse_lammps_thermo(text.replace(row, row.rsplit(maxsplit=1)[0]))
    with pytest.raises(ValueError, match="non-finite"):
        parse_lammps_thermo(text.replace(row, row.replace(row.split()[3], "nan", 1)))


def test_input_parser_reads_continuations_and_energy_fixes() -> None:
    settings = parse_lammps_input(_PLAIN_INPUT)

    assert settings.energy_fixes == ("cellulose_position_restraints",)
    assert settings.thermo_labels[-1] == "f_cellulose_position_restraints"
    assert settings.tail_correction is True
    assert settings.kspace is True
    assert settings.charmm is False
    assert parse_lammps_input(_CHARMM_INPUT).charmm is True
