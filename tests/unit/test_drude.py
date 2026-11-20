from __future__ import annotations

import json
import os
import shutil
import struct
import subprocess
from pathlib import Path

import numpy as np
import pytest
from drude_fixtures import MAIN, carbohydrate

from cellulose_forge.cli import main
from cellulose_forge.drude import (
    DrudeSystem,
    _write_run_configuration,
    lone_pair_position,
    render_drude_namd_configuration,
    render_drude_openmm_runner,
    render_drude_pdb,
    render_drude_psf,
    validate_drude_toppar_directory,
    write_namd_binary_coordinates,
)
from cellulose_forge.drude_topology import (
    DRUDE_CARBOHYDRATE_FILE,
    DRUDE_MAIN_FILE,
    build_drude_cellulose_chain,
    parse_drude_cellulose_definitions,
)
from cellulose_forge.polymer_profiles import capability_error, engine_capability_error
from cellulose_forge.run_inputs import dispersion_correction


def _angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    u, v = a - b, c - b
    return float(np.degrees(np.arccos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v))))


def _dihedral(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> float:
    b1, b2, b3 = b - a, c - b, d - c
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    # IUPAC sign: positive for a clockwise rotation looking along b -> c.
    y = np.dot(b2 / np.linalg.norm(b2), np.cross(n1, n2))
    return float(np.degrees(np.arctan2(y, np.dot(n1, n2))) % 360.0)


HOST = np.array([0.0, 0.0, 0.0])
SECOND = np.array([1.4, 0.1, -0.2])
THIRD = np.array([1.9, 1.2, 0.4])


@pytest.mark.parametrize(("angle", "dihedral"), [(110.0, 91.0), (110.0, 269.0), (95.0, 10.0)])
def test_relative_lone_pair_reproduces_its_internal_coordinates(
    angle: float, dihedral: float
) -> None:
    site = lone_pair_position("relative", HOST, SECOND, THIRD, 0.35, angle, dihedral)
    assert np.linalg.norm(site - HOST) == pytest.approx(0.35)
    assert _angle(site, HOST, SECOND) == pytest.approx(angle)
    assert _dihedral(site, HOST, SECOND, THIRD) == pytest.approx(dihedral)


def test_bisector_lone_pair_uses_the_midpoint_of_the_two_other_hosts() -> None:
    midpoint = 0.5 * (SECOND + THIRD)
    site = lone_pair_position("bisector", HOST, SECOND, THIRD, 0.35, 110.0, 90.0)
    assert _angle(site, HOST, midpoint) == pytest.approx(110.0)
    assert _dihedral(site, HOST, midpoint, THIRD) == pytest.approx(90.0)
    axis = lone_pair_position("bisector", HOST, SECOND, THIRD, 0.10, 0.0, 0.0)
    assert axis == pytest.approx(0.10 * midpoint / np.linalg.norm(midpoint))


def test_lone_pair_placement_refuses_degenerate_hosts_and_unknown_modes() -> None:
    with pytest.raises(ValueError, match="collinear"):
        lone_pair_position("relative", HOST, SECOND, 2.0 * SECOND, 0.35, 110.0, 90.0)
    with pytest.raises(ValueError, match="coincide"):
        lone_pair_position("relative", HOST, HOST, THIRD, 0.35, 110.0, 90.0)
    with pytest.raises(ValueError, match="unsupported lone-pair mode"):
        lone_pair_position("colinear", HOST, SECOND, THIRD, 0.35, 110.0, 90.0)


def _system(chains: int = 2, *, periodic: bool = False) -> DrudeSystem:
    definitions = parse_drude_cellulose_definitions(MAIN, carbohydrate())
    topology = build_drude_cellulose_chain(definitions, 2, periodic=periodic)
    coordinates = np.arange(chains * len(topology.sites) * 3, dtype=float).reshape(-1, 3) / 7.0
    return DrudeSystem(topology, tuple(range(1, chains + 1)), coordinates, np.eye(3) * 40.0)


def _section(lines: list[str], marker: str) -> tuple[list[int], list[str]]:
    start = next(index for index, line in enumerate(lines) if marker in line)
    counts = [int(value) for value in lines[start].split("!")[0].split()]
    body: list[str] = []
    for line in lines[start + 1 :]:
        if "!" in line:
            break
        body.append(line)
    return counts, body


def test_psf_has_drude_lone_pair_and_anisotropy_sections() -> None:  # noqa: PLR0915
    system = _system()
    topology = system.topology
    lines = render_drude_psf(system).splitlines()
    assert lines[0] == "PSF EXT DRUDE"
    (natom,), atom_rows = _section(lines, "!NATOM")
    rows = [row.split() for row in atom_rows if row.strip()]
    assert natom == len(rows) == 2 * len(topology.sites)
    assert {row[1] for row in rows} == {"C001", "C002"}
    by_kind = {"atom": "0", "lonepair": "-1", "drude": "-2"}
    for row, site in zip(rows, topology.sites * 2, strict=True):
        assert row[4:6] == [site.name, site.atom_type]
        assert float(row[6]) == pytest.approx(site.charge, abs=5e-7)
        assert row[8] == by_kind[site.kind]
        assert float(row[9]) == pytest.approx(site.alpha) and float(row[10]) == site.thole
    assert sum(float(row[6]) for row in rows) == pytest.approx(0.0, abs=1e-5)

    for marker, items, size in (
        ("!NTHETA", topology.angles, 3),
        ("!NPHI", topology.dihedrals, 4),
    ):
        (count,), body = _section(lines, marker)
        values = [int(value) for row in body for value in row.split()]
        assert count == 2 * len(items) and len(values) == count * size
        assert values[: len(items) * size] == [index + 1 for item in items for index in item]
        assert min(values) >= 1 and max(values) <= natom

    # Bonds: every chain's ordinary and lone-pair bonds, then every chain's
    # parent--Drude bonds.  NAMD 3.0.3 cannot read the PSF otherwise.
    (bond_count,), bond_rows = _section(lines, "!NBOND")
    values = [int(value) for row in bond_rows for value in row.split()]
    bonds = list(zip(values[::2], values[1::2], strict=True))
    assert bond_count == len(bonds) == 2 * len(topology.bonds)
    flags = [rows[first - 1][8] == "-2" or rows[second - 1][8] == "-2" for first, second in bonds]
    drude_count = 2 * len(topology.drude_parents)
    assert flags == [False] * (len(bonds) - drude_count) + [True] * drude_count
    per_chain = len(topology.sites)
    expected = {
        (first + offset + 1, second + offset + 1)
        for offset in (0, per_chain)
        for first, second in topology.bonds
    }
    assert set(bonds) == expected and len(set(bonds)) == len(bonds)
    for parent, drude in bonds[-drude_count:]:
        assert drude == parent + 1 and rows[parent - 1][8] == "0"

    (numlp, numlph), lone_pair_rows = _section(lines, "!NUMLP")
    assert numlp == 2 * len(topology.lone_pairs) and numlph == 4 * numlp
    records = [row.split() for row in lone_pair_rows if " F" in row]
    assert len(records) == numlp
    assert [int(record[1]) for record in records] == list(range(1, 4 * numlp, 4))
    for record, lone_pair in zip(records, topology.lone_pairs * 2, strict=True):
        sign = -1.0 if lone_pair.mode == "bisector" else 1.0
        assert record[0] == "3" and record[2] == "F"
        assert float(record[3]) == pytest.approx(sign * lone_pair.distance)
        assert float(record[4]) == lone_pair.angle and float(record[5]) == lone_pair.dihedral
    hosts = [int(v) for row in lone_pair_rows if " F" not in row for v in row.split()]
    first = topology.lone_pairs[0]
    assert hosts[:4] == [first.site + 1, *(host + 1 for host in first.hosts)]
    assert len(hosts) == numlph

    (naniso,), anisotropy_rows = _section(lines, "!NUMANISO")
    assert naniso == 2 * len(topology.anisotropies)
    constants = [row.split() for row in anisotropy_rows if "." in row]
    assert len(constants) == naniso
    assert [float(value) for value in constants[0]] == pytest.approx(
        [topology.anisotropies[0].k11, topology.anisotropies[0].k22, topology.anisotropies[0].k33],
        abs=5e-5,
    )
    quads = [
        int(v) for row in anisotropy_rows if row.strip() and "." not in row for v in row.split()
    ]
    assert quads[:4] == [index + 1 for index in topology.anisotropies[0].atoms]
    assert len(quads) == 4 * naniso


def test_coordinates_are_written_at_full_precision(tmp_path: Path) -> None:
    system = _system(chains=1)
    path = tmp_path / "structure.coor"
    write_namd_binary_coordinates(path, system.cartesian_coordinates)
    raw = path.read_bytes()
    assert struct.unpack("=i", raw[:4])[0] == system.atom_count
    assert np.array_equal(
        np.frombuffer(raw[4:], dtype="=f8").reshape(-1, 3), system.cartesian_coordinates
    )
    pdb = [line for line in render_drude_pdb(system).splitlines() if line.startswith("ATOM")]
    assert len(pdb) == system.atom_count
    assert len({len(line) for line in pdb}) == 1
    assert all(len(line[12:16]) == 4 for line in pdb)
    pdb_coordinates = np.array(
        [[float(line[30:38]), float(line[38:46]), float(line[46:54])] for line in pdb]
    )
    assert pdb_coordinates == pytest.approx(system.cartesian_coordinates, abs=5.1e-4)
    # Five-character CHARMM names stay intact in the PSF but never shift PDB columns.
    assert any(len(site.name) > 4 for site in system.topology.sites)


def test_construction_only_output_removes_stale_input_and_marks_periodic_psf(
    tmp_path: Path,
) -> None:
    finite = _system(chains=1)
    configuration = _write_run_configuration(
        tmp_path, finite, lj_correction=False, write_run_input=True
    )
    assert configuration is not None and configuration.is_file()
    periodic = _system(chains=1, periodic=True)
    assert (
        _write_run_configuration(tmp_path, periodic, lj_correction=False, write_run_input=False)
        is None
    )
    assert not configuration.exists()
    assert "PERIODIC CHAINS ARE CONSTRUCTION ONLY; DO NOT RUN" in render_drude_psf(periodic)


@pytest.mark.parametrize("lj_correction", [True, False])
def test_namd_configuration_uses_drude_settings(lj_correction: bool) -> None:
    text = render_drude_namd_configuration(np.eye(3) * 40.0, lj_correction=lj_correction)
    lines = text.splitlines()
    for expected in (
        "bincoordinates structure_namd.coor",
        "drude on",
        "drudeHardWall on",
        "vdwForceSwitching off",
        "exclude scaled1-4",
        "1-4scaling 1.0",
        "timestep 1.0",
        "rigidBonds water",
        "useSettle on",
        f"LJcorrection {'on' if lj_correction else 'off'}",
    ):
        assert expected in lines
    assert "rigidBonds none" not in lines
    assert "vmd" not in text.lower() and "psfgen" not in text.lower()


def test_openmm_runner_uses_drude_integrator_and_native_bundle() -> None:
    text = render_drude_openmm_runner(np.eye(3) * 40.0, lj_correction=False)
    compile(text, "run_drude_openmm.py", "exec")
    for expected in (
        'CharmmPsfFile("system.psf")',
        "CharmmParameterSet",
        'open("structure_namd.coor", "rb")',
        "DrudeLangevinIntegrator",
        "setMaxDrudeDistance(0.025*unit.nanometer)",
        "constraints=None",
        "ewaldErrorTolerance=1.0e-6",
        "system.setDefaultPeriodicBoxVectors(*vectors)",
        "system.getNumParticles() != count",
        "computeVirtualSites()",
    ):
        assert expected in text


def test_capabilities_and_long_range_correction_policy() -> None:
    assert capability_error("cellulose", "charmm-drude", "slab", ("namd",), "finite") is None
    assert capability_error("cellulose", "charmm-drude", "fibril", ("namd",), "periodic") is None
    assert capability_error("cellulose", "charmm-drude", "bulk", ("namd",), "periodic") is None
    assert "require periodic chains" in str(
        capability_error("cellulose", "charmm-drude", "bulk", ("namd",), "finite")
    )
    assert capability_error("cellulose", "charmm-drude", "unit-cell", ("namd",), "periodic")
    assert capability_error("chitin", "charmm-drude", "slab", ("namd",), "finite")
    assert engine_capability_error("cellulose", "charmm-drude", ("openmm",)) is None
    for engine in ("gromacs", "lammps", "amber"):
        assert engine_capability_error("cellulose", "charmm-drude", (engine,))
    assert dispersion_correction("charmm-drude", "bulk") is True
    assert dispersion_correction("charmm-drude", "slab") is False
    assert dispersion_correction("charmm-drude", "fibril") is False


def test_toppar_directory_must_hold_the_audited_streams(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="is missing"):
        validate_drude_toppar_directory(tmp_path)
    (tmp_path / DRUDE_MAIN_FILE).write_text(MAIN, encoding="utf-8")
    (tmp_path / DRUDE_CARBOHYDRATE_FILE).write_text(carbohydrate(), encoding="utf-8")
    with pytest.raises(ValueError, match="not the audited"):
        validate_drude_toppar_directory(tmp_path)


_SLAB = [
    "build", "--allomorph", "I-beta", "--shape", "slab", "--face", "1", "1", "0",
    "--layers", "3", "--chains-per-layer", "6", "--vacuum", "30",
    "--force-field", "charmm-drude", "--md-engine", "namd",
]  # fmt: skip


def test_cli_requires_the_stream_directory_and_refuses_unaudited_files(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    finite = [
        *_SLAB,
        "--chain-boundary", "finite",
        "--degree-polymerization", "6",
        "--padding-longitudinal", "20",
    ]  # fmt: skip
    with pytest.raises(SystemExit):
        main([*finite, "--output", str(tmp_path / "out")])
    assert "requires --drude-toppar-directory" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        main([*finite, "--vmd-executable", "vmd", "--output", str(tmp_path / "out")])
    assert "unrecognized arguments: --vmd-executable" in capsys.readouterr().err
    toppar = tmp_path / "toppar"
    toppar.mkdir()
    (toppar / DRUDE_MAIN_FILE).write_text(MAIN, encoding="utf-8")
    (toppar / DRUDE_CARBOHYDRATE_FILE).write_text(carbohydrate(), encoding="utf-8")
    with pytest.raises(SystemExit):
        main([*finite, "--drude-toppar-directory", str(toppar), "--output", str(tmp_path / "out")])
    assert "not the audited" in capsys.readouterr().err
    # Stream validation happens before anything is written.
    assert not (tmp_path / "out").exists()


def test_cli_refuses_periodic_chains_without_the_construction_only_option(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    periodic = [*_SLAB, "--chain-boundary", "periodic", "--degree-polymerization", "8"]
    with pytest.raises(SystemExit):
        main(
            [*periodic, "--drude-toppar-directory", str(tmp_path), "--output", str(tmp_path / "o")]
        )
    error = capsys.readouterr().err
    assert "construction-only" in error and "--drude-construction-only" in error
    assert not (tmp_path / "o").exists()


def test_cli_refuses_drude_construction_only_for_other_force_fields(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    arguments = [
        "build",
        "--allomorph",
        "I-beta",
        "--shape",
        "unit-cell",
        "--force-field",
        "none",
        "--drude-construction-only",
        "--output",
        str(tmp_path / "out.pdb"),
    ]
    with pytest.raises(SystemExit):
        main(arguments)
    assert "requires --force-field charmm-drude" in capsys.readouterr().err


_TOPPAR = os.environ.get("CELLULOSE_FORGE_DRUDE_TOPPAR")


@pytest.mark.skipif(
    _TOPPAR is None, reason="set CELLULOSE_FORGE_DRUDE_TOPPAR to the official streams"
)
@pytest.mark.parametrize(
    ("boundary", "extra", "per_chain"),
    [
        ("finite", ["--degree-polymerization", "6", "--padding-longitudinal", "20"], 269),
        ("periodic", ["--degree-polymerization", "8", "--drude-construction-only"], 352),
    ],
)
def test_cli_writes_a_native_drude_bundle(
    tmp_path: Path, boundary: str, extra: list[str], per_chain: int
) -> None:
    assert _TOPPAR is not None
    output = tmp_path / "drude"
    arguments = [*_SLAB, *extra, "--chain-boundary", boundary]
    assert main([*arguments, "--drude-toppar-directory", _TOPPAR, "--output", str(output)]) == 0
    report = json.loads((output / "build_report.json").read_text(encoding="utf-8"))
    counts = report["topology_counts"]
    assert report["chain_boundary"] == boundary
    assert report["long_range_lj_correction"] is False
    assert "psfgen" not in json.dumps(report["outputs"])
    assert counts["particles"] == 18 * per_chain
    assert counts["particles"] == (
        counts["real_atoms"] + counts["drude_particles"] + counts["lone_pairs"]
    )
    psf = (output / "system.psf").read_text(encoding="utf-8")
    assert psf.startswith("PSF EXT DRUDE")
    assert f"{counts['particles']:10d} !NATOM" in psf
    raw = (output / "structure_namd.coor").read_bytes()
    assert struct.unpack("=i", raw[:4])[0] == counts["particles"]
    coordinates = np.frombuffer(raw[4:], dtype="=f8").reshape(-1, 3)
    assert np.isfinite(coordinates).all()
    for name in ("structure_namd.pdb", DRUDE_MAIN_FILE):
        assert (output / name).is_file()
    assert not (output / "build_drude_psf.tcl").exists()
    # Periodic chains are construction-only: no run input is written for them.
    construction_only = boundary == "periodic"
    assert report["construction_only"] is construction_only
    assert (output / "run_drude_namd.conf").exists() is not construction_only
    assert (report["outputs"]["namd_configuration"] is None) is construction_only


_NAMD = os.environ.get("CELLULOSE_FORGE_NAMD") or shutil.which("namd3")
_OPENMM_PYTHON = os.environ.get("CELLULOSE_FORGE_OPENMM_PYTHON")


@pytest.mark.skipif(
    _TOPPAR is None or _NAMD is None,
    reason="needs CELLULOSE_FORGE_DRUDE_TOPPAR and a namd3 executable (or CELLULOSE_FORGE_NAMD)",
)
def test_namd_reads_and_runs_the_native_finite_bundle(tmp_path: Path) -> None:
    """Real NAMD regression: the PSF must load and a short run must complete."""
    assert _TOPPAR is not None and _NAMD is not None
    output = tmp_path / "drude"
    arguments = [
        *_SLAB,
        "--chain-boundary", "finite",
        "--degree-polymerization", "6",
        "--padding-longitudinal", "20",
    ]  # fmt: skip
    assert main([*arguments, "--drude-toppar-directory", _TOPPAR, "--output", str(output)]) == 0
    counts = json.loads((output / "build_report.json").read_text(encoding="utf-8"))[
        "topology_counts"
    ]
    configuration = (output / "run_drude_namd.conf").read_text(encoding="utf-8")
    assert "minimize 1000" in configuration and "run 1000" in configuration
    # NAMD requires minimization and dynamics lengths to be multiples of its
    # default 20-step cycle.
    short = (
        configuration.replace("minimize 1000", "minimize 100")
        .replace("run 1000", "run 20")
        .replace("outputEnergies 100", "outputEnergies 20")
    )
    (output / "short.conf").write_text(short, encoding="utf-8")
    result = subprocess.run(
        [_NAMD, "+p2", "short.conf"],
        cwd=output,
        capture_output=True,
        text=True,
        check=False,
        timeout=900,
    )
    log = result.stdout + result.stderr
    assert "FATAL ERROR" not in log, log[-2000:]
    assert result.returncode == 0, log[-2000:]
    assert f"Ignored {counts['lone_pairs']} bonds with lone pairs" in log
    assert f"Ignored {counts['drude_particles']} bonds with Drude particles" in log
    energies = [line.split() for line in log.splitlines() if line.startswith("ENERGY:")]
    assert energies
    assert energies[-1][1] == "120"
    assert all(np.isfinite(float(value)) for value in energies[-1][2:])


@pytest.mark.skipif(
    _TOPPAR is None or _OPENMM_PYTHON is None,
    reason="needs CELLULOSE_FORGE_DRUDE_TOPPAR and CELLULOSE_FORGE_OPENMM_PYTHON",
)
@pytest.mark.parametrize("face", [(1, 1, 0), (1, 0, 0)])
def test_openmm_reads_and_runs_the_native_finite_bundle(
    tmp_path: Path, face: tuple[int, int, int]
) -> None:
    """Real OpenMM regression: load the native PSF and complete 1 ps dynamics."""
    assert _TOPPAR is not None and _OPENMM_PYTHON is not None
    output = tmp_path / "drude"
    geometry = list(_SLAB[:-2])
    face_index = geometry.index("--face") + 1
    geometry[face_index : face_index + 3] = [str(value) for value in face]
    arguments = [
        *geometry,
        "--md-engine", "openmm",
        "--chain-boundary", "finite",
        "--degree-polymerization", "6",
        "--padding-longitudinal", "20",
    ]  # fmt: skip
    assert main([*arguments, "--drude-toppar-directory", _TOPPAR, "--output", str(output)]) == 0
    report = json.loads((output / "build_report.json").read_text(encoding="utf-8"))
    assert report["md_engines"] == ["openmm"]
    assert report["outputs"]["openmm_runner"] == "run_drude_openmm.py"
    assert report["outputs"]["namd_configuration"] is None
    result = subprocess.run(
        [_OPENMM_PYTHON, "run_drude_openmm.py"],
        cwd=output,
        capture_output=True,
        text=True,
        check=False,
        timeout=900,
    )
    log = result.stdout + result.stderr
    assert result.returncode == 0, log[-2000:]
    assert "Initial potential energy:" in log
    assert (output / "openmm_drude.chk").is_file()
    assert (output / "openmm_drude.dcd").is_file()
    rows = (output / "openmm_drude.log").read_text(encoding="utf-8").splitlines()
    assert rows[-1].split(",", 1)[0].strip('"') == "1000"
    assert all("nan" not in row.lower() and "inf" not in row.lower() for row in rows[1:])
