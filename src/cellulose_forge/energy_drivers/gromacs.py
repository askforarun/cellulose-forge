"""Normalize GROMACS single-point energies into canonical energy reports."""

from __future__ import annotations

import math
import subprocess
from collections.abc import Mapping
from pathlib import Path

from cellulose_forge.energy_drivers.coordinates import (
    assert_coordinate_identity,
    load_coordinate_reference,
    read_gro_frame,
)
from cellulose_forge.energy_drivers.report import EnergyReport, EnergyValue

# GROMACS may report more than one functional form in a category.  Keeping the
# names explicit makes an upstream addition fail closed instead of disappearing
# into a broad string match.
_TERM_CATEGORIES: dict[str, str] = {
    # Bonds.
    "Bond": "bond",
    "G96Bond": "bond",
    "Morse": "bond",
    "Cubic Bonds": "bond",
    "FENE": "bond",
    "Tab. Bonds": "bond",
    # Angles, including CHARMM Urey--Bradley terms.
    "Angle": "angle",
    "G96Angle": "angle",
    "Restricted Angles": "angle",
    "Linear Angles": "angle",
    "Cross Bond-Bond": "angle",
    "Cross Bond-Angle": "angle",
    "U-B": "angle",
    "Quartic Angles": "angle",
    "Tab. Angles": "angle",
    # Proper and improper torsions remain separate.
    "Proper Dih.": "dihedral",
    "Ryckaert-Bell.": "dihedral",
    "Fourier Dih.": "dihedral",
    "Tab. Dih.": "dihedral",
    "Improper Dih.": "improper",
    "Per. Imp. Dih.": "improper",
    # Aggregate nonbonded categories include their 1--4 contributions.
    "LJ-14": "vdw",
    "LJ (SR)": "vdw",
    "LJ recip.": "vdw",
    "Coulomb-14": "electrostatics",
    "Coulomb (SR)": "electrostatics",
    "Coul. recip.": "electrostatics",
    # Restraints are visible rather than folded into bonded terms.
    "Position Rest.": "restraint",
    "Flat-bottom posres": "restraint",
    "Distance Rest.": "restraint",
    "Orientation Rest.": "restraint",
    "Angle Rest.": "restraint",
    "Dihedral Rest.": "restraint",
    "Disper. corr.": "dispersion_correction",
    "Potential": "potential",
}

_COMPONENT_TERMS = (
    "bond",
    "angle",
    "dihedral",
    "improper",
    "vdw",
    "electrostatics",
    "restraint",
    "dispersion_correction",
)

_SELECTION_TERMS = (
    *tuple(_TERM_CATEGORIES),
    # Selected deliberately so unsupported contributions cannot disappear
    # from the canonical decomposition.
    "CMAP",
    "Pull",
    "AWH",
)
_LEGEND_PART_COUNT = 3


def _finite_energy(value: object, term: str) -> EnergyValue:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise TypeError(f"GROMACS term {term!r} is not numeric")
    try:
        finite = math.isfinite(value)
    except OverflowError as error:
        raise ValueError(f"GROMACS term {term!r} is outside the supported range") from error
    if not finite:
        raise ValueError(f"GROMACS term {term!r} is not finite")
    return value


def _record_legend(line: str, legends: dict[int, str]) -> None:
    parts = line.split('"')
    prefix = parts[0].split()
    if (
        len(parts) < _LEGEND_PART_COUNT
        or len(prefix) < _LEGEND_PART_COUNT
        or not prefix[1].startswith("s")
        or prefix[2] != "legend"
    ):
        return
    try:
        index = int(prefix[1][1:])
    except ValueError:
        return
    if index in legends:
        raise ValueError(f"duplicate GROMACS XVG legend s{index}")
    legends[index] = parts[1]


def _parse_data_row(line: str, line_number: int) -> list[float]:
    try:
        return [float(value) for value in line.split()]
    except ValueError as error:
        raise ValueError(f"invalid GROMACS XVG data on line {line_number}") from error


def parse_gromacs_xvg(source: str | Path) -> dict[str, float]:
    """Parse a one-frame ``gmx energy`` XVG file or its text.

    A single-point adapter must not silently choose one row from a trajectory,
    so multiple data rows are rejected.
    """
    if isinstance(source, Path):
        text = source.read_text(encoding="utf-8")
    else:
        text = source

    legends: dict[int, str] = {}
    rows: list[list[float]] = []
    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("@"):
            _record_legend(line, legends)
            continue
        rows.append(_parse_data_row(line, line_number))

    if not legends:
        raise ValueError("GROMACS XVG contains no energy legends")
    expected_indices = list(range(len(legends)))
    if sorted(legends) != expected_indices:
        raise ValueError("GROMACS XVG legends are not contiguous from s0")
    if len(rows) != 1:
        raise ValueError(f"GROMACS single-point XVG must contain one data row, found {len(rows)}")
    row = rows[0]
    if len(row) != len(legends) + 1:
        raise ValueError("GROMACS XVG data column count does not match its energy legends")
    if not all(math.isfinite(value) for value in row):
        raise ValueError("GROMACS XVG contains a non-finite value")
    return {legends[index]: row[index + 1] for index in expected_indices}


def normalize_gromacs_energies(
    raw_terms: Mapping[str, EnergyValue], *, coordinate_sha256: str
) -> EnergyReport:
    """Return a complete canonical decomposition of GROMACS potential energy.

    Unknown zero-valued diagnostics are harmless, but every unknown nonzero
    contribution is rejected because otherwise the canonical decomposition
    would not describe the reported potential.
    """
    if not isinstance(raw_terms, Mapping) or not raw_terms:
        raise ValueError("GROMACS energy mapping is empty")

    grouped: dict[str, EnergyValue] = dict.fromkeys(_COMPONENT_TERMS, 0.0)
    potential: EnergyValue | None = None
    unsupported: list[str] = []
    for raw_name, raw_value in raw_terms.items():
        if not isinstance(raw_name, str) or not raw_name:
            raise TypeError("GROMACS energy term names must be non-empty strings")
        value = _finite_energy(raw_value, raw_name)
        category = _TERM_CATEGORIES.get(raw_name)
        if category is None:
            if value != 0:
                unsupported.append(raw_name)
            continue
        if category == "potential":
            potential = value
        else:
            grouped[category] += value

    if unsupported:
        raise ValueError(
            f"GROMACS output contains unsupported nonzero energy terms: {sorted(unsupported)}"
        )
    if potential is None:
        raise ValueError("GROMACS output is missing Potential")
    grouped["potential"] = potential
    return EnergyReport(
        engine="gromacs",
        coordinate_sha256=coordinate_sha256,
        terms=grouped,
        extra={"raw_terms": dict(raw_terms)},
    )


def _run(command: list[str], *, directory: Path, stdin: str | None = None) -> str:
    completed = subprocess.run(
        command,
        cwd=directory,
        input=stdin,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"GROMACS command failed ({' '.join(command)}):\n{completed.stdout}{completed.stderr}"
        )
    return completed.stdout + completed.stderr


def _validate_single_point_mdp(path: Path) -> None:
    settings: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split(";", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        settings[key.strip().lower()] = value.strip().lower()
    if settings.get("nsteps") != "0":
        raise ValueError("GROMACS single-point MDP must set nsteps = 0")
    if settings.get("constraints", "none") != "none":
        raise ValueError("GROMACS single-point MDP must set constraints = none")


def run_gromacs_single_point(  # noqa: PLR0913 - scientific inputs stay explicit
    directory: Path,
    *,
    coordinate_file: Path,
    coordinate_reference_file: Path,
    topology_file: Path,
    mdp_file: Path,
    coordinate_sha256: str,
    executable: str = "gmx",
) -> EnergyReport:
    """Run one GROMACS zero-step evaluation and return canonical energies.

    The caller supplies the MDP because cutoff, switching, PME, dispersion,
    periodic-molecule, and constraint settings must be matched across engines;
    the adapter must never invent those scientific settings.
    """
    directory = directory.resolve()
    coordinate_file = coordinate_file.resolve()
    coordinate_reference_file = coordinate_reference_file.resolve()
    topology_file = topology_file.resolve()
    mdp_file = mdp_file.resolve()
    for path, description in (
        (coordinate_file, "coordinate"),
        (coordinate_reference_file, "coordinate reference"),
        (topology_file, "topology"),
        (mdp_file, "MDP"),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"GROMACS {description} file does not exist: {path}")
    positions, box = read_gro_frame(coordinate_file)
    reference = load_coordinate_reference(coordinate_reference_file, coordinate_sha256)
    assert_coordinate_identity(reference, positions, engine="GROMACS", box_vectors=box)
    _validate_single_point_mdp(mdp_file)

    tpr = directory / "single_point.tpr"
    prefix = directory / "single_point"
    xvg = directory / "single_point_energy.xvg"
    _run(
        [
            executable,
            "grompp",
            "-f",
            str(mdp_file),
            "-c",
            str(coordinate_file),
            "-p",
            str(topology_file),
            "-o",
            str(tpr),
            "-maxwarn",
            "0",
        ],
        directory=directory,
    )
    _run(
        [
            executable,
            "mdrun",
            "-s",
            str(tpr),
            "-deffnm",
            str(prefix),
            "-ntmpi",
            "1",
            "-ntomp",
            "1",
        ],
        directory=directory,
    )
    selections = "\n".join(name.replace(" ", "-") for name in _SELECTION_TERMS)
    _run(
        [
            executable,
            "energy",
            "-f",
            str(prefix.with_suffix(".edr")),
            "-o",
            str(xvg),
        ],
        directory=directory,
        stdin=f"{selections}\n0\n",
    )
    return normalize_gromacs_energies(parse_gromacs_xvg(xvg), coordinate_sha256=coordinate_sha256)
