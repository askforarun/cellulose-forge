"""Normalize GROMACS single-point energies into canonical energy reports."""

from __future__ import annotations

import math
from collections.abc import Mapping
from pathlib import Path

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
