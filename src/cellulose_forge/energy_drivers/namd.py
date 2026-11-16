"""Normalize NAMD CHARMM single-point energies into canonical energy reports."""

from __future__ import annotations

import math
from collections.abc import Mapping
from pathlib import Path

from cellulose_forge.energy_drivers.report import EnergyReport, EnergyValue

KJ_PER_KCAL = 4.184

# NAMD prints one ETITLE row and then ENERGY rows whose columns follow those
# names.  Mapping is by name, never by fixed column index: CUDA/CPU builds and
# optional features (CMAP, extra constraints) can insert or reorder fields.
#
# CHARMM term semantics: ANGLE includes Urey--Bradley 1--3 terms; DIHED is
# proper torsions only; IMPRP is harmonic impropers; ELECT/VDW already include
# their scaled 1--4 contributions.  There is no separate CHARMM UB or 1--4
# column in a standard NAMD log.
_TERM_CATEGORIES: dict[str, str] = {
    "BOND": "bond",
    "ANGLE": "angle",
    "DIHED": "dihedral",
    "IMPRP": "improper",
    "ELECT": "electrostatics",
    "VDW": "vdw",
    "BOUNDARY": "restraint",
    "MISC": "restraint",
    "POTENTIAL": "potential",
}

_REQUIRED_CHARMM_TERMS = (
    "BOND",
    "ANGLE",
    "DIHED",
    "IMPRP",
    "ELECT",
    "VDW",
    "POTENTIAL",
)

# Thermodynamic / mechanical columns that NAMD always prints beside the
# potential decomposition.  They may be nonzero without contributing to
# POTENTIAL, so they are diagnostics rather than unknown energy terms.
_DIAGNOSTICS = frozenset(
    {
        "TS",
        "KINETIC",
        "TOTAL",
        "TEMP",
        "TOTAL3",
        "TEMPAVG",
        "PRESSURE",
        "GPRESSURE",
        "VOLUME",
        "PRESSAVG",
        "GPRESSAVG",
    }
)

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


def _finite_energy(value: object, term: str) -> EnergyValue:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise TypeError(f"NAMD term {term!r} is not numeric")
    try:
        finite = math.isfinite(value)
    except OverflowError as error:
        raise ValueError(f"NAMD term {term!r} is outside the supported range") from error
    if not finite:
        raise ValueError(f"NAMD term {term!r} is not finite")
    return value


def parse_namd_energy_log(source: str | Path) -> dict[str, float]:
    """Parse ETITLE-aligned ENERGY fields from a NAMD log or its text.

    Columns are identified by the ETITLE names, not by position.  The first
    ENERGY row after an ETITLE is the single-point / first reported step.
    """
    if isinstance(source, Path):
        text = source.read_text(encoding="utf-8", errors="replace")
    else:
        text = source

    fields: list[str] | None = None
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("ETITLE:"):
            titles = stripped.split()[1:]
            if not titles:
                raise ValueError("NAMD ETITLE contains no fields")
            if len(titles) != len(set(titles)):
                raise ValueError("NAMD ETITLE contains duplicate field names")
            fields = titles
            continue
        if not stripped.startswith("ENERGY:"):
            continue
        if fields is None:
            raise ValueError("NAMD log contains ENERGY before ETITLE")
        tokens = stripped.split()[1:]
        if len(tokens) != len(fields):
            raise ValueError("NAMD ENERGY field count differs from ETITLE")
        values: dict[str, float] = {}
        for name, token in zip(fields, tokens, strict=True):
            try:
                value = float(token)
            except ValueError as error:
                raise ValueError(f"NAMD ENERGY field {name!r} is not numeric") from error
            if not math.isfinite(value):
                raise ValueError(f"NAMD ENERGY field {name!r} is not finite")
            values[name] = value
        return values
    raise ValueError("NAMD log contains no ENERGY row")


def normalize_namd_energies(
    raw_terms: Mapping[str, EnergyValue], *, coordinate_sha256: str
) -> EnergyReport:
    """Return a complete canonical decomposition of NAMD CHARMM potential energy.

    Raw NAMD values are kcal/mol.  Canonical terms are kJ/mol.  NAMD does not
    report a long-range Lennard-Jones tail correction for the CHARMM settings
    CelluloseForge writes, so ``dispersion_correction`` is an explicit zero
    unless a future ETITLE field maps to it.  Unknown zero diagnostics are
    ignored; unknown nonzero potential-like contributions fail closed.
    """
    if not isinstance(raw_terms, Mapping) or not raw_terms:
        raise ValueError("NAMD energy mapping is empty")

    missing = [name for name in _REQUIRED_CHARMM_TERMS if name not in raw_terms]
    if missing:
        raise ValueError(f"NAMD output is missing required CHARMM energy terms: {missing}")

    grouped_kcal: dict[str, EnergyValue] = dict.fromkeys(_COMPONENT_TERMS, 0.0)
    potential_kcal: EnergyValue | None = None
    unsupported: list[str] = []
    for raw_name, raw_value in raw_terms.items():
        if not isinstance(raw_name, str) or not raw_name:
            raise TypeError("NAMD energy term names must be non-empty strings")
        value = _finite_energy(raw_value, raw_name)
        if raw_name in _DIAGNOSTICS:
            continue
        category = _TERM_CATEGORIES.get(raw_name)
        if category is None:
            if value != 0:
                unsupported.append(raw_name)
            continue
        if category == "potential":
            potential_kcal = value
        else:
            grouped_kcal[category] += value

    if unsupported:
        raise ValueError(
            f"NAMD output contains unsupported nonzero energy terms: {sorted(unsupported)}"
        )
    if potential_kcal is None:
        raise ValueError("NAMD output is missing POTENTIAL")

    terms = {name: value * KJ_PER_KCAL for name, value in grouped_kcal.items()}
    terms["potential"] = potential_kcal * KJ_PER_KCAL
    return EnergyReport(
        engine="namd",
        coordinate_sha256=coordinate_sha256,
        terms=terms,
        extra={"raw_terms": dict(raw_terms), "raw_units": "kcal/mol"},
    )
