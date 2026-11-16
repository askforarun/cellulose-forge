"""Normalize NAMD CHARMM single-point energies into canonical energy reports.

Live NAMD evaluation must set ``mergeCrossterms no`` so a CHARMM CMAP
contribution is printed as ``CROSSTERM`` rather than folded into ``DIHED``.
This module does not write NAMD configuration.
"""

from __future__ import annotations

import hashlib
import math
import struct
import subprocess
from collections.abc import Mapping
from pathlib import Path

from cellulose_forge.energy_drivers.report import CANONICAL_TERMS, EnergyReport, EnergyValue

KJ_PER_KCAL = 4.184
_GRO_MINIMUM_LINES = 3
_NAMD_BINARY_HEADER_BYTES = 4
_COORDINATE_TOLERANCE_ANGSTROM = 1.0e-9

# NAMD prints one ETITLE row and then ENERGY rows whose columns follow those
# names.  Mapping is by name, never by fixed column index: CUDA/CPU builds and
# optional features (CMAP, extra constraints) can insert or reorder fields.
#
# CHARMM term semantics: ANGLE includes Urey--Bradley 1--3 terms; DIHED is
# proper torsions only; IMPRP is harmonic impropers; ELECT/VDW already include
# their scaled 1--4 contributions.  There is no separate CHARMM UB or 1--4
# column in a standard NAMD log.
#
# BOUNDARY (spherical boundary) and MISC (external field / steering / similar)
# are not automatically restraints.  Zero values are inactive; a nonzero value
# must be classified by the caller.
_TERM_CATEGORIES: dict[str, str] = {
    "BOND": "bond",
    "ANGLE": "angle",
    "DIHED": "dihedral",
    "IMPRP": "improper",
    "ELECT": "electrostatics",
    "VDW": "vdw",
    "POTENTIAL": "potential",
}
_CLASSIFIABLE_TERMS = frozenset({"BOUNDARY", "MISC"})
_CLASSIFICATION_TARGETS = CANONICAL_TERMS - {"potential"}

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


def parse_namd_energy_log(source: str | Path) -> dict[str, float]:  # noqa: PLR0912
    """Parse ETITLE-aligned ENERGY fields from a NAMD log or its text.

    Columns are identified by the ETITLE names, not by position. Exactly one
    ENERGY row is required for a single-point calculation.
    """
    if isinstance(source, Path):
        text = source.read_text(encoding="utf-8", errors="replace")
    else:
        text = source

    fields: list[str] | None = None
    rows: list[dict[str, float]] = []
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
        rows.append(values)
    if not rows:
        raise ValueError("NAMD log contains no ENERGY row")
    if len(rows) != 1:
        raise ValueError(f"NAMD single-point log must contain one ENERGY row, found {len(rows)}")
    return rows[0]


def _configuration_value(text: str, name: str) -> str | None:
    for raw_line in text.splitlines():
        line = raw_line.split("#", 1)[0].strip()
        fields = line.split()
        if fields and fields[0].lower() == name.lower():
            return " ".join(fields[1:]).lower()
    return None


def _gro_positions_angstrom(path: Path) -> list[tuple[float, float, float]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < _GRO_MINIMUM_LINES:
        raise ValueError("shared GRO coordinate file is truncated")
    try:
        count = int(lines[1].strip())
    except ValueError as error:
        raise ValueError("shared GRO atom count is invalid") from error
    if len(lines) < count + 3:
        raise ValueError("shared GRO coordinate file is truncated")
    try:
        return [
            (
                10.0 * float(line[20:28]),
                10.0 * float(line[28:36]),
                10.0 * float(line[36:44]),
            )
            for line in lines[2 : 2 + count]
        ]
    except ValueError as error:
        raise ValueError("shared GRO positions are invalid") from error


def _namd_binary_positions(path: Path) -> list[tuple[float, float, float]]:
    payload = path.read_bytes()
    if len(payload) < _NAMD_BINARY_HEADER_BYTES:
        raise ValueError("NAMD binary coordinate file is truncated")
    count = struct.unpack_from("=i", payload)[0]
    expected = _NAMD_BINARY_HEADER_BYTES + count * 3 * 8
    if count < 0 or len(payload) != expected:
        raise ValueError("NAMD binary coordinate size does not match its atom count")
    values = struct.unpack_from(f"={count * 3}d", payload, _NAMD_BINARY_HEADER_BYTES)
    return [
        (values[index], values[index + 1], values[index + 2]) for index in range(0, len(values), 3)
    ]


def run_namd_single_point(  # noqa: PLR0913
    directory: Path,
    *,
    configuration_file: Path,
    shared_coordinate_file: Path,
    binary_coordinate_file: Path,
    coordinate_sha256: str,
    executable: str = "namd3",
    threads: int = 1,
    classified_terms: Mapping[str, str] | None = None,
) -> EnergyReport:
    """Run a matched NAMD single point and return canonical CHARMM energies."""
    directory = directory.resolve()
    configuration_file = configuration_file.resolve()
    shared_coordinate_file = shared_coordinate_file.resolve()
    binary_coordinate_file = binary_coordinate_file.resolve()
    for path, description in (
        (configuration_file, "configuration"),
        (shared_coordinate_file, "shared coordinate"),
        (binary_coordinate_file, "binary coordinate"),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"NAMD {description} file does not exist: {path}")
    if type(threads) is not int or threads < 1:
        raise ValueError("NAMD thread count must be a positive integer")
    actual_sha256 = hashlib.sha256(shared_coordinate_file.read_bytes()).hexdigest()
    if actual_sha256 != coordinate_sha256:
        raise ValueError("NAMD shared coordinate file does not match the shared checksum")
    gro_positions = _gro_positions_angstrom(shared_coordinate_file)
    binary_positions = _namd_binary_positions(binary_coordinate_file)
    if len(gro_positions) != len(binary_positions) or any(
        abs(left - right) > _COORDINATE_TOLERANCE_ANGSTROM
        for gro, binary in zip(gro_positions, binary_positions, strict=True)
        for left, right in zip(gro, binary, strict=True)
    ):
        raise ValueError("NAMD binary coordinates differ from the shared GRO coordinates")

    configuration = configuration_file.read_text(encoding="utf-8")
    required = {
        "mergeCrossterms": "no",
        "LJcorrection": "off",
        "rigidBonds": "none",
        "run": "0",
    }
    for name, expected in required.items():
        if _configuration_value(configuration, name) != expected.lower():
            raise ValueError(f"NAMD single-point configuration must set {name} {expected}")
    completed = subprocess.run(
        [executable, f"+p{threads}", str(configuration_file)],
        cwd=directory,
        text=True,
        capture_output=True,
        check=False,
    )
    log = directory / "single_point_namd.log"
    log.write_text(completed.stdout + completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"NAMD single-point command failed with exit code {completed.returncode}"
        )
    return normalize_namd_energies(
        parse_namd_energy_log(log),
        coordinate_sha256=coordinate_sha256,
        merge_crossterms=False,
        lj_correction=False,
        classified_terms=classified_terms,
    )


def _classification_map(classified_terms: Mapping[str, str] | None) -> dict[str, str]:
    if classified_terms is None:
        return {}
    invalid_names = sorted(set(classified_terms) - _CLASSIFIABLE_TERMS)
    if invalid_names:
        raise ValueError(
            f"NAMD classified_terms may only name BOUNDARY or MISC, not {invalid_names}"
        )
    invalid_targets = sorted(set(classified_terms.values()) - _CLASSIFICATION_TARGETS)
    if invalid_targets:
        raise ValueError(f"NAMD classified_terms use unknown canonical terms: {invalid_targets}")
    return dict(classified_terms)


def _namd_category(
    raw_name: str,
    value: EnergyValue,
    classification: Mapping[str, str],
    unclassified: list[str],
    unsupported: list[str],
) -> str | None:
    if raw_name in _DIAGNOSTICS:
        return None
    if raw_name in classification:
        return classification[raw_name]
    if raw_name in _CLASSIFIABLE_TERMS:
        if value != 0:
            unclassified.append(raw_name)
        return None
    category = _TERM_CATEGORIES.get(raw_name)
    if category is None and value != 0:
        unsupported.append(raw_name)
    return category


def normalize_namd_energies(  # noqa: PLR0912
    raw_terms: Mapping[str, EnergyValue],
    *,
    coordinate_sha256: str,
    merge_crossterms: bool,
    lj_correction: bool,
    classified_terms: Mapping[str, str] | None = None,
) -> EnergyReport:
    """Return a complete canonical decomposition of NAMD CHARMM potential energy.

    Raw NAMD values are kcal/mol.  Canonical terms are kJ/mol.  NAMD does not
    report a long-range Lennard-Jones tail correction for the CHARMM settings
    CelluloseForge writes, so ``dispersion_correction`` is an explicit zero
    unless a future ETITLE field maps to it.  Unknown zero diagnostics are
    ignored; unknown nonzero potential-like contributions fail closed.

    ``BOUNDARY`` and ``MISC`` are classified only when the caller passes
    ``classified_terms`` (typically ``{"MISC": "restraint"}``).  Zero values
    are inactive.  Nonzero unclassified ``BOUNDARY``/``MISC`` fail closed.

    Callers must confirm ``mergeCrossterms no`` so CMAP cannot hide inside
    ``DIHED``, and ``LJcorrection off`` so an invisible tail correction cannot
    be folded into ``VDW``.
    """
    if not isinstance(raw_terms, Mapping) or not raw_terms:
        raise ValueError("NAMD energy mapping is empty")
    if type(merge_crossterms) is not bool:
        raise TypeError("merge_crossterms must be a boolean attestation")
    if type(lj_correction) is not bool:
        raise TypeError("lj_correction must be a boolean attestation")
    if merge_crossterms:
        raise ValueError("NAMD single-point comparison requires mergeCrossterms no")
    if lj_correction:
        raise ValueError("NAMD single-point comparison requires LJcorrection off")

    missing = [name for name in _REQUIRED_CHARMM_TERMS if name not in raw_terms]
    if missing:
        raise ValueError(f"NAMD output is missing required CHARMM energy terms: {missing}")
    classification = _classification_map(classified_terms)

    grouped_kcal: dict[str, EnergyValue] = dict.fromkeys(_COMPONENT_TERMS, 0.0)
    potential_kcal: EnergyValue | None = None
    unclassified: list[str] = []
    unsupported: list[str] = []
    for raw_name, raw_value in raw_terms.items():
        if not isinstance(raw_name, str) or not raw_name:
            raise TypeError("NAMD energy term names must be non-empty strings")
        value = _finite_energy(raw_value, raw_name)
        category = _namd_category(raw_name, value, classification, unclassified, unsupported)
        if category is None:
            continue
        if category == "potential":
            potential_kcal = value
        else:
            grouped_kcal[category] += value

    if unclassified:
        raise ValueError(
            f"NAMD output contains unclassified nonzero BOUNDARY/MISC terms: {sorted(unclassified)}"
        )
    if unsupported:
        raise ValueError(
            f"NAMD output contains unsupported nonzero energy terms: {sorted(unsupported)}"
        )
    if potential_kcal is None:
        raise ValueError("NAMD output is missing POTENTIAL")

    terms = {name: value * KJ_PER_KCAL for name, value in grouped_kcal.items()}
    terms["potential"] = potential_kcal * KJ_PER_KCAL
    extra: dict[str, object] = {
        "raw_terms": dict(raw_terms),
        "raw_units": "kcal/mol",
        "settings": {"merge_crossterms": False, "lj_correction": False},
    }
    if classification:
        extra["classified_terms"] = classification
    return EnergyReport(
        engine="namd",
        coordinate_sha256=coordinate_sha256,
        terms=terms,
        extra=extra,
    )
