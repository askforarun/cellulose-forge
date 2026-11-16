"""Normalize finite-chain GLYCAM Amber energies into canonical energy reports.

Native Amber output in CelluloseForge is finite-chain GLYCAM only.  This
adapter does not accept CHARMM-specific Amber terms (Urey--Bradley, harmonic
impropers, CMAP) or rewrite sander inputs.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
from pathlib import Path

from cellulose_forge.energy_drivers.report import EnergyReport, EnergyValue

KJ_PER_KCAL = 4.184
# Documented GLYCAM Amber vs GROMACS Coulomb constants, in kcal A / mol / e^2.
# Canonical electrostatics stay on Amber's native convention; the GROMACS
# convention value is recorded as a correction so comparisons can apply it.
AMBER_COULOMB_KCAL_ANGSTROM = 332.0522
GROMACS_COULOMB_KCAL_ANGSTROM = 332.0637
COULOMB_SCALE_TO_GROMACS = GROMACS_COULOMB_KCAL_ANGSTROM / AMBER_COULOMB_KCAL_ANGSTROM

_ASSIGNMENT = re.compile(r"([A-Za-z0-9()][A-Za-z0-9() .-]*?)\s+=\s+([-+0-9.Ee]+)")

_TERM_CATEGORIES: dict[str, str] = {
    "BOND": "bond",
    "ANGLE": "angle",
    "DIHED": "dihedral",
    "1-4 NB": "vdw",
    "1-4 VDW": "vdw",
    "VDWAALS": "vdw",
    "1-4 EEL": "electrostatics",
    "EELEC": "electrostatics",
    "EEL": "electrostatics",
    "RESTRAINT": "restraint",
    "EDISPER": "dispersion_correction",
    "Edisper": "dispersion_correction",
    "EPtot": "potential",
}

_REQUIRED_GLYCAM_TERMS = ("BOND", "ANGLE", "DIHED", "VDWAALS", "EPtot")
_FOURTEEN_VDW = ("1-4 NB", "1-4 VDW")
_ELECTROSTATICS = ("EELEC", "EEL")
_CONFLICTING_ALIASES = (
    _FOURTEEN_VDW,
    _ELECTROSTATICS,
)

# Printed beside the potential but not part of its decomposition, or known
# GLYCAM-zero terms that fail closed when they become nonzero.
_DIAGNOSTICS = frozenset(
    {
        "NSTEP",
        "TIME(PS)",
        "TEMP(K)",
        "PRESS",
        "Etot",
        "EKtot",
        "EAMBER (non-restraint)",
        "EKCMT",
        "VIRIAL",
        "VOLUME",
        "Density",
    }
)
_UNSUPPORTED_IF_NONZERO = frozenset({"EHBOND", "UB", "IMP", "CMAP", "IMPROP", "CROSSTERM"})

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
        raise TypeError(f"Amber term {term!r} is not numeric")
    try:
        finite = math.isfinite(value)
    except OverflowError as error:
        raise ValueError(f"Amber term {term!r} is outside the supported range") from error
    if not finite:
        raise ValueError(f"Amber term {term!r} is not finite")
    return value


def _normalize_label(name: str) -> str:
    return " ".join(name.split())


def _parse_assignment(name: str, token: str) -> tuple[str, float]:
    key = _normalize_label(name)
    try:
        value = float(token)
    except ValueError as error:
        raise ValueError(f"Amber energy field {key!r} is not numeric") from error
    if not math.isfinite(value):
        raise ValueError(f"Amber energy field {key!r} is not finite")
    return key, value


def parse_amber_energy(source: str | Path) -> dict[str, float]:
    """Parse the first sander ``mdout`` / ``mdinfo`` potential-energy block.

    Assignments without spaces around ``=`` (namelist echoes such as ``imin=1``)
    are ignored.  Later average or step blocks are not mixed into the first
    single-point evaluation.
    """
    if isinstance(source, Path):
        text = source.read_text(encoding="utf-8")
    else:
        text = source

    blocks: list[dict[str, float]] = []
    current: dict[str, float] = {}

    def flush() -> None:
        nonlocal current
        if current:
            blocks.append(current)
            current = {}

    for line in text.splitlines():
        pairs = _ASSIGNMENT.findall(line)
        if not pairs:
            flush()
            continue
        for name, token in pairs:
            key, value = _parse_assignment(name, token)
            if (key == "NSTEP" and current) or key in current:
                flush()
            current[key] = value
    flush()

    for block in blocks:
        if "EPtot" in block:
            return block
    raise ValueError("Amber output contains no EPtot energy block")


def _coulomb_correction(
    electrostatics_kcal: EnergyValue,
) -> dict[str, object]:
    raw_kj = electrostatics_kcal * KJ_PER_KCAL
    return {
        "kind": "coulomb_constant",
        "term": "electrostatics",
        "raw_kcal_mol": electrostatics_kcal,
        "raw_kj_mol": raw_kj,
        "comparable_kj_mol": raw_kj * COULOMB_SCALE_TO_GROMACS,
        "engine_constant_kcal_angstrom_per_mol_e2": AMBER_COULOMB_KCAL_ANGSTROM,
        "reference_constant_kcal_angstrom_per_mol_e2": GROMACS_COULOMB_KCAL_ANGSTROM,
        "reference": "gromacs",
        "scale": COULOMB_SCALE_TO_GROMACS,
        "reason": (
            "Amber's Coulomb constant (332.0522 kcal A/mol/e^2) differs from "
            "GROMACS (332.0637). Canonical electrostatics remain on Amber's "
            "native convention; comparable_kj_mol is the GROMACS-convention value."
        ),
    }


def _require_glycam_terms(raw_terms: Mapping[str, EnergyValue]) -> None:
    if not isinstance(raw_terms, Mapping) or not raw_terms:
        raise ValueError("Amber energy mapping is empty")
    missing = [name for name in _REQUIRED_GLYCAM_TERMS if name not in raw_terms]
    if missing:
        raise ValueError(f"Amber output is missing required GLYCAM energy terms: {missing}")
    if not any(name in raw_terms for name in _FOURTEEN_VDW):
        raise ValueError("Amber output is missing 1-4 van der Waals (1-4 NB or 1-4 VDW)")
    if "1-4 EEL" not in raw_terms:
        raise ValueError("Amber output is missing 1-4 EEL")
    if not any(name in raw_terms for name in _ELECTROSTATICS):
        raise ValueError("Amber output is missing EELEC")
    for left, right in _CONFLICTING_ALIASES:
        if left in raw_terms and right in raw_terms:
            raise ValueError(
                f"Amber output contains conflicting energy labels: {left!r} and {right!r}"
            )


def _group_amber_terms(
    raw_terms: Mapping[str, EnergyValue],
    *,
    improper_kcal: EnergyValue,
) -> tuple[dict[str, EnergyValue], EnergyValue | None, EnergyValue, list[str]]:
    grouped_kcal: dict[str, EnergyValue] = dict.fromkeys(_COMPONENT_TERMS, 0.0)
    grouped_kcal["improper"] = improper_kcal
    potential_kcal: EnergyValue | None = None
    unsupported: list[str] = []
    electrostatics_kcal: EnergyValue = 0.0
    for raw_name, raw_value in raw_terms.items():
        if not isinstance(raw_name, str) or not raw_name:
            raise TypeError("Amber energy term names must be non-empty strings")
        value = _finite_energy(raw_value, raw_name)
        if raw_name in _DIAGNOSTICS:
            continue
        category = _TERM_CATEGORIES.get(raw_name)
        if raw_name in _UNSUPPORTED_IF_NONZERO or category is None:
            if value != 0:
                unsupported.append(raw_name)
            continue
        if category == "potential":
            potential_kcal = value
            continue
        if category == "dihedral":
            grouped_kcal["dihedral"] += value - improper_kcal
            continue
        grouped_kcal[category] += value
        if category == "electrostatics":
            electrostatics_kcal += value
    return grouped_kcal, potential_kcal, electrostatics_kcal, unsupported


def normalize_glycam_amber_energies(
    raw_terms: Mapping[str, EnergyValue],
    *,
    coordinate_sha256: str,
    periodic_improper: EnergyValue = 0,
) -> EnergyReport:
    """Return a complete canonical decomposition of finite-chain GLYCAM Amber energy.

    ``periodic_improper`` is the topology-classified periodic improper energy in
    kcal/mol.  Amber prints those impropers inside ``DIHED``; when the value is
    omitted it must be zero (cellulose).  Canonical electrostatics stay on
    Amber's Coulomb constant; the GROMACS-convention value is in ``corrections``.
    """
    _require_glycam_terms(raw_terms)
    improper_kcal = _finite_energy(periodic_improper, "periodic_improper")
    grouped_kcal, potential_kcal, electrostatics_kcal, unsupported = _group_amber_terms(
        raw_terms, improper_kcal=improper_kcal
    )
    if unsupported:
        raise ValueError(
            "Amber output contains unsupported nonzero energy terms for "
            f"finite-chain GLYCAM: {sorted(unsupported)}"
        )
    if potential_kcal is None:
        raise ValueError("Amber output is missing EPtot")

    terms = {name: value * KJ_PER_KCAL for name, value in grouped_kcal.items()}
    terms["potential"] = potential_kcal * KJ_PER_KCAL
    return EnergyReport(
        engine="amber",
        coordinate_sha256=coordinate_sha256,
        terms=terms,
        corrections=(_coulomb_correction(electrostatics_kcal),),
        extra={
            "raw_terms": dict(raw_terms),
            "raw_units": "kcal/mol",
            "scope": "finite-chain GLYCAM",
        },
    )
