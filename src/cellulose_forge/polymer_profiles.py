"""Declarative polymer identities and supported build combinations."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass


@dataclass(frozen=True)
class ForceFieldCapability:
    """Shapes and MD engines supported by one polymer/force-field pair."""

    shapes: frozenset[str]
    engines: frozenset[str]


@dataclass(frozen=True)
class PolymerProfile:
    """Names, validation reference, and build capabilities for one polymer."""

    name: str
    residue_name: str
    validation_reference: str
    force_fields: Mapping[str, ForceFieldCapability]


_CELLULOSE_FORCE_FIELD_SHAPES = frozenset({"unit-cell", "bulk", "crystal", "slab", "fibril"})
_CHITIN_FORCE_FIELD_SHAPES = frozenset({"unit-cell", "bulk", "crystal", "slab"})
_ENGINE_ORDER = ("gromacs", "openmm", "lammps", "namd", "amber")
_SHAPE_ORDER = ("unit-cell", "bulk", "crystal", "slab", "fibril")

POLYMER_PROFILES: Mapping[str, PolymerProfile] = {
    "cellulose": PolymerProfile(
        name="cellulose",
        residue_name="BGC",
        validation_reference="independent beta-D-glucopyranose signed parity",
        force_fields={
            "none": ForceFieldCapability(
                frozenset({"unit-cell", "bulk", "crystal", "slab", "fibril"}),
                frozenset(),
            ),
            "glycam": ForceFieldCapability(
                _CELLULOSE_FORCE_FIELD_SHAPES,
                frozenset({"gromacs", "amber", "openmm", "lammps"}),
            ),
            "charmm36": ForceFieldCapability(
                _CELLULOSE_FORCE_FIELD_SHAPES,
                frozenset({"gromacs", "openmm", "lammps", "namd"}),
            ),
            "openff": ForceFieldCapability(
                _CELLULOSE_FORCE_FIELD_SHAPES,
                frozenset({"gromacs", "openmm"}),
            ),
            "opls-aa": ForceFieldCapability(
                _CELLULOSE_FORCE_FIELD_SHAPES,
                frozenset({"gromacs", "openmm", "lammps"}),
            ),
            "opls-aa-sei": ForceFieldCapability(
                _CELLULOSE_FORCE_FIELD_SHAPES,
                frozenset({"gromacs", "openmm", "lammps"}),
            ),
        },
    ),
    "chitin": PolymerProfile(
        name="chitin",
        residue_name="NAG",
        validation_reference="experimental anhydrous alpha/beta-chitin crystal signed parity",
        force_fields={
            "none": ForceFieldCapability(
                frozenset({"unit-cell", "bulk", "crystal", "slab"}), frozenset()
            ),
            "glycam": ForceFieldCapability(
                _CHITIN_FORCE_FIELD_SHAPES,
                frozenset({"gromacs", "amber", "openmm", "lammps"}),
            ),
            "charmm36": ForceFieldCapability(
                _CHITIN_FORCE_FIELD_SHAPES,
                frozenset({"gromacs", "openmm", "lammps", "namd"}),
            ),
            "openff": ForceFieldCapability(
                _CHITIN_FORCE_FIELD_SHAPES,
                frozenset({"gromacs", "openmm"}),
            ),
        },
    ),
    "mannan": PolymerProfile(
        name="mannan",
        residue_name="MAN",
        validation_reference="beta-D-mannose signed parity from the crystal definition",
        force_fields={
            "none": ForceFieldCapability(frozenset({"unit-cell", "bulk", "crystal"}), frozenset())
        },
    ),
}


def polymer_profile(polymer: str) -> PolymerProfile:
    """Return the declared profile for a supported polymer."""
    try:
        return POLYMER_PROFILES[polymer]
    except KeyError as error:
        raise ValueError(f"unsupported polymer {polymer!r}") from error


def capability_error(  # noqa: PLR0911 - each capability refusal stays explicit
    polymer: str,
    force_field: str,
    shape: str,
    engines: tuple[str, ...],
    chain_boundary: str | None,
) -> str | None:
    """Return an actionable shape or boundary error before construction."""
    profile = polymer_profile(polymer)
    capability = profile.force_fields.get(force_field)
    if polymer == "mannan" and shape not in profile.force_fields["none"].shapes:
        return "mannan currently supports --shape unit-cell and bulk"
    if capability is None:
        if polymer == "chitin" and force_field in {"opls-aa", "opls-aa-sei"}:
            name = "OPLS-AA-SEI" if force_field == "opls-aa-sei" else "OPLS-AA"
            return (
                f"experimental {name} chitin output is unavailable because no "
                "published OPLS GlcNAc parameters have been adopted"
            )
        if polymer == "mannan":
            return "mannan force-field output is not implemented; use --force-field none"
        supported = ", ".join(profile.force_fields)
        return f"{polymer} supports --force-field {supported}"
    if shape not in capability.shapes:
        if polymer == "chitin" and shape == "fibril" and force_field == "openff":
            return "experimental OpenFF Sage chitin fibrils are not implemented"
        if polymer == "chitin":
            return "chitin currently supports --shape unit-cell, bulk, and slab"
        if force_field in {"opls-aa", "opls-aa-sei"}:
            name = "OPLS-AA-SEI" if force_field == "opls-aa-sei" else "OPLS-AA"
            return f"experimental {name} supports slabs, fibrils, and periodic bulk crystals"
        supported = ", ".join(item for item in _SHAPE_ORDER if item in capability.shapes)
        return f"{polymer} with --force-field {force_field} supports --shape {supported}"
    if force_field == "openff" and shape in {"bulk", "crystal"} and chain_boundary != "periodic":
        return "experimental OpenFF Sage bulk crystals require periodic chains"
    if force_field in {"opls-aa", "opls-aa-sei"}:
        name = "OPLS-AA-SEI" if force_field == "opls-aa-sei" else "OPLS-AA"
        if shape in {"bulk", "crystal"} and chain_boundary != "periodic":
            return f"experimental {name} bulk crystals require periodic chains"
    return None


def engine_capability_error(  # noqa: PLR0911 - preserve specific historical diagnostics
    polymer: str, force_field: str, engines: tuple[str, ...]
) -> str | None:
    """Return the historical, actionable MD-engine refusal for a request."""
    profile = polymer_profile(polymer)
    capability = profile.force_fields.get(force_field)
    if force_field == "none" and engines:
        return "--md-engine requires a supported force-field build"
    if force_field == "glycam" and not engines:
        return "GLYCAM builds require at least one --md-engine"
    if force_field in {"opls-aa", "opls-aa-sei"} and not engines:
        return (
            "experimental OPLS-AA and OPLS-AA-SEI support --md-engine gromacs, openmm, and lammps"
        )
    declared_engines = {
        "glycam": frozenset({"gromacs", "amber", "openmm", "lammps"}),
        "charmm36": frozenset({"gromacs", "openmm", "lammps", "namd"}),
        "openff": frozenset({"gromacs", "openmm"}),
        "opls-aa": frozenset({"gromacs", "openmm", "lammps"}),
        "opls-aa-sei": frozenset({"gromacs", "openmm", "lammps"}),
    }.get(force_field, capability.engines if capability is not None else frozenset())
    unsupported = [engine for engine in engines if engine not in declared_engines]
    if not unsupported:
        return None
    passed = ", ".join(unsupported)
    if force_field == "charmm36":
        return (
            f"CHARMM36 does not support --md-engine {passed}; use gromacs, openmm, lammps, or namd"
        )
    if force_field == "glycam" and polymer == "chitin":
        return (
            f"chitin GLYCAM does not support --md-engine {passed}; use gromacs, openmm, "
            "lammps, or amber for finite chains"
        )
    if force_field == "openff":
        return (
            f"experimental OpenFF does not support --md-engine {passed}; use gromacs or openmm only"
        )
    if force_field in {"opls-aa", "opls-aa-sei"}:
        return (
            f"experimental OPLS-AA and OPLS-AA-SEI do not support --md-engine {passed}; "
            "they support --md-engine gromacs, openmm, and lammps"
        )
    ordered = [engine for engine in _ENGINE_ORDER if engine in declared_engines]
    if len(ordered) == 1:
        supported = ordered[0]
    elif len(ordered) == 2:  # noqa: PLR2004 - natural-language list grammar
        supported = f"{ordered[0]} and {ordered[1]}"
    else:
        supported = f"{', '.join(ordered[:-1])}, and {ordered[-1]}"
    return (
        f"{polymer} with --force-field {force_field} does not support --md-engine {passed}; "
        f"use {supported}"
    )
