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


_CELLULOSE_FORCE_FIELD_SHAPES = frozenset({"bulk", "crystal", "slab", "fibril"})
_CHITIN_FORCE_FIELD_SHAPES = frozenset({"bulk", "crystal", "slab"})
_ENGINE_ORDER = ("gromacs", "openmm", "lammps", "amber")

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
                frozenset({"gromacs", "openmm", "lammps"}),
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
                frozenset({"gromacs", "openmm", "lammps"}),
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
    """Return an actionable error before any molecular construction begins."""
    profile = polymer_profile(polymer)
    capability = profile.force_fields.get(force_field)
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
    if force_field == "none" and engines:
        return "--md-engine requires a supported force-field build"
    if shape not in capability.shapes:
        if force_field in {"opls-aa", "opls-aa-sei"}:
            name = "OPLS-AA-SEI" if force_field == "opls-aa-sei" else "OPLS-AA"
            return f"experimental {name} supports slabs, fibrils, and periodic bulk crystals"
        supported = ", ".join(sorted(capability.shapes))
        return f"{polymer} with --force-field {force_field} supports --shape {supported}"
    unsupported_engines = set(engines) - capability.engines
    if unsupported_engines:
        if force_field in {"opls-aa", "opls-aa-sei"}:
            return (
                "experimental OPLS-AA and OPLS-AA-SEI support --md-engine "
                "gromacs, openmm, and lammps"
            )
        ordered = [engine for engine in _ENGINE_ORDER if engine in capability.engines]
        supported = (
            f"{ordered[0]} and {ordered[1]}"
            if len(ordered) == 2  # noqa: PLR2004 - natural-language list grammar
            else f"{', '.join(ordered[:-1])}, and {ordered[-1]}"
        )
        return f"{polymer} with --force-field {force_field} supports --md-engine {supported}"
    if force_field != "none" and not engines:
        return f"{force_field} builds require at least one --md-engine"
    if "amber" in engines and chain_boundary == "periodic":
        return "native Amber output is available only for finite chains"
    return None
