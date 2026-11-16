"""Normalize classified OpenMM single-point components.

OpenMM force classes do not map one-to-one onto the canonical report.  In
particular, CHARMM Urey--Bradley terms are stored in a bond force, periodic
impropers can share a torsion force with propers, and dispersion corrections
are included inside nonbonded forces.  The evaluator therefore classifies and
isolates those contributions before calling this normalization layer.
"""

from __future__ import annotations

from collections.abc import Mapping

from cellulose_forge.energy_drivers.report import EnergyReport, EnergyValue

OPENMM_COMPONENTS = frozenset(
    {
        "bond",
        "angle",
        "urey_bradley",
        "proper_dihedral",
        "periodic_improper",
        "harmonic_improper",
        "vdw",
        "vdw_14",
        "electrostatics",
        "electrostatics_14",
        "restraint",
        "dispersion_correction",
        "potential",
    }
)


def openmm_energy_report(
    components: Mapping[str, EnergyValue],
    *,
    coordinate_sha256: str,
    extra: Mapping[str, object] | None = None,
) -> EnergyReport:
    """Return a canonical report from already classified OpenMM components.

    All nonbonded 1--4 contributions are folded into their corresponding
    canonical aggregate.  Urey--Bradley contributes to ``angle`` and both
    periodic and harmonic impropers contribute to ``improper``.
    """
    unknown = set(components) - OPENMM_COMPONENTS
    if unknown:
        raise ValueError(f"OpenMM components contain unsupported terms: {sorted(unknown)}")
    if "potential" not in components:
        raise ValueError("OpenMM components are missing potential")

    def total(*names: str) -> EnergyValue:
        return sum((components.get(name, 0) for name in names), start=0)

    terms: dict[str, EnergyValue] = {
        "bond": total("bond"),
        "angle": total("angle", "urey_bradley"),
        "dihedral": total("proper_dihedral"),
        "improper": total("periodic_improper", "harmonic_improper"),
        "vdw": total("vdw", "vdw_14"),
        "electrostatics": total("electrostatics", "electrostatics_14"),
        "restraint": total("restraint"),
        "dispersion_correction": total("dispersion_correction"),
        "potential": components["potential"],
    }
    return EnergyReport(
        engine="openmm",
        coordinate_sha256=coordinate_sha256,
        terms=terms,
        extra={} if extra is None else extra,
    )
