from __future__ import annotations

import pytest

from cellulose_forge.energy_drivers.openmm import openmm_energy_report
from cellulose_forge.energy_drivers.report import validate_potential_sum


def test_openmm_components_map_to_canonical_terms() -> None:
    report = openmm_energy_report(
        {
            "bond": 1.0,
            "angle": 2.0,
            "urey_bradley": 3.0,
            "proper_dihedral": 4.0,
            "periodic_improper": 5.0,
            "harmonic_improper": 6.0,
            "vdw": 7.0,
            "vdw_14": 8.0,
            "electrostatics": 9.0,
            "electrostatics_14": 10.0,
            "restraint": 11.0,
            "dispersion_correction": 12.0,
            "potential": 78.0,
        },
        coordinate_sha256="abc",
    )

    assert report.terms == {
        "bond": 1.0,
        "angle": 5.0,
        "dihedral": 4.0,
        "improper": 11.0,
        "vdw": 15.0,
        "electrostatics": 19.0,
        "restraint": 11.0,
        "dispersion_correction": 12.0,
        "potential": 78.0,
    }
    validate_potential_sum(report, absolute_kj_mol=0.0, relative=0.0)


def test_openmm_missing_components_are_explicit_zeros() -> None:
    report = openmm_energy_report({"potential": 0}, coordinate_sha256="abc")

    assert set(report.terms.values()) == {0}


@pytest.mark.parametrize("components", [{"bond": 1.0}, {"potential": 0, "cmap": 1.0}])
def test_openmm_refuses_incomplete_or_unclassified_components(
    components: dict[str, float],
) -> None:
    with pytest.raises(ValueError):
        openmm_energy_report(components, coordinate_sha256="abc")


def test_openmm_extra_metadata_is_copied() -> None:
    metadata: dict[str, object] = {"platform": "CPU"}
    report = openmm_energy_report({"potential": 0}, coordinate_sha256="abc", extra=metadata)
    metadata["platform"] = "Reference"

    assert report.extra["platform"] == "CPU"
