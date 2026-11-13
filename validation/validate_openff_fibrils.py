"""Validate experimental OpenFF Sage cellulose fibrils.

Covers the literature 36-chain fibril for I-alpha, I-beta and II with finite
and covalently periodic chains, plus the I-beta 18- and 24-chain presets and a
custom row profile.  Like the slab and bulk matrices, this is a
model-construction validation: graph mapping, stereochemistry, charges, the
vacuum dispersion policy, detailed OpenMM/GROMACS energies, and short
minimization/NVT smoke tests in both engines.  It does not assess Sage's
physical accuracy for cellulose fibrils.
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from importlib import metadata
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "validation"))

from _provenance import provenance  # noqa: E402
from validate_openff_finite_slab import (  # noqa: E402
    bonded_term_counts,
    compare_energy_terms,
)
from validate_openff_periodic_bulk import (  # noqa: E402
    _gromacs_energies,
    _interchange,
    _openmm_energies,
    _run_gromacs,
    _run_openmm,
)

from cellulose_forge.crystallography import load_allomorph  # noqa: E402
from cellulose_forge.geometry import (  # noqa: E402
    I_BETA_FIBRIL_PROFILES,
    BulkCrystal,
    build_36_chain_fibril,
    build_i_beta_profile_fibril,
)
from cellulose_forge.openff_sage import (  # noqa: E402
    build_sage_system,
    dispersion_correction_for_shape,
    periodic_box_heights,
)

DEGREE_POLYMERIZATION = 6
VACUUM_ANGSTROM = 20.0
FINITE_PADDING_ANGSTROM = 12.0
MINIMUM_BOX_HEIGHT_ANGSTROM = 18.0
CUSTOM_PROFILE = (2, 3, 3, 2)

# name -> (allomorph, model, row profile or None, chain boundary)
CASES: dict[str, tuple[str, str, tuple[int, ...] | None, str]] = {
    "i_alpha_36_finite": ("i_alpha", "36-chain", None, "finite"),
    "i_alpha_36_periodic": ("i_alpha", "36-chain", None, "periodic"),
    "i_beta_36_finite": ("i_beta", "36-chain", None, "finite"),
    "i_beta_36_periodic": ("i_beta", "36-chain", None, "periodic"),
    "ii_36_finite": ("ii", "36-chain", None, "finite"),
    "ii_36_periodic": ("ii", "36-chain", None, "periodic"),
    "i_beta_18_finite": ("i_beta", "18-chain", I_BETA_FIBRIL_PROFILES["18-chain"], "finite"),
    "i_beta_24_finite": ("i_beta", "24-chain", I_BETA_FIBRIL_PROFILES["24-chain"], "finite"),
    "i_beta_24_periodic": ("i_beta", "24-chain", I_BETA_FIBRIL_PROFILES["24-chain"], "periodic"),
    "i_beta_custom_finite": ("i_beta", "custom", CUSTOM_PROFILE, "finite"),
}


def _crystal(
    allomorph: str, model: str, profile: tuple[int, ...] | None, boundary: str
) -> BulkCrystal:
    padding = FINITE_PADDING_ANGSTROM if boundary == "finite" else None
    definition = load_allomorph(allomorph)
    if profile is None:
        return build_36_chain_fibril(
            definition,
            degree_polymerization=DEGREE_POLYMERIZATION,
            vacuum=VACUUM_ANGSTROM,
            padding_z=padding,
            chain_boundary=boundary,
        )
    return build_i_beta_profile_fibril(
        definition,
        row_counts=profile,
        model_name=model,
        degree_polymerization=DEGREE_POLYMERIZATION,
        vacuum=VACUUM_ANGSTROM,
        padding_z=padding,
        chain_boundary=boundary,
    )


def validate_case(name: str) -> dict[str, object]:
    allomorph, model, profile, boundary = CASES[name]
    periodic = boundary == "periodic"
    dispersion = dispersion_correction_for_shape("fibril")
    crystal = _crystal(allomorph, model, profile, boundary)
    heights = periodic_box_heights(crystal.basis_vectors)
    if min(heights) <= MINIMUM_BOX_HEIGHT_ANGSTROM:
        raise RuntimeError(f"{name}: box does not exceed twice the 9 A cutoff")
    build = build_sage_system(crystal, boundary)
    interchange = _interchange(build)
    with tempfile.TemporaryDirectory(prefix=f"openff_fibril_{name}_") as temporary:
        directory = Path(temporary)
        openmm_values = _openmm_energies(interchange, periodic=periodic)
        gromacs_values = _gromacs_energies(interchange, directory, periodic=periodic)
        comparison = compare_energy_terms(
            openmm_values,
            gromacs_values,
            bonded_term_counts(
                interchange.to_openmm_system(combine_nonbonded_forces=False)  # type: ignore[attr-defined]
            ),
        )
        gromacs_smoke = _run_gromacs(
            interchange, directory, periodic=periodic, dispersion=dispersion
        )
        openmm_smoke = _run_openmm(build, interchange, periodic=periodic, dispersion=dispersion)
    return {
        "allomorph": allomorph,
        "fibril_model": model,
        "row_profile": list(profile) if profile else None,
        "chain_boundary": boundary,
        "degree_polymerization": DEGREE_POLYMERIZATION,
        "chain_count": build.mapped.chain_count,
        "atom_count": build.mapped.atom_count,
        "perpendicular_box_heights_angstrom": [float(value) for value in heights],
        "box_edges_angstrom": [float(np.linalg.norm(vector)) for vector in crystal.basis_vectors],
        "stereochemistry": build.stereochemistry,
        "charge_provenance": build.charge_provenance,
        "dispersion_correction": dispersion,
        "energy_comparison": comparison,
        "gromacs_smoke_test": gromacs_smoke,
        "openmm_smoke_test": openmm_smoke,
        "passed": bool(
            build.stereochemistry["passed"]
            and all(item["passed"] for item in comparison.values())
            and gromacs_smoke["passed"]
            and openmm_smoke["passed"]
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=sorted(CASES), action="append")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation/results/openff_fibril_nagl_matrix.json",
    )
    args = parser.parse_args()
    names = tuple(args.case or CASES)
    cases = {}
    for name in names:
        cases[name] = validate_case(name)
        print(f"{name}: passed={cases[name]['passed']}", file=sys.stderr, flush=True)
    result = {
        "matrix": "cellulose fibrils, OpenFF Sage 2.3/NAGL",
        "provenance": provenance(
            ("openff-toolkit", "openff-interchange", "openff-nagl", "openff-nagl-models", "openmm")
        ),
        "validation_scope": (
            "model construction validated; force-field accuracy for cellulose not assessed"
        ),
        "software_versions": {
            package: metadata.version(package)
            for package in (
                "openff-toolkit",
                "openff-interchange",
                "openff-nagl",
                "openff-nagl-models",
                "openmm",
            )
        },
        "openmm_platform": {"name": "CPU", "threads": 1},
        "pme_settings": {"ewald_rtol": 1e-6, "pme_order": 8, "fourierspacing_nm": 0.05},
        "dispersion_correction": "disabled: fibrils are surrounded by transverse vacuum",
        "cases": cases,
        "passed": all(case["passed"] for case in cases.values()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
