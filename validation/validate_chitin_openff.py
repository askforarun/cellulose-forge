"""Validate experimental OpenFF Sage support for anhydrous alpha/beta chitin.

This is a construction and engine-consistency validation, not an assessment of
Sage's physical accuracy for crystalline chitin.  It covers periodic bulk
crystals and finite/periodic (100) slabs, with NAGL charges, independent
stereochemical negative controls, term-resolved GROMACS/OpenMM energies, and
one-ps smoke tests in both engines.
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from dataclasses import replace
from importlib import metadata
from pathlib import Path
from typing import Any

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

from cellulose_forge.crystallography import MillerIndex, load_allomorph  # noqa: E402
from cellulose_forge.geometry import (  # noqa: E402
    build_bulk_crystal,
    build_chain_parallel_slab,
)
from cellulose_forge.openff_sage import (  # noqa: E402
    build_sage_system,
    dispersion_correction_for_shape,
    periodic_box_heights,
    validate_sage_stereochemistry,
)
from cellulose_forge.openmm import validation_platform_for_box  # noqa: E402

ALLOMORPHS = ("alpha", "beta")
SLAB_FACE = MillerIndex(1, 0, 0)
LAYERS = 2
CHAINS_PER_LAYER = 4
DEGREE_POLYMERIZATION = 6
CHARGE_TOLERANCE_E = 1.0e-6
MINIMUM_SAGE_BOX_HEIGHT_ANGSTROM = 19.0
MATRIX = tuple(
    [(allomorph, "bulk", "periodic") for allomorph in ALLOMORPHS]
    + [
        (allomorph, "slab", boundary)
        for allomorph in ALLOMORPHS
        for boundary in ("finite", "periodic")
    ]
)


def _replication(allomorph: str) -> tuple[int, int, int]:
    vectors = load_allomorph(allomorph, polymer="chitin").unit_cell.vectors()
    candidates = ((i, j, k) for i in range(1, 9) for j in range(1, 9) for k in range(1, 9))
    valid = [
        replication
        for replication in candidates
        if min(periodic_box_heights(vectors * np.asarray(replication)[:, None]))
        > MINIMUM_SAGE_BOX_HEIGHT_ANGSTROM
    ]
    if not valid:
        raise RuntimeError(f"could not size an OpenFF chitin box for {allomorph}")
    return min(valid, key=lambda item: (item[0] * item[1] * item[2], item))


def _negative_controls(build: Any) -> dict[str, bool]:
    mirrored = build.mapped.cartesian_coordinates.copy()
    mirrored[:, 0] *= -1
    mirror_rejected = False
    try:
        validate_sage_stereochemistry(
            build.molecule,
            replace(build.mapped, cartesian_coordinates=mirrored),
            "chitin",
        )
    except ValueError:
        mirror_rejected = True

    inverted = build.mapped.cartesian_coordinates.copy()
    centre = next(
        atom
        for atom in build.molecule.atoms
        if atom.name == "C2" and int(atom.metadata["residue_number"]) == 1
    )
    first, second = sorted(centre.bonded_atoms, key=lambda atom: (atom.symbol, atom.name))[:2]
    inverted[[first.molecule_atom_index, second.molecule_atom_index]] = inverted[
        [second.molecule_atom_index, first.molecule_atom_index]
    ]
    inversion_rejected = False
    try:
        validate_sage_stereochemistry(
            build.molecule,
            replace(build.mapped, cartesian_coordinates=inverted),
            "chitin",
        )
    except ValueError:
        inversion_rejected = True
    return {
        "mirrored_structure_rejected": mirror_rejected,
        "local_c2_inversion_rejected": inversion_rejected,
        "passed": mirror_rejected and inversion_rejected,
    }


def _crystal(allomorph: str, shape: str, boundary: str) -> tuple[Any, list[int] | None]:
    definition = load_allomorph(allomorph, polymer="chitin")
    if shape == "bulk":
        replication = _replication(allomorph)
        return build_bulk_crystal(definition, replication), list(replication)
    return (
        build_chain_parallel_slab(
            definition,
            SLAB_FACE,
            layers=LAYERS,
            chains_per_layer=CHAINS_PER_LAYER,
            degree_polymerization=DEGREE_POLYMERIZATION,
            vacuum=30.0,
            padding_z=12.0 if boundary == "finite" else 0.0,
            chain_boundary=boundary,
        ),
        None,
    )


def validate_case(allomorph: str, shape: str, boundary: str) -> dict[str, object]:
    """Validate one anhydrous chitin geometry in GROMACS and OpenMM."""
    crystal, replication = _crystal(allomorph, shape, boundary)
    periodic = boundary == "periodic"
    dispersion = dispersion_correction_for_shape(shape)
    build = build_sage_system(crystal, boundary, "chitin")
    interchange: Any = _interchange(build)
    with tempfile.TemporaryDirectory(
        prefix=f"chitin_openff_{allomorph}_{shape}_{boundary}_"
    ) as temporary:
        directory = Path(temporary)
        openmm_values = _openmm_energies(interchange, periodic=periodic)
        gromacs_values = _gromacs_energies(interchange, directory, periodic=periodic)
        comparison = compare_energy_terms(
            openmm_values,
            gromacs_values,
            bonded_term_counts(interchange.to_openmm_system(combine_nonbonded_forces=False)),
        )
        gromacs_smoke = _run_gromacs(
            interchange,
            directory,
            periodic=periodic,
            dispersion=dispersion,
        )
        openmm_smoke = _run_openmm(
            build,
            interchange,
            periodic=periodic,
            dispersion=dispersion,
        )
    controls = _negative_controls(build)
    passed = bool(
        build.stereochemistry["passed"]
        and controls["passed"]
        and build.charge_provenance["periodic_equivalent_atom_charges"]
        and abs(float(str(build.charge_provenance["net_charge_e"]))) <= CHARGE_TOLERANCE_E
        and all(values["passed"] for values in comparison.values())
        and gromacs_smoke["passed"]
        and openmm_smoke["passed"]
    )
    return {
        "allomorph": allomorph,
        "shape": shape,
        "chain_boundary": boundary,
        "face": [1, 0, 0] if shape == "slab" else None,
        "replication": replication,
        "atom_count": build.mapped.atom_count,
        "chain_count": build.mapped.chain_count,
        "degree_polymerization": build.degree_polymerization,
        "perpendicular_box_heights_angstrom": list(periodic_box_heights(crystal.basis_vectors)),
        "charge_provenance": build.charge_provenance,
        "stereochemistry": build.stereochemistry,
        "stereochemistry_negative_controls": controls,
        "openmm_build_validation_platform": validation_platform_for_box(crystal.basis_vectors),
        "dispersion_correction": dispersion,
        "energy_comparison": comparison,
        "gromacs_smoke_test": gromacs_smoke,
        "openmm_smoke_test": openmm_smoke,
        "passed": passed,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation" / "results" / "chitin_openff_matrix.json",
    )
    parser.add_argument(
        "--case",
        choices=("all", *("_".join(case) for case in MATRIX)),
        default="all",
        help="run the full matrix or one named case",
    )
    args = parser.parse_args()
    matrix = MATRIX if args.case == "all" else (tuple(args.case.split("_")),)
    cases = {"_".join(case): validate_case(*case) for case in matrix}
    packages = (
        "openff-toolkit",
        "openff-interchange",
        "openff-nagl",
        "openff-nagl-models",
        "openmm",
    )
    result = {
        "matrix": "anhydrous alpha/beta chitin, OpenFF Sage 2.3/NAGL",
        "provenance": provenance(packages),
        "software_versions": {package: metadata.version(package) for package in packages},
        "validation_scope": (
            "model construction and cross-engine consistency validated; "
            "force-field accuracy for crystalline chitin not assessed"
        ),
        "pme_settings": {
            "ewald_rtol": 1.0e-6,
            "pme_order": 8,
            "fourierspacing_nm": 0.05,
        },
        "cases": cases,
        "passed": all(case["passed"] for case in cases.values()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
