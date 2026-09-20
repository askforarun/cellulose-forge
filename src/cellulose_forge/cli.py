"""Command-line interface for CelluloseForge."""

from __future__ import annotations

import argparse
import json
import tempfile
from collections.abc import Sequence
from pathlib import Path

from cellulose_forge import __version__
from cellulose_forge.coordinates import write_coordinates
from cellulose_forge.crystallography import AllomorphDefinition, MillerIndex, load_allomorph
from cellulose_forge.doctor import render_json, render_text, run_diagnostics
from cellulose_forge.geometry import (
    BulkCrystal,
    build_36_chain_fibril,
    build_bulk_crystal,
    build_chain_parallel_slab,
    build_slab_by_dimensions,
    build_unit_cell,
)
from cellulose_forge.glycam import GLYCAM_LEAPRC, generate_glycam_gromacs_reference
from cellulose_forge.glycam_mapping import (
    map_crystal_to_finite_glycam,
    map_crystal_to_periodic_glycam,
)
from cellulose_forge.gromacs import write_glycam_gromacs_system
from cellulose_forge.opls_cm5 import (
    derive_periodic_opls_cm5_topology,
    map_finite_dp20_opls_cm5,
    map_periodic_opls_cm5,
    periodic_interaction_counts,
    write_opls_cm5_system,
    write_periodic_opls_cm5_system,
)


def build_parser() -> argparse.ArgumentParser:
    """Create the top-level command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cellulose-forge",
        description="Build crystalline cellulose structures for GROMACS simulations.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    subparsers = parser.add_subparsers(dest="command")
    build = subparsers.add_parser("build", help="build a cellulose structure")
    build.add_argument(
        "--allomorph",
        choices=["I-alpha", "I-beta", "II", "III-I"],
        required=True,
    )
    build.add_argument(
        "--shape",
        choices=["unit-cell", "bulk", "crystal", "slab", "fibril"],
        required=True,
    )
    build.add_argument(
        "--cells",
        nargs=3,
        type=int,
        metavar=("NA", "NB", "NC"),
        help="unit-cell repetitions along crystallographic a, b, and c",
    )
    build.add_argument(
        "--face",
        nargs=3,
        type=int,
        metavar=("H", "K", "L"),
        help="signed Miller index of the exposed slab face",
    )
    build.add_argument("--layers", type=int, help="number of complete molecular slab layers")
    build.add_argument("--chains-per-layer", type=int, help="cellulose chains in each slab layer")
    build.add_argument(
        "--degree-polymerization",
        type=int,
        help="anhydroglucose units (AGUs) in each cellulose chain",
    )
    build.add_argument(
        "--fibril-model",
        choices=["36-chain"],
        help="named literature-defined fibril cross-section",
    )
    build.add_argument(
        "--chain-boundary",
        choices=["finite", "periodic"],
        help="finite termini or a covalent beta(1->4) bond across the chain boundary",
    )
    build.add_argument("--slab-width", type=float, help="minimum slab width in angstrom")
    build.add_argument("--slab-thickness", type=float, help="minimum slab thickness in angstrom")
    build.add_argument(
        "--slab-length", type=float, help="minimum chain-direction length in angstrom"
    )
    build.add_argument(
        "--vacuum",
        type=float,
        default=0.0,
        help="total slab-normal or fibril-transverse padding in angstrom (default: 0)",
    )
    build.add_argument(
        "--padding-x",
        type=float,
        help="total fibril padding along Cartesian x; overrides --vacuum for x",
    )
    build.add_argument(
        "--padding-y",
        type=float,
        help="total fibril padding along Cartesian y; overrides --vacuum for y",
    )
    build.add_argument(
        "--padding-z",
        type=float,
        help="total chain-direction padding for finite fibrils; periodic fibrils require zero",
    )
    build.add_argument(
        "--force-field",
        choices=["none", "glycam", "opls-cm5"],
        default="none",
        help="write coordinates only or a simulation-ready GROMACS system",
    )
    build.add_argument(
        "--opls-template",
        type=Path,
        help="published pristine DP20 OPLS-CM5 CNC.itp reference",
    )
    build.add_argument(
        "--opls-reference-gro",
        type=Path,
        help="published pristine CNC GRO containing the matching DP20 chain",
    )
    build.add_argument(
        "--tleap-executable",
        default="tleap",
        help="tleap executable used for GLYCAM reference generation",
    )
    build.add_argument(
        "--parmed-python",
        type=Path,
        help="optional Python interpreter containing ParmEd",
    )
    build.add_argument("--output", type=Path, required=True)

    doctor = subparsers.add_parser(
        "doctor",
        help="check the environment for coordinate-only, GLYCAM, and GROMACS readiness",
    )
    doctor.add_argument(
        "--tleap-executable",
        default="tleap",
        help="tleap executable to check for GLYCAM reference generation",
    )
    doctor.add_argument(
        "--parmed-python",
        type=Path,
        help="Python interpreter to check for an importable ParmEd",
    )
    doctor.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    return parser


def main(argv: Sequence[str] | None = None) -> int:  # noqa: PLR0912
    """Run the command-line interface."""
    args = build_parser().parse_args(argv)
    if args.command == "doctor":
        checks = run_diagnostics(
            tleap_executable=args.tleap_executable,
            parmed_python=str(args.parmed_python) if args.parmed_python else None,
        )
        print(render_json(checks) if args.json else render_text(checks))
        return 0 if all(check.ok for check in checks if check.name != "gmx") else 1
    if args.command == "build":
        allomorph_key = {
            "I-alpha": "i_alpha",
            "I-beta": "i_beta",
            "II": "ii",
            "III-I": "iii_i",
        }[args.allomorph]
        definition = load_allomorph(allomorph_key)
        if args.shape != "fibril" and any(
            value is not None for value in (args.padding_x, args.padding_y, args.padding_z)
        ):
            build_parser().error(
                "--padding-x, --padding-y, and --padding-z are valid only with --shape fibril"
            )
        if args.shape == "unit-cell":
            if args.cells is not None:
                build_parser().error("--cells is not valid with --shape unit-cell")
            crystal = build_unit_cell(definition)
        elif args.shape in {"bulk", "crystal"}:
            if args.cells is None:
                build_parser().error("--cells NA NB NC is required for bulk crystals")
            crystal = build_bulk_crystal(definition, tuple(args.cells))
        elif args.shape == "slab":
            crystal = _build_slab(args, definition, build_parser())
        else:
            crystal = _build_fibril(args, definition, build_parser())
        if args.force_field == "glycam":
            _write_glycam_build(args, crystal, build_parser())
        elif args.force_field == "opls-cm5":
            _write_opls_cm5_build(args, crystal, build_parser())
        else:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            written_path = write_coordinates(crystal, args.output)
            if written_path != args.output:
                print(
                    f"Requested structure has {crystal.atom_count} atoms; "
                    f"wrote mmCIF instead of PDB: {written_path}"
                )
    return 0


def _write_opls_cm5_build(
    args: argparse.Namespace,
    crystal: BulkCrystal,
    parser: argparse.ArgumentParser,
) -> None:
    """Map pristine I-beta to a finite or derived periodic OPLS-CM5 topology."""
    if args.output.suffix:
        parser.error("--output must be a directory when --force-field opls-cm5 is selected")
    if args.allomorph != "I-beta":
        parser.error("the experimental OPLS-CM5 backend currently supports I-beta only")
    if args.opls_template is None or args.opls_reference_gro is None:
        parser.error(
            "--opls-template and --opls-reference-gro are required with "
            "--force-field opls-cm5"
        )
    chain_boundary = args.chain_boundary
    if args.shape in {"bulk", "crystal", "unit-cell"}:
        if chain_boundary != "periodic":
            parser.error("OPLS-CM5 bulk and unit-cell builds require --chain-boundary periodic")
    elif chain_boundary is None:
        parser.error("--chain-boundary finite|periodic is required with OPLS-CM5")
    assert chain_boundary is not None
    degree_polymerization = _crystal_degree_polymerization(crystal)
    periodic_topology = None
    if chain_boundary == "periodic":
        periodic_topology = derive_periodic_opls_cm5_topology(
            args.opls_template,
            args.opls_reference_gro,
            degree_polymerization=degree_polymerization,
        )
        mapped = map_periodic_opls_cm5(crystal, periodic_topology)
    else:
        if args.shape not in {"slab", "fibril"}:
            parser.error("finite OPLS-CM5 currently supports slabs and fibrils only")
        mapped = map_finite_dp20_opls_cm5(
            crystal,
            args.opls_template,
            args.opls_reference_gro,
        )
    args.output.mkdir(parents=True, exist_ok=True)
    coordinate_path = write_coordinates(crystal, args.output / "structure.pdb")
    if periodic_topology is None:
        gro_path, itp_path, topology_path = write_opls_cm5_system(
            mapped, args.opls_template, args.output
        )
    else:
        gro_path, itp_path, topology_path = write_periodic_opls_cm5_system(
            mapped, periodic_topology, args.output
        )
    report = {
        "cellulose_forge_version": __version__,
        "allomorph": crystal.allomorph,
        "shape": args.shape,
        "force_field": "OPLS-CM5",
        "status": "experimental",
        "chain_boundary": chain_boundary,
        "degree_polymerization": degree_polymerization,
        "chain_count": mapped.chain_count,
        "crystallographic_atom_count": crystal.atom_count,
        "force_field_atom_count": mapped.atom_count,
        "reference_topology": str(args.opls_template),
        "reference_coordinates": str(args.opls_reference_gro),
        "outputs": {
            "coordinates": coordinate_path.name,
            "gromacs_coordinates": gro_path.name,
            "gromacs_molecule_topology": itp_path.name,
            "gromacs_topology": topology_path.name,
        },
        "geometry_annotations": list(crystal.annotations),
    }
    if periodic_topology is not None:
        report["periodic_topology_derivation"] = {
            "source": "internal repeat of published finite DP20 topology",
            "source_repeat_charge": periodic_topology.source_repeat_charge,
            "uniform_charge_correction_per_atom": (
                periodic_topology.charge_correction_per_atom
            ),
            "corrected_repeat_charge": periodic_topology.corrected_repeat_charge,
            "interaction_counts_per_chain": periodic_interaction_counts(
                periodic_topology
            ),
        }
    (args.output / "build_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_glycam_build(
    args: argparse.Namespace,
    crystal: BulkCrystal,
    parser: argparse.ArgumentParser,
) -> None:
    """Generate a GLYCAM reference, map coordinates, and write GROMACS files."""
    if args.output.suffix:
        parser.error("--output must be a directory when --force-field glycam is selected")
    chain_boundary = args.chain_boundary
    if args.shape in {"bulk", "crystal", "unit-cell"}:
        if chain_boundary is None:
            parser.error(
                "--chain-boundary periodic is required for GLYCAM bulk and unit-cell builds"
            )
        if chain_boundary != "periodic":
            parser.error("GLYCAM bulk and unit-cell builds require --chain-boundary periodic")
    elif chain_boundary is None:
        parser.error("--chain-boundary finite|periodic is required with --force-field glycam")
    assert chain_boundary is not None
    degree_polymerization = _crystal_degree_polymerization(crystal)
    output_directory = args.output

    with tempfile.TemporaryDirectory(prefix="cellulose_forge_glycam_") as temporary:
        files, reference_metadata = generate_glycam_gromacs_reference(
            degree_polymerization=degree_polymerization,
            chain_boundary=chain_boundary,
            directory=temporary,
            tleap_executable=args.tleap_executable,
            parmed_python=args.parmed_python,
        )
        if chain_boundary == "periodic":
            mapped = map_crystal_to_periodic_glycam(crystal, files.pdb)
        else:
            mapped = map_crystal_to_finite_glycam(crystal, files.pdb)
        output_directory.mkdir(parents=True, exist_ok=True)
        coordinate_path = write_coordinates(crystal, output_directory / "structure.pdb")
        gro_path, topology_path = write_glycam_gromacs_system(
            mapped,
            files.gromacs_topology,
            output_directory,
        )

    report = {
        "cellulose_forge_version": __version__,
        "allomorph": crystal.allomorph,
        "shape": args.shape,
        "force_field": "GLYCAM",
        "force_field_release": GLYCAM_LEAPRC,
        "chain_boundary": chain_boundary,
        "degree_polymerization": degree_polymerization,
        "chain_count": mapped.chain_count,
        "crystallographic_atom_count": crystal.atom_count,
        "force_field_atom_count": mapped.atom_count,
        "outputs": {
            "coordinates": coordinate_path.name,
            "gromacs_coordinates": gro_path.name,
            "gromacs_topology": topology_path.name,
        },
        "geometry_annotations": list(crystal.annotations),
        "reference_generation": reference_metadata,
    }
    (output_directory / "build_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _crystal_degree_polymerization(crystal: BulkCrystal) -> int:
    """Return the common residue count for every generated cellulose chain."""
    degrees = {
        chain_number: len(
            {
                residue_number
                for atom_chain, residue_number in zip(
                    crystal.chain_numbers,
                    crystal.residue_numbers,
                    strict=True,
                )
                if atom_chain == chain_number
            }
        )
        for chain_number in dict.fromkeys(crystal.chain_numbers)
    }
    if not degrees or len(set(degrees.values())) != 1:
        msg = "GLYCAM output requires equal, nonempty cellulose chain lengths"
        raise ValueError(msg)
    return next(iter(degrees.values()))


def _build_fibril(
    args: argparse.Namespace,
    definition: AllomorphDefinition,
    parser: argparse.ArgumentParser,
) -> BulkCrystal:
    """Validate named-fibril options and build its crystallographic geometry."""
    if args.cells is not None or args.face is not None:
        parser.error("--cells and --face are not valid with --shape fibril")
    if any(
        value is not None
        for value in (
            args.layers,
            args.chains_per_layer,
            args.slab_width,
            args.slab_thickness,
            args.slab_length,
        )
    ):
        parser.error("slab sizing options are not valid with --shape fibril")
    if args.fibril_model is None:
        parser.error("--fibril-model 36-chain is required for fibrils")
    if args.degree_polymerization is None:
        parser.error("--degree-polymerization N is required for fibrils")
    if args.chain_boundary is None:
        parser.error("--chain-boundary finite|periodic is required for fibrils")
    if args.chain_boundary == "periodic" and args.padding_z not in {None, 0}:
        parser.error("--padding-z must be zero with --chain-boundary periodic")
    return build_36_chain_fibril(
        definition,
        degree_polymerization=args.degree_polymerization,
        vacuum=args.vacuum,
        padding_x=args.padding_x,
        padding_y=args.padding_y,
        padding_z=args.padding_z,
        chain_boundary=args.chain_boundary,
    )


def _build_slab(
    args: argparse.Namespace,
    definition: AllomorphDefinition,
    parser: argparse.ArgumentParser,
) -> BulkCrystal:
    """Resolve mutually exclusive slab sizing modes and build the slab."""
    if args.cells is not None:
        parser.error("--cells is not valid with --shape slab")
    if args.face is None:
        parser.error("--face H K L is required for slabs")
    if args.chain_boundary is None:
        parser.error("--chain-boundary finite|periodic is required for slabs")
    count_values = (args.layers, args.chains_per_layer, args.degree_polymerization)
    size_values = (args.slab_width, args.slab_thickness, args.slab_length)
    has_count_value = any(value is not None for value in count_values)
    has_size_value = any(value is not None for value in size_values)
    if has_count_value and has_size_value:
        parser.error("exact-count and physical-size slab options cannot be combined")
    if has_count_value:
        if not all(value is not None for value in count_values):
            parser.error(
                "--layers N, --chains-per-layer N, and --degree-polymerization N "
                "must be supplied together"
            )
        assert args.layers is not None
        assert args.chains_per_layer is not None
        assert args.degree_polymerization is not None
        return build_chain_parallel_slab(
            definition,
            MillerIndex(*args.face),
            layers=args.layers,
            chains_per_layer=args.chains_per_layer,
            degree_polymerization=args.degree_polymerization,
            vacuum=args.vacuum,
            chain_boundary=args.chain_boundary,
        )
    if has_size_value:
        if not all(value is not None for value in size_values):
            parser.error(
                "--slab-width A, --slab-thickness A, and --slab-length A must be supplied together"
            )
        assert args.slab_width is not None
        assert args.slab_thickness is not None
        assert args.slab_length is not None
        return build_slab_by_dimensions(
            definition,
            MillerIndex(*args.face),
            width=args.slab_width,
            thickness=args.slab_thickness,
            length=args.slab_length,
            vacuum=args.vacuum,
            chain_boundary=args.chain_boundary,
        )
    parser.error(
        "supply either exact counts (--layers, --chains-per-layer, "
        "--degree-polymerization) or physical dimensions (--slab-width, "
        "--slab-thickness, --slab-length)"
    )


if __name__ == "__main__":
    raise SystemExit(main())
