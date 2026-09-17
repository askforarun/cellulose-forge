"""Command-line interface for CelluloseForge."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from cellulose_forge import __version__
from cellulose_forge.coordinates import write_coordinates
from cellulose_forge.crystallography import AllomorphDefinition, MillerIndex, load_allomorph
from cellulose_forge.geometry import (
    BulkCrystal,
    build_36_chain_fibril,
    build_bulk_crystal,
    build_chain_parallel_slab,
    build_slab_by_dimensions,
    build_unit_cell,
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
        help="total vacuum normal to the slab in angstrom (default: 0)",
    )
    build.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line interface."""
    args = build_parser().parse_args(argv)
    if args.command == "build":
        allomorph_key = {
            "I-alpha": "i_alpha",
            "I-beta": "i_beta",
            "II": "ii",
            "III-I": "iii_i",
        }[args.allomorph]
        definition = load_allomorph(allomorph_key)
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
        args.output.parent.mkdir(parents=True, exist_ok=True)
        written_path = write_coordinates(crystal, args.output)
        if written_path != args.output:
            print(
                f"Requested structure has {crystal.atom_count} atoms; "
                f"wrote mmCIF instead of PDB: {written_path}"
            )
    return 0


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
    return build_36_chain_fibril(
        definition,
        degree_polymerization=args.degree_polymerization,
        vacuum=args.vacuum,
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
