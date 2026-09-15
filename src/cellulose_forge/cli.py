"""Command-line interface for CelluloseForge."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from cellulose_forge import __version__
from cellulose_forge.crystallography import MillerIndex, load_allomorph
from cellulose_forge.geometry import build_bulk_crystal, build_chain_parallel_slab, build_unit_cell
from cellulose_forge.pdb import write_pdb


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
    build.add_argument("--allomorph", choices=["I-beta"], required=True)
    build.add_argument(
        "--shape",
        choices=["unit-cell", "bulk", "crystal", "slab"],
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
        allomorph_key = {"I-beta": "i_beta"}[args.allomorph]
        definition = load_allomorph(allomorph_key)
        if args.shape == "unit-cell":
            if args.cells is not None:
                build_parser().error("--cells is not valid with --shape unit-cell")
            crystal = build_unit_cell(definition)
        elif args.shape in {"bulk", "crystal"}:
            if args.cells is None:
                build_parser().error("--cells NA NB NC is required for bulk crystals")
            crystal = build_bulk_crystal(definition, tuple(args.cells))
        else:
            if args.cells is not None:
                build_parser().error("--cells is not valid with --shape slab")
            if (
                args.face is None
                or args.layers is None
                or args.chains_per_layer is None
                or args.degree_polymerization is None
            ):
                build_parser().error(
                    "--face H K L, --layers N, --chains-per-layer N, and "
                    "--degree-polymerization N are required for slabs"
                )
            crystal = build_chain_parallel_slab(
                definition,
                MillerIndex(*args.face),
                layers=args.layers,
                chains_per_layer=args.chains_per_layer,
                degree_polymerization=args.degree_polymerization,
                vacuum=args.vacuum,
            )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        write_pdb(crystal, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
