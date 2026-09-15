"""Command-line interface for CelluloseForge."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from cellulose_forge import __version__
from cellulose_forge.crystallography import load_allomorph
from cellulose_forge.geometry import build_bulk_crystal, build_unit_cell
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
        choices=["unit-cell", "bulk", "crystal"],
        required=True,
    )
    build.add_argument(
        "--cells",
        nargs=3,
        type=int,
        metavar=("NA", "NB", "NC"),
        help="unit-cell repetitions along crystallographic a, b, and c",
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
        else:
            if args.cells is None:
                build_parser().error("--cells NA NB NC is required for bulk crystals")
            crystal = build_bulk_crystal(definition, tuple(args.cells))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        write_pdb(crystal, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
