"""Command-line interface for CelluloseForge."""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from cellulose_forge import __version__


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
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line interface."""
    build_parser().parse_args(argv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
