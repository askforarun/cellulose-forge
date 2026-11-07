"""Crystalline cellulose structure builder."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("cellulose-forge")
except PackageNotFoundError:
    __version__ = "1.1.2"

__all__ = ["__version__"]
