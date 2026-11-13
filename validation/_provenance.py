"""Shared provenance block for validation result files."""

from __future__ import annotations

import subprocess
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def provenance(packages: tuple[str, ...]) -> dict[str, object]:
    """Record the CelluloseForge commit, tree state, and package versions."""

    def git(*arguments: str) -> str:
        result = subprocess.run(
            ["git", *arguments], cwd=ROOT, text=True, capture_output=True, check=False
        )
        return result.stdout.strip()

    versions: dict[str, str] = {}
    for package in packages:
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = "not installed"
    return {
        "cellulose_forge_commit": git("rev-parse", "HEAD"),
        "tracked_tree_clean": git("status", "--porcelain", "--untracked-files=no") == "",
        "software_versions": versions,
    }
