"""Tests for force-field installation and discovery."""

from __future__ import annotations

import hashlib
import tarfile
from pathlib import Path

import pytest

from cellulose_forge import forcefields


def _write_archive(path: Path, source: Path) -> None:
    release = source / forcefields.CHARMM36_RELEASE
    release.mkdir(parents=True)
    for name in (
        "forcefield.itp",
        "forcefield.doc",
        "ffbonded.itp",
        "ffnonbonded.itp",
        "nbfix.itp",
        "carb.rtp",
        "carb.hdb",
    ):
        (release / name).write_text(f"test {name}\n", encoding="utf-8")
    with tarfile.open(path, "w:gz") as archive:
        archive.add(release, arcname=forcefields.CHARMM36_RELEASE)


def test_installs_and_discovers_local_archive(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    archive = tmp_path / "source.tgz"
    _write_archive(archive, tmp_path / "content")
    checksum = hashlib.sha256(archive.read_bytes()).hexdigest()
    monkeypatch.setenv("CELLULOSE_FORGE_DATA_DIR", str(tmp_path / "data"))
    monkeypatch.setattr(forcefields, "CHARMM36_SHA256", checksum)

    installed = forcefields.install_charmm36(archive_path=archive)

    assert installed == forcefields.charmm36_directory()
    assert forcefields.discover_charmm36_directory() == installed
    assert (installed / "carb.rtp").is_file()


def test_rejects_archive_with_wrong_checksum(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    archive = tmp_path / "source.tgz"
    _write_archive(archive, tmp_path / "content")
    monkeypatch.setenv("CELLULOSE_FORGE_DATA_DIR", str(tmp_path / "data"))
    monkeypatch.setattr(forcefields, "CHARMM36_SHA256", "0" * 64)

    with pytest.raises(ValueError, match="checksum mismatch"):
        forcefields.install_charmm36(archive_path=archive)


def test_reports_missing_installation(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("CELLULOSE_FORGE_DATA_DIR", str(tmp_path / "data"))

    with pytest.raises(FileNotFoundError, match="forcefield install charmm36"):
        forcefields.discover_charmm36_directory()
