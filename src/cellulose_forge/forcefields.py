"""Installation and discovery of external force-field distributions."""

from __future__ import annotations

import hashlib
import os
import shutil
import tarfile
import tempfile
import urllib.request
from pathlib import Path, PurePosixPath

CHARMM36_RELEASE = "charmm36-feb2026_cgenff-5.0.ff"
CHARMM36_ARCHIVE = f"{CHARMM36_RELEASE}.tgz"
CHARMM36_URL = (
    "https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/"
    f"{CHARMM36_ARCHIVE}"
)
CHARMM36_SHA256 = "c2b0c70e21cced150528ee25d205e5bd3367a5e55e85c2ff0403911980528151"


def forcefield_data_directory() -> Path:
    """Return the user-local CelluloseForge force-field directory."""
    override = os.environ.get("CELLULOSE_FORGE_DATA_DIR")
    if override:
        return Path(override).expanduser().resolve() / "forcefields"
    xdg_data_home = os.environ.get("XDG_DATA_HOME")
    base = Path(xdg_data_home).expanduser() if xdg_data_home else Path.home() / ".local/share"
    return base / "cellulose-forge/forcefields"


def charmm36_directory() -> Path:
    """Return the expected installation directory for the supported release."""
    return forcefield_data_directory() / CHARMM36_RELEASE


def validate_charmm36_directory(directory: Path) -> Path:
    """Validate the files required by the native CHARMM36 adapter."""
    required = (
        "forcefield.itp",
        "forcefield.doc",
        "ffbonded.itp",
        "ffnonbonded.itp",
        "nbfix.itp",
        "carb.rtp",
        "carb.hdb",
    )
    missing = [name for name in required if not (directory / name).is_file()]
    if missing:
        msg = f"incomplete CHARMM36 installation at {directory}: missing {', '.join(missing)}"
        raise ValueError(msg)
    return directory


def discover_charmm36_directory() -> Path:
    """Locate and validate the supported user-local CHARMM36 release."""
    directory = charmm36_directory()
    if not directory.is_dir():
        msg = (
            f"CHARMM36 {CHARMM36_RELEASE} is not installed; run "
            "'cellulose-forge forcefield install charmm36'"
        )
        raise FileNotFoundError(msg)
    return validate_charmm36_directory(directory)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_archive_members(archive: tarfile.TarFile) -> None:
    expected_root = PurePosixPath(CHARMM36_RELEASE)
    for member in archive.getmembers():
        member_path = PurePosixPath(member.name)
        if (
            member_path.is_absolute()
            or ".." in member_path.parts
            or not member_path.parts
            or member_path.parts[0] != expected_root.name
        ):
            msg = f"unsafe or unexpected path in CHARMM36 archive: {member.name}"
            raise ValueError(msg)
        if member.issym() or member.islnk() or member.isdev():
            msg = f"unsupported special entry in CHARMM36 archive: {member.name}"
            raise ValueError(msg)


def install_charmm36(*, archive_path: Path | None = None) -> Path:
    """Install the pinned CHARMM36 release from its official archive or a local copy."""
    destination = charmm36_directory()
    if destination.exists():
        return validate_charmm36_directory(destination)

    parent = destination.parent
    parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".charmm36-install-", dir=parent) as temporary:
        temporary_path = Path(temporary)
        downloaded = temporary_path / CHARMM36_ARCHIVE
        if archive_path is None:
            with urllib.request.urlopen(CHARMM36_URL) as response, downloaded.open("wb") as out:
                shutil.copyfileobj(response, out)
        else:
            shutil.copyfile(archive_path, downloaded)

        actual_sha256 = _sha256(downloaded)
        if actual_sha256 != CHARMM36_SHA256:
            msg = (
                "CHARMM36 archive checksum mismatch: "
                f"expected {CHARMM36_SHA256}, received {actual_sha256}"
            )
            raise ValueError(msg)

        extraction = temporary_path / "extracted"
        extraction.mkdir()
        with tarfile.open(downloaded, "r:gz") as archive:
            _validate_archive_members(archive)
            archive.extractall(extraction, filter="data")
        installed = validate_charmm36_directory(extraction / CHARMM36_RELEASE)
        installed.rename(destination)

    return destination
