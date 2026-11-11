"""Environment diagnostics for supported coordinate and simulation outputs."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version

import numpy


@dataclass(frozen=True)
class DiagnosticCheck:
    """One environment check: what was checked, whether it passed, and why."""

    name: str
    required_for: str
    ok: bool
    detail: str
    fix: str | None = None


def _check_python() -> DiagnosticCheck:
    version_info = sys.version_info
    ok = version_info >= (3, 11)
    detail = f"Python {version_info.major}.{version_info.minor}.{version_info.micro}"
    return DiagnosticCheck(
        name="python",
        required_for="coordinate-only and GLYCAM builds",
        ok=ok,
        detail=detail,
        fix=None if ok else "install Python 3.11 or newer",
    )


def _check_numpy() -> DiagnosticCheck:
    return DiagnosticCheck(
        name="numpy",
        required_for="coordinate-only and GLYCAM builds",
        ok=True,
        detail=f"numpy {numpy.__version__}",
    )


def _check_cellulose_forge_version() -> DiagnosticCheck:
    try:
        installed = version("cellulose-forge")
    except PackageNotFoundError:
        installed = "unknown (running from a source checkout without an installed package)"
    return DiagnosticCheck(
        name="cellulose-forge",
        required_for="all builds",
        ok=True,
        detail=f"version {installed}",
    )


def _check_tleap(executable: str) -> DiagnosticCheck:
    resolved = shutil.which(executable)
    if resolved is None:
        return DiagnosticCheck(
            name="tleap",
            required_for="--force-field glycam builds",
            ok=False,
            detail=f"'{executable}' was not found on PATH",
            fix=(
                "conda create -n cellulose-forge -c conda-forge python=3.12 ambertools "
                "parmed pip, then conda activate cellulose-forge"
            ),
        )
    return DiagnosticCheck(
        name="tleap",
        required_for="--force-field glycam builds",
        ok=True,
        detail=f"found at {resolved}",
    )


def _check_parmed(python_executable: str | None) -> DiagnosticCheck:
    interpreter = python_executable or sys.executable
    try:
        probe = subprocess.run(
            [interpreter, "-c", "import parmed; print(parmed.__version__)"],
            capture_output=True,
            text=True,
            check=False,
            timeout=15,
        )
    except OSError as error:
        return DiagnosticCheck(
            name="parmed",
            required_for="--force-field glycam builds",
            ok=False,
            detail=f"could not run interpreter {interpreter}: {error}",
            fix=(
                "pass an existing Python interpreter containing ParmEd with "
                "--parmed-python"
            ),
        )
    except subprocess.TimeoutExpired:
        return DiagnosticCheck(
            name="parmed",
            required_for="--force-field glycam builds",
            ok=False,
            detail=f"ParmEd import probe timed out with interpreter {interpreter}",
            fix="check that the selected Python interpreter starts normally",
        )
    if probe.returncode != 0:
        return DiagnosticCheck(
            name="parmed",
            required_for="--force-field glycam builds",
            ok=False,
            detail=f"not importable with interpreter {interpreter}",
            fix=(
                "install ParmEd in the interpreter passed to --parmed-python "
                "(conda install -c conda-forge parmed, or: python -m pip install parmed)"
            ),
        )
    return DiagnosticCheck(
        name="parmed",
        required_for="--force-field glycam builds",
        ok=True,
        detail=f"parmed {probe.stdout.strip()} via {interpreter}",
    )


def _check_gromacs() -> DiagnosticCheck:
    resolved = shutil.which("gmx")
    if resolved is None:
        return DiagnosticCheck(
            name="gmx",
            required_for="preprocessing and simulating generated GLYCAM systems",
            ok=False,
            detail="'gmx' was not found on PATH",
            fix="install GROMACS, or transfer generated files to a machine that has it",
        )
    return DiagnosticCheck(
        name="gmx",
        required_for="preprocessing and simulating generated GLYCAM systems",
        ok=True,
        detail=f"found at {resolved}",
    )


def _check_openmm(python_executable: str | None) -> DiagnosticCheck:
    interpreter = python_executable or sys.executable
    try:
        probe = subprocess.run(
            [interpreter, "-c", "import openmm; print(openmm.__version__)"],
            capture_output=True,
            text=True,
            check=False,
            timeout=15,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        return DiagnosticCheck(
            name="openmm",
            required_for="--md-engine openmm builds",
            ok=False,
            detail=f"could not probe OpenMM with {interpreter}: {error}",
            fix="install cellulose-forge[openmm] or pass --openmm-python",
        )
    if probe.returncode != 0:
        return DiagnosticCheck(
            name="openmm",
            required_for="--md-engine openmm builds",
            ok=False,
            detail=f"not importable with interpreter {interpreter}",
            fix="install cellulose-forge[openmm] or pass --openmm-python",
        )
    return DiagnosticCheck(
        name="openmm",
        required_for="--md-engine openmm builds",
        ok=True,
        detail=f"OpenMM {probe.stdout.strip()} via {interpreter}",
    )


def run_diagnostics(
    *,
    tleap_executable: str = "tleap",
    parmed_python: str | None = None,
    openmm_python: str | None = None,
) -> tuple[DiagnosticCheck, ...]:
    """Run every environment check and return the results in report order."""
    return (
        _check_python(),
        _check_numpy(),
        _check_cellulose_forge_version(),
        _check_tleap(tleap_executable),
        _check_parmed(parmed_python),
        _check_gromacs(),
        _check_openmm(openmm_python),
    )


def render_text(checks: tuple[DiagnosticCheck, ...]) -> str:
    """Render diagnostics as human- and agent-readable plain text."""
    lines = ["CelluloseForge environment diagnostics", ""]
    for check in checks:
        status = "OK" if check.ok else "MISSING"
        lines.append(f"[{status}] {check.name} ({check.required_for}): {check.detail}")
        if not check.ok and check.fix:
            lines.append(f"        fix: {check.fix}")
    coordinate_ready = all(
        check.ok for check in checks if check.name in {"python", "numpy", "cellulose-forge"}
    )
    glycam_ready = coordinate_ready and all(
        check.ok for check in checks if check.name in {"tleap", "parmed"}
    )
    simulation_ready = glycam_ready and all(check.ok for check in checks if check.name == "gmx")
    openmm_ready = glycam_ready and all(
        check.ok for check in checks if check.name == "openmm"
    )
    lines.extend(
        (
            "",
            f"coordinate-only builds ready: {coordinate_ready}",
            f"GLYCAM/GROMACS-file builds ready: {glycam_ready}",
            f"GROMACS preprocessing/simulation ready: {simulation_ready}",
            f"OpenMM builds ready: {openmm_ready}",
        )
    )
    return "\n".join(lines)


def render_json(checks: tuple[DiagnosticCheck, ...]) -> str:
    """Render diagnostics as a single JSON object for programmatic consumption."""
    coordinate_ready = all(
        check.ok for check in checks if check.name in {"python", "numpy", "cellulose-forge"}
    )
    glycam_ready = coordinate_ready and all(
        check.ok for check in checks if check.name in {"tleap", "parmed"}
    )
    simulation_ready = glycam_ready and all(check.ok for check in checks if check.name == "gmx")
    openmm_ready = glycam_ready and all(
        check.ok for check in checks if check.name == "openmm"
    )
    payload = {
        "checks": [
            {
                "name": check.name,
                "required_for": check.required_for,
                "ok": check.ok,
                "detail": check.detail,
                "fix": check.fix,
            }
            for check in checks
        ],
        "coordinate_only_ready": coordinate_ready,
        "glycam_gromacs_files_ready": glycam_ready,
        "gromacs_simulation_ready": simulation_ready,
        "openmm_ready": openmm_ready,
    }
    return json.dumps(payload, indent=2, sort_keys=True)
