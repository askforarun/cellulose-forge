"""Experimental AMOEBA cellulose output from externally supplied parameters.

The published Welborn-group files are not redistributed.  This module accepts
only the audited cleaned beta-(1->4)-D-glucose parameter file and matching
finite oligomer topology from the cited source commit.  Crystal coordinates
replace the reference coordinates; atom types and connectivity remain exactly
those of the published Tinker topology.
"""

from __future__ import annotations

import hashlib
import math
import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cellulose_forge.crystallography import FloatArray
from cellulose_forge.geometry import BulkCrystal
from cellulose_forge.polymer_mapping import chain_template, map_crystal_to_chain_template

AMOEBA_SOURCE_REPOSITORY = "https://github.com/WelbornGroup/Newman_D-Glucose_2026"
AMOEBA_SOURCE_COMMIT = "a2adb2a3b21f4df9696fcd13cf0cc21d36ae6c71"
AMOEBA_PARAMETER_FILE = "beta_glucuse_1-4.prm"  # Published filename, including spelling.
AMOEBA_PARAMETER_SHA256 = "0f007da962a54e95ebb299e2d978e5a734e62c32ad2d3c8570dc19d4282ae169"
AMOEBA_SUPPORTED_DP = frozenset({2, 4, 6, 8, 10})
_REFERENCE_SHA256 = {
    2: "37627ac77093f1ac415ff43acf7e2f6c3a4a173c3fe8e44d0d256346000c0fc9",
    4: "856729d6677f4cc0b61b48df429d763c5840f6a5e2daff67b3c4d191f566fc5c",
    6: "7cbfc0cca05d6801a571f36916f91358e40dacbc2ef7a0eafd6806c211f1554d",
    8: "68d1cf7e0488082fd94405fbaf88746f53c9d69056f3060baad9fede9ac1ee3c",
    10: "d9cf4ef578476ba2bac1238e0630aa53eefc6196973cb1fa83f6655c0fe96efc",
}

_FIRST = (
    ("HO1", 412),
    ("O1", 400),
    ("C1", 401),
    ("H1", 413),
    ("O5", 402),
    ("C5", 403),
    ("H5", 414),
    ("C6", 404),
    ("H61", 415),
    ("H62", 415),
    ("O6", 405),
    ("HO6", 416),
    ("C4", 406),
    ("H4", 417),
    ("C3", 408),
    ("H3", 419),
    ("O3", 409),
    ("HO3", 420),
    ("C2", 410),
    ("H2", 421),
    ("O2", 411),
    ("HO2", 422),
    ("O4", 499),
)
_INTERNAL = (
    ("C1", 401),
    ("H1", 413),
    ("O5", 402),
    ("C5", 403),
    ("H5", 414),
    ("C6", 404),
    ("H61", 415),
    ("H62", 415),
    ("O6", 405),
    ("HO6", 416),
    ("C4", 406),
    ("H4", 417),
    ("C3", 408),
    ("H3", 419),
    ("O3", 409),
    ("HO3", 420),
    ("C2", 410),
    ("H2", 421),
    ("O2", 411),
    ("HO2", 422),
    ("O4", 499),
)
_LAST = (*_INTERNAL[:-1], ("O4", 407), ("HO4", 418))
_ORTHOGONAL_TOLERANCE = 1.0e-8


@dataclass(frozen=True)
class TinkerAtom:
    """One atom from an audited Tinker XYZ topology."""

    element: str
    atom_type: int
    neighbours: tuple[int, ...]


@dataclass(frozen=True)
class AmoebaSystemFiles:
    """Files written for one experimental AMOEBA system."""

    xyz: Path
    key: Path
    parameter: Path
    openmm_runner: Path | None
    source_sha256: dict[str, str]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _reference_path(root: Path, degree_polymerization: int) -> Path:
    index = degree_polymerization // 2
    return root / f"Oligomers/Dynamics_1/Vacuum/B_{index}/glc.xyz_2"


def validate_amoeba_source_directory(
    directory: str | Path, degree_polymerization: int
) -> tuple[Path, Path]:
    """Return the audited cleaned parameter and matching oligomer topology."""
    if degree_polymerization not in AMOEBA_SUPPORTED_DP:
        supported = ", ".join(str(value) for value in sorted(AMOEBA_SUPPORTED_DP))
        raise ValueError(f"experimental AMOEBA output supports finite DP {supported}")
    root = Path(directory).expanduser().resolve()
    parameter = root / AMOEBA_PARAMETER_FILE
    reference = _reference_path(root, degree_polymerization)
    missing = [str(path.relative_to(root)) for path in (parameter, reference) if not path.is_file()]
    if missing:
        raise ValueError(f"AMOEBA source directory {root} is missing {', '.join(missing)}")
    expected = {
        parameter: AMOEBA_PARAMETER_SHA256,
        reference: _REFERENCE_SHA256[degree_polymerization],
    }
    for path, wanted in expected.items():
        actual = _sha256(path)
        if actual != wanted:
            raise ValueError(
                f"{path.name} is not the audited file from {AMOEBA_SOURCE_COMMIT}: "
                f"expected {wanted}, received {actual}"
            )
    return parameter, reference


def _atom_layout(degree_polymerization: int) -> tuple[tuple[str, int, int], ...]:
    rows: list[tuple[str, int, int]] = []
    for residue in range(1, degree_polymerization + 1):
        atoms = _FIRST if residue == 1 else _LAST if residue == degree_polymerization else _INTERNAL
        rows.extend((name, atom_type, residue) for name, atom_type in atoms)
    return tuple(rows)


def _read_reference_topology(path: Path, degree_polymerization: int) -> tuple[TinkerAtom, ...]:
    lines = path.read_text(encoding="utf-8").splitlines()
    count = int(lines[0].split()[0])
    rows = lines[2:]
    if len(rows) != count or count != 21 * degree_polymerization + 3:
        raise ValueError("audited AMOEBA reference has an unexpected atom count")
    atoms: list[TinkerAtom] = []
    for expected_index, line in enumerate(rows, start=1):
        fields = line.split()
        if int(fields[0]) != expected_index:
            raise ValueError("audited AMOEBA reference atom indices are not consecutive")
        atoms.append(TinkerAtom(fields[1], int(fields[5]), tuple(int(v) - 1 for v in fields[6:])))
    layout = _atom_layout(degree_polymerization)
    if tuple(atom.atom_type for atom in atoms) != tuple(item[1] for item in layout):
        raise ValueError("audited AMOEBA reference atom types do not match the cellulose layout")
    for index, atom in enumerate(atoms):
        for neighbour in atom.neighbours:
            if neighbour < 0 or neighbour >= count or index not in atoms[neighbour].neighbours:
                raise ValueError("audited AMOEBA reference connectivity is not symmetric")
    return tuple(atoms)


def _lengths_angles(box: FloatArray) -> tuple[float, float, float, float, float, float]:
    vectors = np.asarray(box, dtype=np.float64)
    lengths = np.linalg.norm(vectors, axis=1)
    angles = []
    for first, second in ((1, 2), (0, 2), (0, 1)):
        cosine = float(np.dot(vectors[first], vectors[second]) / (lengths[first] * lengths[second]))
        angles.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
    return (
        float(lengths[0]),
        float(lengths[1]),
        float(lengths[2]),
        angles[0],
        angles[1],
        angles[2],
    )


def _require_orthorhombic(box: FloatArray) -> None:
    vectors = np.asarray(box, dtype=np.float64)
    products = [abs(float(np.dot(vectors[i], vectors[j]))) for i, j in ((0, 1), (0, 2), (1, 2))]
    if max(products) > _ORTHOGONAL_TOLERANCE:
        raise ValueError("experimental AMOEBA output currently requires an orthorhombic slab box")


def render_tinker_xyz(
    crystal: BulkCrystal, reference: tuple[TinkerAtom, ...], degree_polymerization: int
) -> str:
    """Render all crystal chains in audited Tinker atom order."""
    _require_orthorhombic(crystal.basis_vectors)
    layout = _atom_layout(degree_polymerization)
    template = chain_template((name, residue) for name, _, residue in layout)
    mapped = map_crystal_to_chain_template(crystal, template, "finite", "BGC")
    chain_numbers = tuple(dict.fromkeys(mapped.chain_numbers))
    per_chain = len(reference)
    if len(mapped.atom_names) != len(chain_numbers) * per_chain:
        raise ValueError("AMOEBA crystal mapping and reference topology have different sizes")
    coordinates = np.asarray(mapped.cartesian_coordinates, dtype=np.float64)
    lengths_angles = _lengths_angles(mapped.basis_vectors)
    lines = [f"{len(coordinates):8d}  CelluloseForge AMOEBA finite cellulose"]
    lines.append("".join(f"{value:13.6f}" for value in lengths_angles))
    for chain_index in range(len(chain_numbers)):
        offset = chain_index * per_chain
        for local, (atom, coordinate) in enumerate(
            zip(reference, coordinates[offset : offset + per_chain], strict=True)
        ):
            index = offset + local + 1
            neighbours = "".join(f"{offset + neighbour + 1:8d}" for neighbour in atom.neighbours)
            lines.append(
                f"{index:8d}  {atom.element:<2s} {coordinate[0]:13.6f} {coordinate[1]:13.6f} "
                f"{coordinate[2]:13.6f} {atom.atom_type:6d}{neighbours}"
            )
    return "\n".join(lines) + "\n"


def render_tinker_key() -> str:
    """Return a conservative AMOEBA single-point/dynamics key."""
    return """parameters amoeba_cellulose.prm
integrator stochastic
neighbor-list
polar-eps 0.00001
vdw-cutoff 12.0
ewald
ewald-cutoff 7.0
polarization mutual
"""


def render_openmm_runner() -> str:
    """Return an OpenMM runner that consumes the same native Tinker bundle."""
    return """from openmm import AmoebaVdwForce, LangevinMiddleIntegrator, unit
from openmm.app import PME, Simulation, StateDataReporter, TinkerFiles

tinker = TinkerFiles("structure.xyz", ["amoeba_cellulose.prm"])
system = tinker.createSystem(
    nonbondedMethod=PME,
    nonbondedCutoff=0.7 * unit.nanometer,
    vdwCutoff=1.2 * unit.nanometer,
    ewaldErrorTolerance=1.0e-6,
)
vdw_forces = [force for force in system.getForces() if isinstance(force, AmoebaVdwForce)]
if len(vdw_forces) != 1:
    raise RuntimeError(f"expected one AmoebaVdwForce, found {len(vdw_forces)}")
vdw_forces[0].setUseDispersionCorrection(False)
integrator = LangevinMiddleIntegrator(
    300.0 * unit.kelvin, 1.0 / unit.picosecond, 0.001 * unit.picoseconds
)
simulation = Simulation(tinker.topology, system, integrator)
simulation.context.setPositions(tinker.positions)
simulation.minimizeEnergy(maxIterations=1000)
simulation.context.setVelocitiesToTemperature(300.0 * unit.kelvin, 2026)
simulation.reporters.append(
    StateDataReporter(
        "amoeba_openmm.csv",
        100,
        step=True,
        potentialEnergy=True,
        temperature=True,
        volume=True,
    )
)
simulation.step(1000)
"""


def write_amoeba_system(
    crystal: BulkCrystal,
    source_directory: str | Path,
    output_directory: str | Path,
    degree_polymerization: int,
    *,
    write_openmm_runner: bool,
) -> AmoebaSystemFiles:
    """Write a finite cellulose Tinker bundle and optional OpenMM runner."""
    parameter, reference_path = validate_amoeba_source_directory(
        source_directory, degree_polymerization
    )
    reference = _read_reference_topology(reference_path, degree_polymerization)
    xyz_text = render_tinker_xyz(crystal, reference, degree_polymerization)
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    xyz = output / "structure.xyz"
    key = output / "amoeba.key"
    copied_parameter = output / "amoeba_cellulose.prm"
    xyz.write_text(xyz_text, encoding="utf-8")
    key.write_text(render_tinker_key(), encoding="utf-8")
    shutil.copyfile(parameter, copied_parameter)
    runner_path = output / "run_amoeba_openmm.py"
    runner = runner_path if write_openmm_runner else None
    if runner is None:
        runner_path.unlink(missing_ok=True)
    else:
        runner.write_text(render_openmm_runner(), encoding="utf-8")
    return AmoebaSystemFiles(
        xyz,
        key,
        copied_parameter,
        runner,
        {
            AMOEBA_PARAMETER_FILE: _sha256(parameter),
            f"beta_dp{degree_polymerization}_topology": _sha256(reference_path),
        },
    )
