"""Native NAMD/X-PLOR PSF output for CHARMM36 cellulose and chitin."""

from __future__ import annotations

import hashlib
import shutil
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cellulose_forge.charmm36 import (
    Charmm36TemplateAtom,
    FiniteCharmm36Topology,
    PeriodicCharmm36Topology,
    charmm36_dihedral_kinds,
    charmm36_interaction_indices,
)
from cellulose_forge.crystallography import FloatArray
from cellulose_forge.glycam_mapping import GlycamMappedSystem

CharmmTopology = PeriodicCharmm36Topology | FiniteCharmm36Topology
_ONE_FOUR_GRAPH_DISTANCE = 3
_MAXIMUM_SEGMENT_COUNT = 999
_MAXIMUM_PDB_ATOM_COUNT = 99_999
_MAXIMUM_PDB_RESIDUE_NAME_LENGTH = 4
_ORTHORHOMBIC_GRO_BOX_VALUE_COUNT = 3
_TRICLINIC_GRO_BOX_VALUE_COUNT = 9
_PDB_RESIDUE_ALIASES = {"BGLCNA": "NAG"}
CHARMM36_CARB_PARAMETER_SHA256 = "9acd5f1d887afaf5e2cdc0e420440ab912cc2b584f787120d6845dd077d58ee6"


@dataclass(frozen=True)
class NamdSystemFiles:
    """Files comprising one self-contained NAMD CHARMM36 input."""

    psf: Path
    pdb: Path
    configuration: Path
    single_point_configuration: Path
    single_point_coordinates: Path
    shared_coordinates: Path
    parameters: Path
    gromacs_single_point: Path
    single_point_runner: Path
    energy_adapter: Path


def _chain_atoms(topology: CharmmTopology) -> tuple[Charmm36TemplateAtom, ...]:
    if isinstance(topology, PeriodicCharmm36Topology):
        return tuple(
            atom
            for _residue in range(1, topology.degree_polymerization + 1)
            for atom in topology.topology_atoms
        )
    return tuple(topology.topology_atoms)


def _rows(values: tuple[tuple[int, ...], ...], *, width: int = 8, offset: int = 1) -> list[str]:
    flattened = [value + offset for item in values for value in item]
    return [
        "".join(f"{value:10d}" for value in flattened[start : start + width])
        for start in range(0, len(flattened), width)
    ]


def validate_charmm_parameter_file(path: str | Path) -> tuple[Path, str]:
    """Require the canonical February 2026 CHARMM carbohydrate parameters."""
    parameter_file = Path(path)
    if not parameter_file.is_file():
        raise ValueError("NAMD CHARMM36 output requires an existing CHARMM parameter file")
    digest = hashlib.sha256()
    with parameter_file.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    checksum = digest.hexdigest()
    if checksum != CHARMM36_CARB_PARAMETER_SHA256:
        raise ValueError(
            "NAMD CHARMM36 output requires the February 2026 par_all36_carb.prm "
            f"(expected SHA-256 {CHARMM36_CARB_PARAMETER_SHA256}, received {checksum})"
        )
    return parameter_file, checksum


def _repeat_interactions(
    interactions: tuple[tuple[int, ...], ...], atom_count: int, chain_count: int
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(index + chain * atom_count for index in interaction)
        for chain in range(chain_count)
        for interaction in interactions
    )


def _graph_one_four_pairs(
    atom_count: int, bonds: tuple[tuple[int, ...], ...]
) -> set[tuple[int, int]]:
    adjacency: list[set[int]] = [set() for _ in range(atom_count)]
    for first, second in bonds:
        adjacency[first].add(second)
        adjacency[second].add(first)
    pairs: set[tuple[int, int]] = set()
    for start in range(atom_count):
        distance = {start: 0}
        frontier = {start}
        for depth in range(1, _ONE_FOUR_GRAPH_DISTANCE + 1):
            frontier = {
                neighbour
                for atom in frontier
                for neighbour in adjacency[atom]
                if neighbour not in distance
            }
            for atom in frontier:
                distance[atom] = depth
        pairs.update(
            (min(start, atom), max(start, atom))
            for atom, depth in distance.items()
            if depth == _ONE_FOUR_GRAPH_DISTANCE
        )
    return pairs


def validate_namd_one_four_pairs(topology: CharmmTopology) -> int:
    """Require explicit CHARMM pairs to equal the bond-graph 1-4 endpoints."""
    graph = charmm36_interaction_indices(topology)
    expected = _graph_one_four_pairs(topology.atom_count, graph["bonds"])
    actual_rows = [tuple(sorted(pair)) for pair in graph["pairs"]]
    actual = set(actual_rows)
    if len(actual) != len(actual_rows):
        raise ValueError("CHARMM36 topology contains duplicate explicit 1-4 pairs")
    if actual != expected:
        missing = len(expected.difference(actual))
        extra = len(actual.difference(expected))
        raise ValueError(
            "NAMD requires CHARMM36 explicit pairs to equal graph-derived 1-4 pairs "
            f"(missing {missing}, extra {extra})"
        )
    return len(actual)


def render_psf(system: GlycamMappedSystem, topology: CharmmTopology) -> str:
    """Render an EXT X-PLOR PSF using the validated CHARMM36 interaction graph."""
    validate_namd_one_four_pairs(topology)
    chain_atoms = _chain_atoms(topology)
    atoms_per_chain = len(chain_atoms)
    if atoms_per_chain * system.chain_count != system.atom_count:
        raise ValueError("NAMD topology atom count does not match the mapped system")
    if system.chain_count > _MAXIMUM_SEGMENT_COUNT:
        raise ValueError("NAMD PSF output currently supports at most 999 chain segments")
    graph = charmm36_interaction_indices(topology)
    proper, improper = charmm36_dihedral_kinds(topology)
    bonds = _repeat_interactions(graph["bonds"], atoms_per_chain, system.chain_count)
    angles = _repeat_interactions(graph["angles"], atoms_per_chain, system.chain_count)
    dihedrals = _repeat_interactions(proper, atoms_per_chain, system.chain_count)
    impropers = _repeat_interactions(improper, atoms_per_chain, system.chain_count)
    lines = ["PSF EXT", "", "         1 !NTITLE", " REMARKS generated by CelluloseForge", ""]
    lines.append(f"{system.atom_count:10d} !NATOM")
    for index, (atom, residue_number, chain_number) in enumerate(
        zip(
            chain_atoms * system.chain_count,
            system.residue_numbers,
            system.chain_numbers,
            strict=True,
        ),
        start=1,
    ):
        local_residue = (residue_number - 1) % topology.degree_polymerization + 1
        segment = f"C{chain_number:03d}"
        lines.append(
            f"{index:10d} {segment:<8s} {local_residue:<8d} {atom.residue_name:<8s} "
            f"{atom.atom_name:<8s} {atom.atom_type:<8s} {atom.charge:14.6f} "
            f"{atom.mass:13.4f}           0"
        )
    for label, values, integers_per_line in (
        ("!NBOND: bonds", bonds, 8),
        ("!NTHETA: angles", angles, 9),
        ("!NPHI: dihedrals", dihedrals, 8),
        ("!NIMPHI: impropers", impropers, 8),
    ):
        lines.extend(("", f"{len(values):10d} {label}", *_rows(values, width=integers_per_line)))
    lines.extend(("", f"{0:10d} !NDON: donors", "", f"{0:10d} !NACC: acceptors", ""))
    lines.append(f"{0:10d} !NNB")
    lines.extend(_rows(tuple((0,) for _ in range(system.atom_count)), offset=0))
    lines.extend(
        (
            "",
            f"{1:10d}{0:10d} !NGRP",
            f"{0:10d}{0:10d}{0:10d}",
            "",
            f"{0:10d}{0:10d} !NUMLP NUMLPH",
            "",
            f"{0:10d} !NCRTERM: cross-terms",
            "",
        )
    )
    return "\n".join(lines)


def render_pdb(system: GlycamMappedSystem, topology: CharmmTopology) -> str:
    """Render mapped CHARMM atom coordinates with stable NAMD segment IDs."""
    if system.atom_count > _MAXIMUM_PDB_ATOM_COUNT:
        raise ValueError("NAMD PDB output currently supports at most 99,999 atoms")
    vectors = np.asarray(system.basis_vectors, dtype=float)
    lengths = np.linalg.norm(vectors, axis=1)
    alpha, beta, gamma = (
        np.degrees(
            np.arccos(
                np.clip(
                    np.dot(vectors[first], vectors[second]) / (lengths[first] * lengths[second]),
                    -1.0,
                    1.0,
                )
            )
        )
        for first, second in ((1, 2), (0, 2), (0, 1))
    )
    lines: list[str] = [
        f"CRYST1{lengths[0]:9.3f}{lengths[1]:9.3f}{lengths[2]:9.3f}"
        f"{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1           1"
    ]
    for index, (name, residue_name, residue_number, chain_number, coordinate) in enumerate(
        zip(
            system.atom_names,
            system.residue_names,
            system.residue_numbers,
            system.chain_numbers,
            system.cartesian_coordinates,
            strict=True,
        ),
        start=1,
    ):
        local_residue = (residue_number - 1) % topology.degree_polymerization + 1
        segment = f"C{chain_number:03d}"
        element = next((character for character in name if character.isalpha()), "X").upper()
        pdb_residue = _PDB_RESIDUE_ALIASES.get(residue_name, residue_name)
        if len(pdb_residue) > _MAXIMUM_PDB_RESIDUE_NAME_LENGTH:
            raise ValueError(
                f"NAMD PDB output has no four-character alias for residue {residue_name}"
            )
        lines.append(
            f"ATOM  {index:5d} {name:<4s} {pdb_residue:<4s}A{local_residue:4d}    "
            f"{coordinate[0]:8.3f}{coordinate[1]:8.3f}{coordinate[2]:8.3f}"
            f"  1.00  0.00      {segment:<4s}{element:>2s}"
        )
    lines.extend(("TER", "END", ""))
    return "\n".join(lines)


def _read_gro_coordinates(path: Path) -> tuple[FloatArray, FloatArray]:
    """Read the rounded coordinates and full triclinic box used by GROMACS."""
    lines = path.read_text(encoding="utf-8").splitlines()
    atom_count = int(lines[1])
    coordinates = np.asarray(
        [
            [float(line[20:28]), float(line[28:36]), float(line[36:44])]
            for line in lines[2 : 2 + atom_count]
        ],
        dtype=float,
    )
    values = [float(value) for value in lines[2 + atom_count].split()]
    if len(values) == _ORTHORHOMBIC_GRO_BOX_VALUE_COUNT:
        vectors = np.diag(values)
    elif len(values) == _TRICLINIC_GRO_BOX_VALUE_COUNT:
        vectors = np.asarray(
            [
                [values[0], values[3], values[4]],
                [values[5], values[1], values[6]],
                [values[7], values[8], values[2]],
            ],
            dtype=float,
        )
    else:
        raise ValueError("GROMACS coordinate box must contain 3 or 9 values")
    return coordinates * 10.0, vectors * 10.0


def _write_namd_binary_coordinates(path: Path, coordinates: FloatArray) -> None:
    """Write NAMD binary coordinates from the same rounded positions as the GRO file."""
    with path.open("wb") as handle:
        handle.write(struct.pack("=i", len(coordinates)))
        for coordinate in coordinates:
            handle.write(struct.pack("=ddd", *coordinate))


def _cell_lines(vectors: FloatArray) -> list[str]:
    origin = 0.5 * vectors.sum(axis=0)
    return [
        *[
            f"cellBasisVector{index} {vector[0]:.8f} {vector[1]:.8f} {vector[2]:.8f}"
            for index, vector in enumerate(vectors, start=1)
        ],
        f"cellOrigin {origin[0]:.8f} {origin[1]:.8f} {origin[2]:.8f}",
    ]


def render_namd_configuration(system: GlycamMappedSystem, *, periodic_chains: bool) -> str:
    """Render a conservative CHARMM36 PME minimization/NVT smoke-test input."""
    vectors = np.asarray(system.basis_vectors, dtype=float)
    return "\n".join(
        [
            "structure system.psf",
            "coordinates structure_namd.pdb",
            "paraTypeCharmm on",
            "parameters par_all36_carb.prm",
            "exclude scaled1-4",
            "oneFourScaling 1.0",
            "switching on",
            "switchdist 10.0",
            "cutoff 12.0",
            "pairlistdist 14.0",
            "PME yes",
            "PMEGridSpacing 1.0",
            "PMETolerance 1.0e-6",
            *_cell_lines(vectors),
            "vdwForceSwitching on",
            "LJcorrection off",
            "wrapAll off",
            f"# covalently periodic chains: {'yes' if periodic_chains else 'no'}",
            "temperature 300",
            "timestep 2.0",
            "rigidBonds all",
            "nonbondedFreq 1",
            "fullElectFrequency 2",
            "stepspercycle 10",
            "langevin on",
            "langevinTemp 300",
            "langevinDamping 1.0",
            "outputEnergies 1",
            "outputName namd",
            "binaryoutput no",
            "minimize 1000",
            "reinitvels 300",
            "run 500",
            "",
        ]
    )


def render_namd_single_point_configuration(vectors: FloatArray, *, periodic_chains: bool) -> str:
    """Render an unconstrained NAMD single-point calculation."""
    return "\n".join(
        [
            "structure system.psf",
            "coordinates structure_namd.pdb",
            "binCoordinates single_point.coor",
            "temperature 0",
            "outputName single_point",
            "paraTypeCharmm on",
            "parameters par_all36_carb.prm",
            "exclude scaled1-4",
            "oneFourScaling 1.0",
            "switching on",
            "switchdist 10.0",
            "cutoff 12.0",
            "pairlistdist 14.0",
            "PME yes",
            "PMEGridSpacing 1.0",
            "PMETolerance 1.0e-6",
            *_cell_lines(vectors),
            "vdwForceSwitching on",
            "LJcorrection off",
            "wrapAll off",
            f"# covalently periodic chains: {'yes' if periodic_chains else 'no'}",
            "rigidBonds none",
            "outputEnergies 1",
            "run 0",
            "",
        ]
    )


def render_gromacs_single_point_mdp(*, periodic_chains: bool) -> str:
    """Render a single-point GROMACS reference using the NAMD CHARMM settings."""
    return "\n".join(
        [
            "integrator = md",
            "nsteps = 0",
            "dt = 0.001",
            "continuation = yes",
            "constraints = none",
            "cutoff-scheme = Verlet",
            "nstlist = 10",
            "rlist = 1.2",
            "coulombtype = PME",
            "rcoulomb = 1.2",
            "fourierspacing = 0.10",
            "pme-order = 4",
            "ewald-rtol = 1e-6",
            "vdwtype = Cut-off",
            "vdw-modifier = Force-switch",
            "rvdw-switch = 1.0",
            "rvdw = 1.2",
            "DispCorr = no",
            "pbc = xyz",
            f"periodic-molecules = {'yes' if periodic_chains else 'no'}",
            "nstcalcenergy = 1",
            "nstenergy = 1",
            "nstxout = 0",
            "nstvout = 0",
            "nstfout = 0",
            "nstlog = 1",
            "gen-vel = no",
            "",
        ]
    )


def render_single_point_runner() -> str:
    """Render the sole portable NAMD/GROMACS single-point runner."""
    return """#!/bin/bash
set -euo pipefail

NAMD_EXE="${NAMD_EXE:-namd3}"
GMX_EXE="${GMX_EXE:-gmx}"
PYTHON_EXE="${PYTHON_EXE:-python3}"
NAMD_THREADS="${NAMD_THREADS:-8}"

if [[ ! -f system.top ]]; then
  echo "cross-engine comparison requires system.top; generate GROMACS output alongside NAMD" >&2
  exit 2
fi

"${NAMD_EXE}" +p"${NAMD_THREADS}" single_point_namd.conf > single_point_namd.log

"${GMX_EXE}" grompp -f single_point_gromacs.mdp -c single_point_coordinates.gro \
  -p system.top -o single_point_gromacs.tpr -po single_point_gromacs_processed.mdp -maxwarn 0
"${GMX_EXE}" mdrun -s single_point_gromacs.tpr -deffnm single_point_gromacs -nsteps 0
printf '%s\n' Bond U-B Proper-Dih. Improper-Dih. LJ-14 Coulomb-14 \
  'LJ-(SR)' 'Coulomb-(SR)' 'Coul.-recip.' Potential 0 | \
  "${GMX_EXE}" energy -f single_point_gromacs.edr -o single_point_gromacs.xvg
"${PYTHON_EXE}" single_point_energy.py single_point_namd.log \
  single_point_gromacs.xvg single_point_coordinates.gro
"""


def render_single_point_energy_script() -> str:
    """Render the dependency-free, failing NAMD/GROMACS energy adapter."""
    return '''"""Normalize and compare NAMD/GROMACS single-point energies."""
from __future__ import annotations
import hashlib
import json
import re
import sys
from pathlib import Path

KJ_PER_KCAL = 4.184
ABSOLUTE_TOLERANCE = 0.1
RELATIVE_TOLERANCE = 5.0e-5
def namd_energies(path: Path) -> dict[str, float]:
    fields: list[str] | None = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ETITLE:"):
            fields = line.split()[1:]
        if line.startswith("ENERGY:"):
            if fields is None:
                raise ValueError("NAMD log contains ENERGY before ETITLE")
            values = [float(value) for value in line.split()[1:]]
            if len(values) != len(fields):
                raise ValueError("NAMD ENERGY field count differs from ETITLE")
            return dict(zip(fields, values, strict=True))
    raise ValueError("NAMD log contains no ENERGY row")

def gromacs_energies(path: Path) -> dict[str, float]:
    legends: dict[int, str] = {}
    data: list[float] | None = None
    pattern = re.compile(r'@ s(\\d+) legend "([^"]+)"')
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match:
            legends[int(match.group(1))] = match.group(2)
        elif line and not line.startswith(("#", "@")):
            data = [float(value) for value in line.split()][1:]
            break
    if data is None:
        raise ValueError("GROMACS XVG contains no data row")
    return {legends[index]: value for index, value in enumerate(data)}

def require(values: dict[str, float], *names: str) -> float:
    missing = [name for name in names if name not in values]
    if missing:
        raise ValueError(f"GROMACS energy output is missing required terms: {missing}")
    return sum(values[name] for name in names)

def psf_count(marker: str) -> int:
    for line in Path("system.psf").read_text().splitlines():
        if marker in line:
            return int(line.split()[0])
    raise ValueError(f"PSF is missing {marker}")

namd = namd_energies(Path(sys.argv[1]))
gmx = gromacs_energies(Path(sys.argv[2]))
checksum = hashlib.sha256(Path(sys.argv[3]).read_bytes()).hexdigest()
gmx_groups = {
    "bond": require(gmx, "Bond"),
    "angle": require(gmx, "U-B"),
    "dihedral": require(gmx, "Proper Dih."),
    "improper": (
        require(gmx, "Improper Dih.")
        if "Improper Dih." in gmx
        else (0.0 if psf_count("!NIMPHI") == 0 else require(gmx, "Improper Dih."))
    ),
    "electrostatics": require(gmx, "Coulomb-14", "Coulomb (SR)", "Coul. recip."),
    "vdw": require(gmx, "LJ-14", "LJ (SR)"),
    "potential": require(gmx, "Potential"),
}
namd_groups = {
    "bond": namd["BOND"] * KJ_PER_KCAL,
    "angle": namd["ANGLE"] * KJ_PER_KCAL,
    "dihedral": namd["DIHED"] * KJ_PER_KCAL,
    "improper": namd["IMPRP"] * KJ_PER_KCAL,
    "electrostatics": namd["ELECT"] * KJ_PER_KCAL,
    "vdw": namd["VDW"] * KJ_PER_KCAL,
    "potential": namd["POTENTIAL"] * KJ_PER_KCAL,
}
for engine, terms in (("namd", namd_groups), ("gromacs", gmx_groups)):
    component_sum = sum(value for name, value in terms.items() if name != "potential")
    sum_tolerance = max(0.05, 1.0e-6 * abs(terms["potential"]))
    if abs(component_sum - terms["potential"]) > sum_tolerance:
        raise ValueError(f"{engine} components do not sum to Potential")
    directory = Path("single_point") / engine
    directory.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1, "engine": engine, "units": "kJ/mol",
        "coordinate_sha256": checksum, "terms": terms, "corrections": [],
    }
    (directory / "energy.json").write_text(json.dumps(payload, indent=2) + "\\n")

comparison = {}
passed = True
for name, reference in gmx_groups.items():
    observed = namd_groups[name]
    difference = observed - reference
    allowed = max(ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE * max(abs(reference), abs(observed), 1.0))
    term_passed = abs(difference) <= allowed
    passed &= term_passed
    comparison[name] = {"gromacs_kj_mol": reference, "namd_kj_mol": observed,
                        "difference_kj_mol": difference, "allowed_kj_mol": allowed,
                        "passed": term_passed}
report = {"schema_version": 1, "coordinate_sha256": checksum,
          "reference_engine": "gromacs", "comparison": comparison, "passed": passed}
Path("single_point_energy_comparison.json").write_text(json.dumps(report, indent=2) + "\\n")
if not passed:
    raise SystemExit("single-point energy comparison failed")
'''


def write_namd_system(  # noqa: PLR0913
    system: GlycamMappedSystem,
    topology: CharmmTopology,
    parameter_file: str | Path,
    gromacs_coordinates: str | Path,
    output_directory: str | Path,
    *,
    periodic_chains: bool,
) -> NamdSystemFiles:
    """Write PSF, PDB, parameter file, and a runnable NAMD smoke-test input."""
    source_parameters, _checksum = validate_charmm_parameter_file(parameter_file)
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=True)
    psf = directory / "system.psf"
    pdb = directory / "structure_namd.pdb"
    configuration = directory / "run_namd.conf"
    single_point_configuration = directory / "single_point_namd.conf"
    single_point_coordinates = directory / "single_point.coor"
    shared_coordinates = directory / "single_point_coordinates.gro"
    parameters = directory / "par_all36_carb.prm"
    gromacs_single_point = directory / "single_point_gromacs.mdp"
    single_point_runner = directory / "run_single_point_comparison.sh"
    energy_adapter = directory / "single_point_energy.py"
    shutil.copy2(gromacs_coordinates, shared_coordinates)
    rounded_coordinates, rounded_vectors = _read_gro_coordinates(shared_coordinates)
    _write_namd_binary_coordinates(single_point_coordinates, rounded_coordinates)
    psf.write_text(render_psf(system, topology), encoding="utf-8")
    pdb.write_text(render_pdb(system, topology), encoding="utf-8")
    configuration.write_text(
        render_namd_configuration(system, periodic_chains=periodic_chains), encoding="utf-8"
    )
    single_point_configuration.write_text(
        render_namd_single_point_configuration(rounded_vectors, periodic_chains=periodic_chains),
        encoding="utf-8",
    )
    gromacs_single_point.write_text(
        render_gromacs_single_point_mdp(periodic_chains=periodic_chains), encoding="utf-8"
    )
    single_point_runner.write_text(render_single_point_runner(), encoding="utf-8")
    single_point_runner.chmod(single_point_runner.stat().st_mode | 0o111)
    energy_adapter.write_text(render_single_point_energy_script(), encoding="utf-8")
    shutil.copy2(source_parameters, parameters)
    return NamdSystemFiles(
        psf=psf,
        pdb=pdb,
        configuration=configuration,
        single_point_configuration=single_point_configuration,
        single_point_coordinates=single_point_coordinates,
        shared_coordinates=shared_coordinates,
        parameters=parameters,
        gromacs_single_point=gromacs_single_point,
        single_point_runner=single_point_runner,
        energy_adapter=energy_adapter,
    )
