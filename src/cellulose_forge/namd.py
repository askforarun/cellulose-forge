"""Native NAMD/X-PLOR PSF output for CHARMM36 cellulose and chitin."""

from __future__ import annotations

import hashlib
import shutil
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
from cellulose_forge.glycam_mapping import GlycamMappedSystem

CharmmTopology = PeriodicCharmm36Topology | FiniteCharmm36Topology
_ONE_FOUR_GRAPH_DISTANCE = 3
_MAXIMUM_SEGMENT_COUNT = 999
_MAXIMUM_PDB_ATOM_COUNT = 99_999
_MAXIMUM_PDB_RESIDUE_NAME_LENGTH = 4
_PDB_RESIDUE_ALIASES = {"BGLCNA": "NAG"}
CHARMM36_CARB_PARAMETER_SHA256 = "9acd5f1d887afaf5e2cdc0e420440ab912cc2b584f787120d6845dd077d58ee6"


@dataclass(frozen=True)
class NamdSystemFiles:
    """Files comprising one self-contained NAMD CHARMM36 input."""

    psf: Path
    pdb: Path
    configuration: Path
    parameters: Path
    gromacs_reference: Path
    slurm: Path
    energy_comparison: Path


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


def render_namd_configuration(system: GlycamMappedSystem, *, periodic_chains: bool) -> str:
    """Render a conservative CHARMM36 PME minimization/NVT smoke-test input."""
    vectors = np.asarray(system.basis_vectors, dtype=float)
    origin = 0.5 * vectors.sum(axis=0)
    vector_lines = [
        f"cellBasisVector{index} {vector[0]:.8f} {vector[1]:.8f} {vector[2]:.8f}"
        for index, vector in enumerate(vectors, start=1)
    ]
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
            *vector_lines,
            f"cellOrigin {origin[0]:.8f} {origin[1]:.8f} {origin[2]:.8f}",
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
            "run 0",
            "minimize 1000",
            "reinitvels 300",
            "run 500",
            "",
        ]
    )


def render_gromacs_t0_mdp(*, periodic_chains: bool) -> str:
    """Render a zero-step GROMACS reference using the NAMD CHARMM settings."""
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


def render_slurm_script() -> str:
    """Render a portable Slurm job for NAMD and optional GROMACS t=0 energies."""
    return """#!/bin/bash
#SBATCH --job-name=cf_namd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

# Set module names only when the target cluster needs them, for example:
# export NAMD_MODULE=apps/namd/3
# export GROMACS_MODULE=apps/gromacs_cuda/2026.1
[[ -z "${NAMD_MODULE:-}" ]] || module load "${NAMD_MODULE}"
[[ -z "${GROMACS_MODULE:-}" ]] || module load "${GROMACS_MODULE}"

NAMD_EXE="${NAMD_EXE:-namd3}"
GMX_EXE="${GMX_EXE:-gmx}"

"${NAMD_EXE}" +p"${SLURM_CPUS_PER_TASK}" run_namd.conf > namd.log

if [[ -f system.top && -f structure.gro ]]; then
  "${GMX_EXE}" grompp -f gromacs_namd_t0.mdp -c structure_namd.pdb \
    -p system.top -o gromacs_t0.tpr -po gromacs_t0_processed.mdp -maxwarn 0
  "${GMX_EXE}" mdrun -s gromacs_t0.tpr -deffnm gromacs_t0 -nsteps 0
  printf '%s\n' Bond U-B Angle Proper-Dih. Improper-Dih. LJ-14 Coulomb-14 \
    LJ-SR Coulomb-SR Coul.-recip. Potential 0 | \
    "${GMX_EXE}" energy -f gromacs_t0.edr -o gromacs_t0_energy.xvg
  python compare_t0.py namd.log gromacs_t0_energy.xvg > t0_energy_comparison.json
else
  echo 'GROMACS files absent: NAMD completed; cross-engine t=0 comparison skipped.'
fi
"""


def render_t0_comparison_script() -> str:
    """Render the dependency-free NAMD/GROMACS t=0 energy comparison script."""
    return '''"""Compare first NAMD ENERGY row with a GROMACS energy XVG."""
from __future__ import annotations
import json
import re
import sys
from pathlib import Path

KJ_PER_KCAL = 4.184
NAMD_FIELDS = (
    "TS", "BOND", "ANGLE", "DIHED", "IMPRP", "ELECT", "VDW", "BOUNDARY",
    "MISC", "KINETIC", "TOTAL", "TEMP", "POTENTIAL", "TOTAL3", "TEMPAVG",
    "PRESSURE", "GPRESSURE", "VOLUME", "PRESSAVG", "GPRESSAVG",
)

def namd_energies(path: Path) -> dict[str, float]:
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ENERGY:"):
            values = [float(value) for value in line.split()[1:]]
            return dict(zip(NAMD_FIELDS, values, strict=False))
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
    return {legends[index]: value / KJ_PER_KCAL for index, value in enumerate(data)}

def pick(values: dict[str, float], *names: str) -> float:
    normalized = {key.replace(" ", "").lower(): value for key, value in values.items()}
    return sum(normalized.get(name.replace(" ", "").lower(), 0.0) for name in names)

namd = namd_energies(Path(sys.argv[1]))
gmx = gromacs_energies(Path(sys.argv[2]))
gmx_groups = {
    "BOND": pick(gmx, "Bond"),
    "ANGLE": pick(gmx, "Angle", "U-B"),
    "DIHED": pick(gmx, "Proper Dih.", "Proper-Dih."),
    "IMPRP": pick(gmx, "Improper Dih.", "Improper-Dih."),
    "ELECT": pick(gmx, "Coulomb-14", "Coulomb (SR)", "Coulomb-SR", "Coul. recip."),
    "VDW": pick(gmx, "LJ-14", "LJ (SR)", "LJ-SR"),
    "POTENTIAL": pick(gmx, "Potential"),
}
comparison = {
    key: {
        "namd_kcal_mol": namd[key],
        "gromacs_kcal_mol": value,
        "difference_kcal_mol": namd[key] - value,
    }
    for key, value in gmx_groups.items()
}
print(json.dumps({"namd": namd, "gromacs_kcal_mol": gmx, "comparison": comparison}, indent=2))
'''


def write_namd_system(
    system: GlycamMappedSystem,
    topology: CharmmTopology,
    parameter_file: str | Path,
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
    parameters = directory / "par_all36_carb.prm"
    gromacs_reference = directory / "gromacs_namd_t0.mdp"
    slurm = directory / "run_cross_engine.sbatch"
    energy_comparison = directory / "compare_t0.py"
    psf.write_text(render_psf(system, topology), encoding="utf-8")
    pdb.write_text(render_pdb(system, topology), encoding="utf-8")
    configuration.write_text(
        render_namd_configuration(system, periodic_chains=periodic_chains), encoding="utf-8"
    )
    gromacs_reference.write_text(
        render_gromacs_t0_mdp(periodic_chains=periodic_chains), encoding="utf-8"
    )
    slurm.write_text(render_slurm_script(), encoding="utf-8")
    energy_comparison.write_text(render_t0_comparison_script(), encoding="utf-8")
    shutil.copy2(source_parameters, parameters)
    return NamdSystemFiles(
        psf=psf,
        pdb=pdb,
        configuration=configuration,
        parameters=parameters,
        gromacs_reference=gromacs_reference,
        slurm=slurm,
        energy_comparison=energy_comparison,
    )
