"""Experimental CHARMM Drude polarizable cellulose output for NAMD and OpenMM.

The topology is generated natively by :mod:`cellulose_forge.drude_topology`;
no CHARMM, VMD or psfgen executable is used.
"""

from __future__ import annotations

import shutil
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cellulose_forge.crystallography import FloatArray
from cellulose_forge.drude_topology import (
    DRUDE_CARBOHYDRATE_FILE,
    DRUDE_MAIN_FILE,
    DrudeChainTopology,
    build_drude_cellulose_chain,
    read_drude_cellulose_definitions,
    verify_audited_streams,
)
from cellulose_forge.geometry import BulkCrystal
from cellulose_forge.polymer_mapping import (
    chain_template,
    map_crystal_to_chain_template,
    minimum_image,
)

DRUDE_RELEASE = "2023a"
DRUDE_RESIDUE_NAME = "BGLC"
_MAXIMUM_SEGMENT_COUNT = 999
_PDB_SERIAL_MODULUS = 100_000
_PDB_RESIDUE_MODULUS = 10_000
_PSF_KIND_FLAG = {"atom": 0, "lonepair": -1, "drude": -2}
_DEGENERATE_GEOMETRY = 1.0e-8
_KJ_NM2_TO_NAMD_KCAL_A2 = 1.0 / 836.8


@dataclass(frozen=True)
class DrudeSystem:
    """Every particle of a Drude system: repeated chains in PSF order."""

    topology: DrudeChainTopology
    chain_numbers: tuple[int, ...]
    cartesian_coordinates: FloatArray
    basis_vectors: FloatArray

    @property
    def chain_count(self) -> int:
        """Return the number of chains."""
        return len(self.chain_numbers)

    @property
    def atom_count(self) -> int:
        """Return the particle count including Drude particles and lone pairs."""
        return self.chain_count * len(self.topology.sites)


@dataclass(frozen=True)
class DrudeSystemFiles:
    """Files written for one experimental CHARMM Drude system."""

    psf: Path
    pdb: Path
    binary_coordinates: Path
    configuration: Path | None
    openmm_runner: Path | None
    restraint_reference: Path | None
    main_stream: Path
    carbohydrate_stream: Path
    source_sha256: dict[str, str]
    counts: dict[str, int]


def drude_heavy_atom_restraint_indices(
    system: DrudeSystem, chain_ids: tuple[int, ...]
) -> tuple[int, ...]:
    """Return zero-based real-heavy-atom indices for selected crystal chains."""
    selected = set(chain_ids)
    per_chain = len(system.topology.sites)
    return tuple(
        chain_index * per_chain + site_index
        for chain_index, chain_id in enumerate(system.chain_numbers)
        if chain_id in selected
        for site_index, site in enumerate(system.topology.sites)
        if site.kind == "atom" and not site.name.upper().startswith("H")
    )


def validate_drude_toppar_directory(directory: str | Path) -> tuple[Path, Path]:
    """Return the two audited official 2023a stream files required by the backend."""
    root = Path(directory).expanduser().resolve()
    main = root / DRUDE_MAIN_FILE
    carbohydrate = root / DRUDE_CARBOHYDRATE_FILE
    missing = [path.name for path in (main, carbohydrate) if not path.is_file()]
    if missing:
        raise ValueError(
            f"CHARMM Drude directory {root} is missing {', '.join(missing)}; "
            "use the official MacKerell drude_toppar 2023 release"
        )
    verify_audited_streams(main, carbohydrate)
    return main, carbohydrate


def lone_pair_position(  # noqa: PLR0913, PLR0917 - mirrors the CHARMM LONEPAIR record
    mode: str,
    host: FloatArray,
    second: FloatArray,
    third: FloatArray,
    distance: float,
    angle_degrees: float,
    dihedral_degrees: float,
) -> FloatArray:
    """Place one CHARMM ``relative`` or ``bisector`` lone pair.

    The site lies ``distance`` from ``host``; ``angle`` is site--host--second
    and ``dihedral`` is site--host--second--third.  For ``bisector`` the second
    point is the midpoint of ``second`` and ``third``.  NAMD and OpenMM
    recompute lone pairs from the PSF record, so this fixes only the written
    starting coordinates.
    """
    if mode == "bisector":
        second = 0.5 * (second + third)
    elif mode != "relative":
        raise ValueError(f"unsupported lone-pair mode {mode!r}")
    axis = host - second
    length = float(np.linalg.norm(axis))
    if length < _DEGENERATE_GEOMETRY:
        raise ValueError("lone-pair hosts coincide")
    axis = axis / length
    normal = np.cross(second - third, axis)
    normal_length = float(np.linalg.norm(normal))
    if normal_length < _DEGENERATE_GEOMETRY:
        raise ValueError("lone-pair hosts are collinear")
    normal = normal / normal_length
    in_plane = np.cross(normal, axis)
    theta = np.radians(angle_degrees)
    phi = np.radians(dihedral_degrees)
    offset = (
        -np.cos(theta) * axis
        + np.sin(theta) * np.cos(phi) * in_plane
        + np.sin(theta) * np.sin(phi) * normal
    )
    return np.asarray(host + distance * offset, dtype=np.float64)


def build_drude_system(
    crystal: BulkCrystal, topology: DrudeChainTopology, chain_boundary: str
) -> DrudeSystem:
    """Map crystal atoms onto the Drude chain and place Drudes and lone pairs."""
    if (chain_boundary == "periodic") != topology.periodic:
        raise ValueError("Drude topology and chain boundary disagree")
    real = [index for index, site in enumerate(topology.sites) if site.kind == "atom"]
    template = chain_template(
        (topology.sites[index].name, topology.sites[index].residue) for index in real
    )
    mapped = map_crystal_to_chain_template(crystal, template, chain_boundary, DRUDE_RESIDUE_NAME)
    chain_numbers = tuple(dict.fromkeys(mapped.chain_numbers))
    box = np.asarray(mapped.basis_vectors, dtype=np.float64)
    per_chain = len(topology.sites)
    coordinates = np.zeros((len(chain_numbers) * per_chain, 3), dtype=np.float64)
    physical = np.asarray(mapped.cartesian_coordinates, dtype=np.float64).reshape(
        len(chain_numbers), len(real), 3
    )
    for chain in range(len(chain_numbers)):
        local = coordinates[chain * per_chain : (chain + 1) * per_chain]
        local[real] = physical[chain]
        for drude, parent in topology.drude_parents:
            local[drude] = local[parent]
        for lone_pair in topology.lone_pairs:
            host = local[lone_pair.hosts[0]]
            second = host + minimum_image(local[lone_pair.hosts[1]] - host, box)
            third = host + minimum_image(local[lone_pair.hosts[2]] - host, box)
            local[lone_pair.site] = lone_pair_position(
                lone_pair.mode,
                host,
                second,
                third,
                lone_pair.distance,
                lone_pair.angle,
                lone_pair.dihedral,
            )
    coordinates.setflags(write=False)
    return DrudeSystem(topology, chain_numbers, coordinates, mapped.basis_vectors)


def _integer_rows(values: list[int], width: int) -> list[str]:
    return [
        "".join(f"{value:10d}" for value in values[start : start + width])
        for start in range(0, len(values), width)
    ]


def render_drude_psf(system: DrudeSystem) -> str:
    """Render a ``PSF EXT DRUDE`` file with lone-pair and anisotropy sections."""
    topology = system.topology
    per_chain = len(topology.sites)
    if system.chain_count > _MAXIMUM_SEGMENT_COUNT:
        raise ValueError("NAMD PSF output currently supports at most 999 chain segments")
    offsets = [chain * per_chain for chain in range(system.chain_count)]

    def repeated(items: tuple[tuple[int, ...], ...]) -> list[int]:
        return [index + offset + 1 for offset in offsets for item in items for index in item]

    title_lines = [
        " REMARKS generated natively by CelluloseForge from CHARMM Drude BGLC and 14bb",
        f" REMARKS Drude charges use a Coulomb constant of {topology.coulomb_constant:.4f}",
    ]
    if topology.periodic:
        title_lines.append(" REMARKS PERIODIC CHAINS ARE CONSTRUCTION ONLY; DO NOT RUN")
    lines = [
        "PSF EXT DRUDE",
        "",
        f"{len(title_lines):10d} !NTITLE",
        *title_lines,
        "",
        f"{system.atom_count:10d} !NATOM",
    ]
    index = 0
    for chain_number in system.chain_numbers:
        segment = f"C{chain_number:03d}"
        for site in topology.sites:
            index += 1
            lines.append(
                f"{index:10d} {segment:<8s} {site.residue:<8d} {DRUDE_RESIDUE_NAME:<8s} "
                f"{site.name:<8s} {site.atom_type:<8s} {site.charge:14.6f} {site.mass:13.4f} "
                f"{_PSF_KIND_FLAG[site.kind]:11d} {site.alpha:14.6f} {site.thole:14.6f}"
            )
    drude_bond_count = len(topology.drude_parents)
    ordinary_bonds = topology.bonds[: len(topology.bonds) - drude_bond_count]
    drude_bonds = topology.bonds[len(topology.bonds) - drude_bond_count :]
    if any(
        "drude" in {topology.sites[first].kind, topology.sites[second].kind}
        for first, second in ordinary_bonds
    ) or any(topology.sites[second].kind != "drude" for _first, second in drude_bonds):
        raise ValueError("Drude topology must list parent--Drude bonds after all other bonds")
    # Across the whole system too: NAMD needs every parent--Drude bond last.
    bond_indices = [*repeated(ordinary_bonds), *repeated(drude_bonds)]
    lines.extend(
        (
            "",
            f"{len(bond_indices) // 2:10d} !NBOND: bonds",
            *_integer_rows(bond_indices, 8),
        )
    )
    for label, items, width in (
        ("!NTHETA: angles", topology.angles, 9),
        ("!NPHI: dihedrals", topology.dihedrals, 8),
    ):
        lines.extend(
            (
                "",
                f"{len(items) * system.chain_count:10d} {label}",
                *_integer_rows(repeated(items), width),
            )
        )
    lines.extend(
        (
            "",
            f"{0:10d} !NIMPHI: impropers",
            "",
            f"{0:10d} !NDON: donors",
            "",
            f"{0:10d} !NACC: acceptors",
            "",
            f"{0:10d} !NNB",
            "",
            *_integer_rows([0] * system.atom_count, 8),
            "",
            f"{1:10d}{0:10d} !NGRP",
            f"{0:10d}{0:10d}{0:10d}",
            "",
        )
    )
    lone_pair_count = len(topology.lone_pairs) * system.chain_count
    lines.append(f"{lone_pair_count:10d}{4 * lone_pair_count:10d} !NUMLP NUMLPH")
    pointer = 1
    for _offset in offsets:
        for lone_pair in topology.lone_pairs:
            # CHARMM and NAMD mark a bisector lone pair by a negative distance.
            distance = -lone_pair.distance if lone_pair.mode == "bisector" else lone_pair.distance
            lines.append(
                f"{3:10d}{pointer:10d}   F{distance:14.4f}"
                f"{lone_pair.angle:14.4f}{lone_pair.dihedral:14.4f}"
            )
            pointer += 4
    lines.extend(
        _integer_rows(repeated(tuple((item.site, *item.hosts) for item in topology.lone_pairs)), 8)
    )
    lines.extend(("", f"{len(topology.anisotropies) * system.chain_count:10d} !NUMANISO"))
    for _offset in offsets:
        lines.extend(
            f"{item.k11:20.4f}{item.k22:14.4f}{item.k33:14.4f}" for item in topology.anisotropies
        )
    lines.extend(_integer_rows(repeated(tuple(item.atoms for item in topology.anisotropies)), 8))
    lines.append("")
    return "\n".join(lines)


def render_drude_pdb(system: DrudeSystem) -> str:
    """Render a NAMD coordinate PDB; full precision is in the binary coordinates."""
    topology = system.topology
    pdb_names: list[str] = []
    names_by_residue: dict[int, set[str]] = {}
    for site in topology.sites:
        # The PDB atom-name field is exactly four columns. CHARMM lone-pair
        # names such as LPDA4 are retained in the PSF but shortened here.
        pdb_name = site.name[:4]
        residue_names = names_by_residue.setdefault(site.residue, set())
        if pdb_name in residue_names:
            raise ValueError(
                f"Drude PDB atom-name shortening is ambiguous in residue {site.residue}: "
                f"{site.name!r} becomes {pdb_name!r}"
            )
        residue_names.add(pdb_name)
        pdb_names.append(pdb_name)
    lines = ["REMARK CHARMM Drude system; structure_namd.coor holds full-precision coordinates"]
    index = 0
    for chain_number in system.chain_numbers:
        segment = f"C{chain_number:03d}"
        for site, pdb_name in zip(topology.sites, pdb_names, strict=True):
            x, y, z = system.cartesian_coordinates[index]
            index += 1
            element = "" if site.kind != "atom" else site.name[0]
            lines.append(
                f"ATOM  {index % _PDB_SERIAL_MODULUS:5d} {pdb_name:<4s} {DRUDE_RESIDUE_NAME:<4s}"
                f"A{site.residue % _PDB_RESIDUE_MODULUS:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
                f"  1.00  0.00      {segment:<4s}{element:>2s}"
            )
    lines.extend(("END", ""))
    return "\n".join(lines)


def write_namd_binary_coordinates(path: Path, coordinates: FloatArray) -> None:
    """Write NAMD binary coordinates in double precision."""
    with path.open("wb") as handle:
        handle.write(struct.pack("=i", len(coordinates)))
        for coordinate in coordinates:
            handle.write(struct.pack("=ddd", *coordinate))


def _write_run_configuration(
    directory: Path,
    system: DrudeSystem,
    *,
    lj_correction: bool,
    write_run_input: bool,
    restraint_force: tuple[float, float, float] | None = None,
) -> Path | None:
    """Write the NAMD input, or remove a stale one for construction-only output."""
    path = directory / "run_drude_namd.conf"
    if not write_run_input:
        path.unlink(missing_ok=True)
        return None
    path.write_text(
        render_drude_namd_configuration(
            system.basis_vectors,
            lj_correction=lj_correction,
            restraint_force=restraint_force,
        ),
        encoding="utf-8",
    )
    return path


def _namd_restraint_lines(
    force_constants: tuple[float, float, float] | None,
) -> list[str]:
    """Render NAMD harmonic constraints with matched Cartesian stiffness."""
    if force_constants is None:
        return []
    coefficient = validate_namd_restraint_force(force_constants)
    assert coefficient is not None
    return [
        "constraints on",
        "consRef restraints_namd.pdb",
        "consKFile restraints_namd.pdb",
        "consKCol B",
        f"constraintScaling {coefficient:.12g}",
        "selectConstraints on",
        *(
            f"selectConstr{axis} {'on' if value > 0.0 else 'off'}"
            for axis, value in zip("XYZ", force_constants, strict=True)
        ),
    ]


def validate_namd_restraint_force(
    force_constants: tuple[float, float, float] | None,
) -> float | None:
    """Return the exact NAMD scale or reject an unrepresentable request."""
    if force_constants is None:
        return None
    active = [value for value in force_constants if value > 0.0]
    if not active:
        raise ValueError("position restraints require at least one positive force constant")
    if not np.allclose(active, active[0], atol=1.0e-9, rtol=0.0):
        raise ValueError(
            "NAMD position restraints require equal nonzero force constants on selected axes"
        )
    # NAMD uses U=k*d^2, while the public contract is U=1/2*K*d^2.
    # Convert K from kJ mol^-1 nm^-2 to kcal mol^-1 A^-2 and include 1/2.
    return active[0] * _KJ_NM2_TO_NAMD_KCAL_A2


def render_drude_restraint_pdb(
    system: DrudeSystem,
    atom_indices: tuple[int, ...],
    force_constants: tuple[float, float, float],
) -> str:
    """Encode selected real heavy atoms and the NAMD coefficient in the B column."""
    validate_namd_restraint_force(force_constants)
    selected = set(atom_indices)
    particle = 0
    rendered: list[str] = []
    for line in render_drude_pdb(system).splitlines():
        rendered_line = line
        if line.startswith("ATOM"):
            value = 1.0 if particle in selected else 0.0
            rendered_line = f"{line[:60]}{value:6.3f}{line[66:]}"
            particle += 1
        rendered.append(rendered_line)
    return "\n".join(rendered) + "\n"


def render_drude_namd_configuration(
    vectors: FloatArray,
    *,
    lj_correction: bool,
    restraint_force: tuple[float, float, float] | None = None,
) -> str:
    """Render a conservative minimization and 1 ps dual-thermostat smoke test.

    The Drude streams specify a 10--12 A potential switch (``vswitch``), not the
    additive CHARMM36 force switch.  ``lj_correction`` applies NAMD's isotropic
    long-range Lennard-Jones correction and is meaningful only for a
    homogeneous periodic crystal.  With ``drude on`` NAMD assumes the SWM4
    water model and refuses to start unless ``rigidBonds`` is ``water`` or
    ``all``; ``water`` leaves every cellulose bond flexible.
    """
    a, b, c = np.asarray(vectors, dtype=np.float64)
    center = 0.5 * (a + b + c)
    vector_lines = [
        f"cellBasisVector1 {a[0]:.9f} {a[1]:.9f} {a[2]:.9f}",
        f"cellBasisVector2 {b[0]:.9f} {b[1]:.9f} {b[2]:.9f}",
        f"cellBasisVector3 {c[0]:.9f} {c[1]:.9f} {c[2]:.9f}",
        f"cellOrigin {center[0]:.9f} {center[1]:.9f} {center[2]:.9f}",
    ]
    return "\n".join(
        [
            "structure system.psf",
            "coordinates structure_namd.pdb",
            "bincoordinates structure_namd.coor",
            f"parameters {DRUDE_MAIN_FILE}",
            f"parameters {DRUDE_CARBOHYDRATE_FILE}",
            "paraTypeCharmm on",
            "exclude scaled1-4",
            "1-4scaling 1.0",
            "switching on",
            "vdwForceSwitching off",
            "switchdist 10.0",
            "cutoff 12.0",
            "pairlistdist 14.0",
            f"LJcorrection {'on' if lj_correction else 'off'}",
            "PME yes",
            "PMEGridSpacing 1.0",
            "PMETolerance 1.0e-6",
            *vector_lines,
            "wrapAll off",
            "timestep 1.0",
            "rigidBonds water",
            "useSettle on",
            "langevin on",
            "langevinTemp 300",
            "langevinDamping 1.0",
            "langevinHydrogen on",
            "drude on",
            "drudeTemp 1.0",
            "drudeDamping 20.0",
            "drudeHardWall on",
            "drudeBondLen 0.25",
            "drudeNbTholeCut 5.0",
            *_namd_restraint_lines(restraint_force),
            "temperature 300",
            "outputEnergies 100",
            "outputName drude",
            "minimize 1000",
            "reinitvels 300",
            "run 1000",
            "",
        ]
    )


def render_drude_openmm_runner(
    vectors: FloatArray,
    *,
    lj_correction: bool,
    restraint_atom_indices: tuple[int, ...] = (),
    restraint_force: tuple[float, float, float] | None = None,
) -> str:
    """Render an OpenMM Drude minimization and 1 ps dynamics runner."""
    box = [[float(value) / 10.0 for value in row] for row in np.asarray(vectors)]
    return f'''#!/usr/bin/env python3
"""Run a finite-chain CelluloseForge CHARMM Drude system."""
import math
import struct

import openmm
from openmm import unit
from openmm.app import (
    CharmmParameterSet, CharmmPsfFile, DCDReporter, PME, Simulation, StateDataReporter
)

BOX_NM = {box!r}
RESTRAINT_ATOM_INDICES = {list(restraint_atom_indices)!r}
RESTRAINT_FORCE_KJ_MOL_NM2 = {restraint_force!r}
psf = CharmmPsfFile("system.psf")
parameters = CharmmParameterSet("{DRUDE_MAIN_FILE}", "{DRUDE_CARBOHYDRATE_FILE}")
with open("structure_namd.coor", "rb") as handle:
    count = struct.unpack("=i", handle.read(4))[0]
    values = struct.unpack(f"={{count * 3}}d", handle.read())
positions = [openmm.Vec3(*values[index:index+3]) for index in range(0, len(values), 3)]
positions = positions * unit.angstrom
if count != len(positions):
    raise SystemExit(f"binary coordinate count mismatch: {{count}} != {{len(positions)}}")

vectors = tuple(openmm.Vec3(*row) * unit.nanometer for row in BOX_NM)
lengths = [math.sqrt(sum(value * value for value in row)) for row in BOX_NM]
def angle(first, second):
    dot = sum(a*b for a, b in zip(first, second))
    norms = math.sqrt(sum(a*a for a in first) * sum(b*b for b in second))
    return math.degrees(math.acos(dot / norms))
psf.setBox(*(value * unit.nanometer for value in lengths),
           angle(BOX_NM[1], BOX_NM[2])*unit.degree,
           angle(BOX_NM[0], BOX_NM[2])*unit.degree,
           angle(BOX_NM[0], BOX_NM[1])*unit.degree)
system = psf.createSystem(
    parameters, nonbondedMethod=PME, nonbondedCutoff=1.2*unit.nanometer,
    switchDistance=1.0*unit.nanometer, constraints=None,
    ewaldErrorTolerance=1.0e-6,
)
system.setDefaultPeriodicBoxVectors(*vectors)
if system.getNumParticles() != count:
    raise SystemExit(
        f"PSF particle count mismatch: {{system.getNumParticles()}} != {{count}}"
    )
for force in system.getForces():
    if isinstance(force, openmm.NonbondedForce):
        force.setUseDispersionCorrection({lj_correction!r})
if RESTRAINT_FORCE_KJ_MOL_NM2 is not None:
    restraint = openmm.CustomExternalForce(
        "0.5*(kx*(x-x0)^2+ky*(y-y0)^2+kz*(z-z0)^2)"
    )
    for name in ("x0", "y0", "z0", "kx", "ky", "kz"):
        restraint.addPerParticleParameter(name)
    kx, ky, kz = RESTRAINT_FORCE_KJ_MOL_NM2
    for index in RESTRAINT_ATOM_INDICES:
        x0, y0, z0 = positions[index].value_in_unit(unit.nanometer)
        restraint.addParticle(index, [x0, y0, z0, kx, ky, kz])
    restraint.setName("CellulosePositionRestraints")
    system.addForce(restraint)
integrator = openmm.DrudeLangevinIntegrator(
    300*unit.kelvin, 1/unit.picosecond, 1*unit.kelvin,
    20/unit.picosecond, 0.001*unit.picoseconds,
)
integrator.setMaxDrudeDistance(0.025*unit.nanometer)
simulation = Simulation(psf.topology, system, integrator)
simulation.context.setPeriodicBoxVectors(*vectors)
simulation.context.setPositions(positions)
simulation.context.computeVirtualSites()
print("Initial potential energy:", simulation.context.getState(getEnergy=True).getPotentialEnergy())
simulation.minimizeEnergy(maxIterations=1000)
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
simulation.reporters.append(StateDataReporter(
    "openmm_drude.log", 100, step=True, potentialEnergy=True, temperature=True, speed=True
))
simulation.reporters.append(DCDReporter("openmm_drude.dcd", 100))
simulation.step(1000)
simulation.saveCheckpoint("openmm_drude.chk")
'''


def write_drude_namd_system(  # noqa: PLR0913 - explicit build options
    crystal: BulkCrystal,
    toppar_directory: str | Path,
    output_directory: str | Path,
    *,
    chain_boundary: str,
    lj_correction: bool,
    write_run_input: bool = True,
    write_openmm_runner: bool = False,
    restraint_chain_ids: tuple[int, ...] = (),
    restraint_force: tuple[float, float, float] | None = None,
) -> DrudeSystemFiles:
    """Write a native CHARMM Drude PSF and coordinates, and optionally a NAMD input.

    Covalently periodic chains are construction-only: three lone pairs per
    chain have a host across the periodic boundary. NAMD 3.0.3 does not image
    those hosts correctly and stops before energy evaluation, so no run input
    may be written for them.
    """
    if chain_boundary == "periodic" and (write_run_input or write_openmm_runner):
        raise ValueError(
            "periodic CHARMM Drude chains are construction-only; no NAMD or OpenMM run input "
            "is written because neither engine images cross-boundary lone-pair hosts correctly"
        )
    if restraint_force is not None:
        if not any(value > 0.0 for value in restraint_force):
            raise ValueError("position restraints require at least one positive force constant")
        if write_run_input:
            validate_namd_restraint_force(restraint_force)
    main_source, carbohydrate_source = validate_drude_toppar_directory(toppar_directory)
    checksums = verify_audited_streams(main_source, carbohydrate_source)
    definitions = read_drude_cellulose_definitions(main_source, carbohydrate_source)
    chain_residues: dict[int, set[int]] = {}
    for chain, residue in zip(crystal.chain_numbers, crystal.residue_numbers, strict=True):
        chain_residues.setdefault(chain, set()).add(residue)
    lengths = {len(residues) for residues in chain_residues.values()}
    if len(lengths) != 1:
        raise ValueError("CHARMM Drude output requires chains of one degree of polymerization")
    topology = build_drude_cellulose_chain(
        definitions, lengths.pop(), periodic=chain_boundary == "periodic"
    )
    system = build_drude_system(crystal, topology, chain_boundary)
    restraint_indices = drude_heavy_atom_restraint_indices(system, restraint_chain_ids)
    if restraint_force is not None and not restraint_indices:
        raise ValueError("CHARMM Drude position restraints selected no real heavy atoms")
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=True)
    main_stream = directory / DRUDE_MAIN_FILE
    carbohydrate_stream = directory / DRUDE_CARBOHYDRATE_FILE
    shutil.copy2(main_source, main_stream)
    shutil.copy2(carbohydrate_source, carbohydrate_stream)
    psf = directory / "system.psf"
    pdb = directory / "structure_namd.pdb"
    binary_coordinates = directory / "structure_namd.coor"
    psf.write_text(render_drude_psf(system), encoding="utf-8")
    pdb.write_text(render_drude_pdb(system), encoding="utf-8")
    write_namd_binary_coordinates(binary_coordinates, system.cartesian_coordinates)
    configuration = _write_run_configuration(
        directory,
        system,
        lj_correction=lj_correction,
        write_run_input=write_run_input,
        restraint_force=restraint_force,
    )
    restraint_reference = directory / "restraints_namd.pdb"
    if write_run_input and restraint_force is not None:
        restraint_reference.write_text(
            render_drude_restraint_pdb(system, restraint_indices, restraint_force),
            encoding="utf-8",
        )
    else:
        restraint_reference.unlink(missing_ok=True)
    openmm_runner = directory / "run_drude_openmm.py"
    if write_openmm_runner:
        openmm_runner.write_text(
            render_drude_openmm_runner(
                system.basis_vectors,
                lj_correction=lj_correction,
                restraint_atom_indices=restraint_indices,
                restraint_force=restraint_force,
            ),
            encoding="utf-8",
        )
    else:
        openmm_runner.unlink(missing_ok=True)
    chains = system.chain_count
    return DrudeSystemFiles(
        psf=psf,
        pdb=pdb,
        binary_coordinates=binary_coordinates,
        configuration=configuration,
        openmm_runner=openmm_runner if write_openmm_runner else None,
        restraint_reference=(
            restraint_reference if write_run_input and restraint_force is not None else None
        ),
        main_stream=main_stream,
        carbohydrate_stream=carbohydrate_stream,
        source_sha256=checksums,
        counts={
            "particles": system.atom_count,
            "real_atoms": chains * sum(site.kind == "atom" for site in topology.sites),
            "drude_particles": chains * len(topology.drude_parents),
            "lone_pairs": chains * len(topology.lone_pairs),
            "anisotropy_records": chains * len(topology.anisotropies),
            "bonds": chains * len(topology.bonds),
            "angles": chains * len(topology.angles),
            "dihedrals": chains * len(topology.dihedrals),
        },
    )
