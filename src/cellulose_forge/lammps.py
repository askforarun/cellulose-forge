"""Native LAMMPS data and input files for mapped cellulose systems.

The optional OpenMM interpreter is used as a force-field parser.  It reads the
same temporary GROMACS-format topology used by the independently validated
OpenMM backend; GROMACS itself is not required.  The converter deliberately
accepts only force terms and 1--4 semantics that can be represented exactly by
the emitted LAMMPS styles.
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

from cellulose_forge.glycam_mapping import GlycamMappedSystem

LAMMPS_DATA_FILE = "system.data"
LAMMPS_INPUT_FILE = "in.cellulose"


@dataclass(frozen=True)
class LammpsSystemFiles:
    """Files written for a LAMMPS simulation."""

    data: Path
    input: Path


_CONVERTER = r"""
import json
import math
import sys

import numpy
import openmm
from openmm import unit
from openmm.app import GromacsGroFile, GromacsTopFile, PME

request = json.loads(open(sys.argv[1], encoding="utf-8").read())
gro = GromacsGroFile(request["coordinates"])
box_vectors = gro.getPeriodicBoxVectors()
top = GromacsTopFile(
    request["topology"], periodicBoxVectors=box_vectors,
    includeDir=request["include_directory"],
)
system = top.createSystem(
    nonbondedMethod=PME, nonbondedCutoff=1.0*unit.nanometer, constraints=None,
)
if system.getNumParticles() != len(request["molecule_ids"]):
    raise SystemExit("LAMMPS conversion atom count does not match mapped coordinates")

nonbonded = next(
    force for force in system.getForces() if isinstance(force, openmm.NonbondedForce)
)
particles = [nonbonded.getParticleParameters(i) for i in range(system.getNumParticles())]

def type_id(table, value):
    key = tuple(round(float(item), 12) for item in value)
    if key not in table:
        table[key] = len(table) + 1
    return table[key]

atom_types = {}
atom_type_ids = []
for index, (charge, sigma, epsilon) in enumerate(particles):
    mass = system.getParticleMass(index).value_in_unit(unit.dalton)
    atom_type_ids.append(type_id(atom_types, (
        mass,
        epsilon.value_in_unit(unit.kilocalorie_per_mole),
        sigma.value_in_unit(unit.angstrom),
    )))

charmm = request["force_field"] == "charmm36"
exception_pairs = set()
exception_by_type = {}
for index in range(nonbonded.getNumExceptions()):
    first, second, charge, sigma, epsilon = nonbonded.getExceptionParameters(index)
    charge_value = charge.value_in_unit(unit.elementary_charge**2)
    epsilon_value = epsilon.value_in_unit(unit.kilocalorie_per_mole)
    if abs(charge_value) < 1e-12 and abs(epsilon_value) < 1e-12:
        continue
    qi, si, ei = particles[first]
    qj, sj, ej = particles[second]
    expected_charge = (qi*qj).value_in_unit(unit.elementary_charge**2)
    if abs(charge_value-expected_charge) > 1e-6*max(1.0, abs(expected_charge)):
        raise SystemExit(f"unsupported scaled Coulomb 1-4 exception {first+1}-{second+1}")
    expected_sigma = 0.5*(si+sj).value_in_unit(unit.angstrom)
    expected_epsilon = math.sqrt(
        ei.value_in_unit(unit.kilocalorie_per_mole)
        * ej.value_in_unit(unit.kilocalorie_per_mole)
    )
    sigma_value = sigma.value_in_unit(unit.angstrom)
    if not charmm and (
        abs(sigma_value-expected_sigma) > 1e-6*max(1.0, abs(expected_sigma))
        or abs(epsilon_value-expected_epsilon) > 1e-6*max(1.0, abs(expected_epsilon))
    ):
        raise SystemExit(
            f"1-4 exception {first+1}-{second+1} requires pair-specific parameters "
            "that the GLYCAM LAMMPS representation cannot express"
        )
    exception_pairs.add(tuple(sorted((int(first), int(second)))))
    type_pair = tuple(sorted((atom_type_ids[first], atom_type_ids[second])))
    value = (epsilon_value, sigma_value)
    previous = exception_by_type.get(type_pair)
    if previous is not None and any(
        abs(a-b) > 1e-6*max(1.0, abs(b)) for a, b in zip(previous, value)
    ):
        raise SystemExit(f"inconsistent CHARMM 1-4 parameters for atom types {type_pair}")
    exception_by_type[type_pair] = value

bonds, angles, dihedrals = [], [], []
bond_types, angle_types, dihedral_types = {}, {}, {}
covalent_pairs = {
    tuple(sorted((first.index, second.index))) for first, second in top.topology.bonds()
}
if not charmm:
    adjacency = [set() for _ in range(system.getNumParticles())]
    for first, second in covalent_pairs:
        adjacency[first].add(second)
        adjacency[second].add(first)
    graph_one_four = set()
    for root in range(system.getNumParticles()):
        distance = {root: 0}
        frontier = [root]
        while frontier:
            first = frontier.pop(0)
            if distance[first] == 3:
                continue
            for second in adjacency[first]:
                if second not in distance:
                    distance[second] = distance[first] + 1
                    frontier.append(second)
        graph_one_four.update(
            (root, second) for second, depth in distance.items()
            if root < second and depth == 3
        )
    if graph_one_four != exception_pairs:
        missing = sorted(graph_one_four-exception_pairs)
        extra = sorted(exception_pairs-graph_one_four)
        raise SystemExit(
            "GLYCAM topology 1-4 exceptions do not match the covalent graph; "
            f"missing={missing[:5]}, extra={extra[:5]}"
        )
raw_bonds = []
raw_angles = []
raw_dihedrals = []
raw_impropers = []
# OpenMM's GromacsTopFile writes GROMACS function-2 (harmonic) impropers as this
# CustomTorsionForce. Only this exact form is converted.
harmonic_torsion_expression = (
    "0.5*k*(thetap-theta0)^2; thetap = step(-(theta-theta0+pi))*2*pi+theta"
    "+step(theta-theta0-pi)*(-2*pi); pi = %.15g" % math.pi
)
for force in system.getForces():
    if isinstance(force, openmm.HarmonicBondForce):
        for i in range(force.getNumBonds()):
            raw_bonds.append(force.getBondParameters(i))
    elif isinstance(force, openmm.HarmonicAngleForce):
        for i in range(force.getNumAngles()):
            raw_angles.append(force.getAngleParameters(i))
    elif isinstance(force, openmm.PeriodicTorsionForce):
        for i in range(force.getNumTorsions()):
            raw_dihedrals.append(force.getTorsionParameters(i))
    elif isinstance(force, openmm.CustomTorsionForce):
        names = [
            force.getPerTorsionParameterName(i)
            for i in range(force.getNumPerTorsionParameters())
        ]
        if (
            force.getEnergyFunction() != harmonic_torsion_expression
            or names != ["theta0", "k"]
            or force.getNumGlobalParameters() != 0
        ):
            raise SystemExit(
                "unsupported CustomTorsionForce for LAMMPS conversion: "
                f"{force.getEnergyFunction()}"
            )
        for i in range(force.getNumTorsions()):
            raw_impropers.append(force.getTorsionParameters(i))
    elif isinstance(force, (openmm.CMMotionRemover, openmm.NonbondedForce)):
        pass
    else:
        raise SystemExit(f"unsupported OpenMM force for LAMMPS conversion: {type(force).__name__}")

urey_bradley = {}
for a, b, length, k in raw_bonds:
    pair = tuple(sorted((int(a), int(b))))
    params = (
        0.5*k.value_in_unit(unit.kilocalorie_per_mole/unit.angstrom**2),
        length.value_in_unit(unit.angstrom),
    )
    if charmm and pair not in covalent_pairs:
        urey_bradley[pair] = params
    else:
        bonds.append((type_id(bond_types, params), int(a)+1, int(b)+1))

for a, b, c, theta, k in raw_angles:
    params = [
        0.5*k.value_in_unit(unit.kilocalorie_per_mole/unit.radian**2),
        theta.value_in_unit(unit.degree),
    ]
    if charmm:
        params.extend(urey_bradley.get(tuple(sorted((int(a), int(c)))), (0.0, 0.0)))
    parameters = tuple(params)
    angles.append((type_id(angle_types, parameters), int(a)+1, int(b)+1, int(c)+1))

weighted_pairs = set()
for a, b, c, d, periodicity, phase, k in raw_dihedrals:
    endpoint_pair = tuple(sorted((int(a), int(d))))
    weight = 0.0
    if charmm and endpoint_pair in exception_pairs and endpoint_pair not in weighted_pairs:
        weight = 1.0
        weighted_pairs.add(endpoint_pair)
    params = [
        k.value_in_unit(unit.kilocalorie_per_mole),
        int(periodicity),
        phase.value_in_unit(unit.degree),
    ]
    phase_degrees = params[2]
    if abs(phase_degrees-round(phase_degrees)) > 1e-8:
        raise SystemExit(
            f"torsion phase {phase_degrees} degrees is not an integer "
            "and cannot be represented by the LAMMPS CHARMM dihedral styles"
        )
    # Both CHARMM styles use K, n, d, weight coefficients.  GLYCAM's 1-4
    # interactions are supplied by special_bonds, so a zero dihedral weight
    # avoids double-counting while requiring only standard LAMMPS packages.
    params.append(weight)
    parameters = tuple(params)
    dihedrals.append((
        type_id(dihedral_types, parameters), int(a)+1, int(b)+1, int(c)+1, int(d)+1
    ))
# CHARMM dihedral types with zero force constant (for example HCP1-NC2D1-
# CC3161-CC3162 in BGLCNA) are dropped by OpenMM, but their 1-4 pairs remain.
# Carry each such pair on a zero-force-constant dihedral with 1-4 weight 1.
zero_dihedral_pairs = []
if charmm and weighted_pairs != exception_pairs:
    adjacency = [set() for _ in range(system.getNumParticles())]
    for first, second in covalent_pairs:
        adjacency[first].add(second)
        adjacency[second].add(first)
    for first, last in sorted(exception_pairs - weighted_pairs):
        path = next(
            (
                (first, second, third, last)
                for second in sorted(adjacency[first])
                for third in sorted(adjacency[second] & adjacency[last])
                if len({first, second, third, last}) == 4
            ),
            None,
        )
        if path is None:
            continue
        parameters = (0.0, 1, 0.0, 1.0)
        dihedrals.append((type_id(dihedral_types, parameters), *(index + 1 for index in path)))
        weighted_pairs.add((first, last))
        zero_dihedral_pairs.append((first, last))

impropers, improper_types = [], {}
for a, b, c, d, (theta0, k) in raw_impropers:
    if abs(theta0) > 1e-12:
        raise SystemExit(
            f"improper {a+1}-{b+1}-{c+1}-{d+1} has a nonzero equilibrium angle; "
            "LAMMPS improper_style harmonic is used only for planar (chi0 = 0) impropers"
        )
    # GROMACS/OpenMM: E = (1/2) k xi^2 (kJ/mol/rad^2); LAMMPS harmonic: E = K chi^2.
    # Hence K = k/2, then kJ -> kcal.
    parameters = (0.5*k/4.184, 0.0)
    impropers.append((
        type_id(improper_types, parameters), int(a)+1, int(b)+1, int(c)+1, int(d)+1
    ))

if charmm and weighted_pairs != exception_pairs:
    missing = sorted(exception_pairs-weighted_pairs)
    raise SystemExit(f"CHARMM 1-4 exceptions without a matching dihedral: {missing[:5]}")

positions = gro.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
box = [[float(x) for x in v.value_in_unit(unit.angstrom)] for v in box_vectors]
a, b, c = box
if any(abs(value) > 1e-6 for value in (a[1], a[2], b[2])):
    raise SystemExit("LAMMPS output requires a restricted triclinic box")

# Wrap atoms and derive consistent LAMMPS image flags from the actual bond
# graph.  This keeps covalent bonds crossing a periodic face contiguous when
# coordinates are unwrapped, without relying on LAMMPS reset_atoms heuristics.
box_matrix = numpy.asarray(box)
fractional = numpy.asarray(positions) @ numpy.linalg.inv(box_matrix)
wrapped_fractional = fractional - numpy.floor(fractional)
positions = wrapped_fractional @ box_matrix
neighbors = [[] for _ in range(system.getNumParticles())]
for _, first, second in bonds:
    neighbors[first-1].append(second-1)
    neighbors[second-1].append(first-1)
images = numpy.zeros((system.getNumParticles(), 3), dtype=int)
visited = numpy.zeros(system.getNumParticles(), dtype=bool)
for root in range(system.getNumParticles()):
    if visited[root]:
        continue
    visited[root] = True
    stack = [root]
    while stack:
        first = stack.pop()
        for second in neighbors[first]:
            candidate = images[first] - numpy.rint(
                wrapped_fractional[second] - wrapped_fractional[first]
            ).astype(int)
            if not visited[second]:
                images[second] = candidate
                visited[second] = True
                stack.append(second)
            elif not numpy.array_equal(images[second], candidate):
                # A periodic polymer graph can contain a cycle whose net image
                # displacement is nonzero. LAMMPS uses the individual bond's
                # minimum image for forces; retain the first spanning-tree flag.
                pass

lines = [
    "CelluloseForge LAMMPS data", "",
    f"{system.getNumParticles()} atoms", f"{len(bonds)} bonds",
    f"{len(angles)} angles", f"{len(dihedrals)} dihedrals", f"{len(impropers)} impropers", "",
    f"{len(atom_types)} atom types", f"{len(bond_types)} bond types",
    f"{len(angle_types)} angle types", f"{len(dihedral_types)} dihedral types",
    *([f"{len(improper_types)} improper types"] if improper_types else []), "",
    f"0.0 {a[0]:.12g} xlo xhi", f"0.0 {b[1]:.12g} ylo yhi",
    f"0.0 {c[2]:.12g} zlo zhi", f"{b[0]:.12g} {c[0]:.12g} {c[1]:.12g} xy xz yz",
    "", "Masses", "",
]
for params, index in sorted(atom_types.items(), key=lambda item: item[1]):
    lines.append(f"{index} {params[0]:.12g}")
if charmm:
    lines += ["", "PairIJ Coeffs", ""]
    types_by_id = {index: params for params, index in atom_types.items()}
    for first in range(1, len(atom_types)+1):
        for second in range(first, len(atom_types)+1):
            left, right = types_by_id[first], types_by_id[second]
            epsilon = math.sqrt(left[1]*right[1])
            sigma = 0.5*(left[2]+right[2])
            epsilon14, sigma14 = exception_by_type.get(
                (first, second), (epsilon, sigma)
            )
            lines.append(
                f"{first} {second} {epsilon:.12g} {sigma:.12g} "
                f"{epsilon14:.12g} {sigma14:.12g}"
            )
else:
    lines += ["", "Pair Coeffs", ""]
    for params, index in sorted(atom_types.items(), key=lambda item: item[1]):
        lines.append(f"{index} {params[1]:.12g} {params[2]:.12g}")
lines += ["", "Bond Coeffs", ""]
for params, index in sorted(bond_types.items(), key=lambda item: item[1]):
    lines.append(f"{index} {params[0]:.12g} {params[1]:.12g}")
lines += ["", "Angle Coeffs", ""]
for params, index in sorted(angle_types.items(), key=lambda item: item[1]):
    lines.append(f"{index} " + " ".join(f"{value:.12g}" for value in params))
lines += ["", "Dihedral Coeffs", ""]
for params, index in sorted(dihedral_types.items(), key=lambda item: item[1]):
    lines.append(
        f"{index} {params[0]:.12g} {int(params[1])} "
        f"{int(params[2])} {params[3]:.1f}"
    )
if improper_types:
    lines += ["", "Improper Coeffs", ""]
    for params, index in sorted(improper_types.items(), key=lambda item: item[1]):
        lines.append(f"{index} {params[0]:.12g} {params[1]:.12g}")
lines += ["", "Atoms # full", ""]
for index, ((charge, _, _), atom_type, molecule, xyz, image) in enumerate(zip(
    particles, atom_type_ids, request["molecule_ids"], positions, images
), start=1):
    q = charge.value_in_unit(unit.elementary_charge)
    lines.append(
        f"{index} {molecule} {atom_type} {q:.12g} "
        f"{xyz[0]:.12g} {xyz[1]:.12g} {xyz[2]:.12g} "
        f"{image[0]} {image[1]} {image[2]}"
    )
sections = [("Bonds", bonds), ("Angles", angles), ("Dihedrals", dihedrals)]
if impropers:
    sections.append(("Impropers", impropers))
for title, records in sections:
    lines += ["", title, ""]
    for index, record in enumerate(records, start=1):
        lines.append(f"{index} " + " ".join(str(value) for value in record))
open(request["data_path"], "w", encoding="utf-8").write("\n".join(lines) + "\n")
print(json.dumps({
    "atom_count": system.getNumParticles(), "bond_count": len(bonds),
    "angle_count": len(angles), "dihedral_count": len(dihedrals),
    "improper_count": len(impropers),
    "zero_force_constant_one_four_dihedrals": len(zero_dihedral_pairs),
    "atom_type_count": len(atom_types),
    "special_bonds": "charmm" if charmm else "lj/coul 0 0 1",
    "urey_bradley_count": len(urey_bradley),
    "weighted_one_four_count": len(weighted_pairs),
}))
"""


def _resolve_python(executable: str | None) -> str:
    return executable or sys.executable


def write_lammps_system(  # noqa: PLR0913
    system: GlycamMappedSystem,
    gromacs_topology: str | Path,
    gromacs_coordinates: str | Path,
    output_directory: str | Path,
    *,
    force_field: str = "glycam",
    dispersion_correction: bool = False,
    python_executable: str | None = None,
) -> tuple[LammpsSystemFiles, dict[str, object]]:
    """Write an exact LAMMPS representation or reject unsupported semantics."""
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=True)
    files = LammpsSystemFiles(directory / LAMMPS_DATA_FILE, directory / LAMMPS_INPUT_FILE)
    molecule_ids = [int(value) for value in system.chain_numbers]
    with tempfile.TemporaryDirectory(prefix="cellulose_forge_lammps_") as temporary:
        request_path = Path(temporary) / "request.json"
        request_path.write_text(
            json.dumps(
                {
                    "topology": str(Path(gromacs_topology).resolve()),
                    "coordinates": str(Path(gromacs_coordinates).resolve()),
                    "include_directory": str(Path(gromacs_topology).resolve().parent),
                    "data_path": str(files.data.resolve()),
                    "molecule_ids": molecule_ids,
                    "force_field": force_field.lower(),
                }
            ),
            encoding="utf-8",
        )
        result = subprocess.run(
            [_resolve_python(python_executable), "-c", _CONVERTER, str(request_path)],
            check=False,
            capture_output=True,
            text=True,
        )
    if result.returncode != 0:
        raise RuntimeError(f"LAMMPS conversion failed:\n{result.stdout}\n{result.stderr}")
    if not files.data.is_file():
        raise RuntimeError("LAMMPS conversion did not create system.data")
    metadata = json.loads(result.stdout.strip().splitlines()[-1])
    if force_field.lower() == "charmm36":
        angle_style = "charmm"
        dihedral_style = "charmmfsw"
        pair_style = "lj/charmmfsw/coul/long 10.0 12.0"
        special_bonds = "charmm"
    else:
        angle_style = "harmonic"
        dihedral_style = "charmm"
        pair_style = "lj/cut/coul/long 10.0"
        special_bonds = "lj/coul 0.0 0.0 1.0"
    has_impropers = int(metadata.get("improper_count", 0)) > 0
    improper_line = "improper_style harmonic\n" if has_impropers else ""
    improper_thermo = " eimp" if has_impropers else ""
    tail_correction = "yes" if dispersion_correction and force_field.lower() == "glycam" else "no"
    files.input.write_text(
        f"""units real
atom_style full
boundary p p p
bond_style harmonic
angle_style {angle_style}
dihedral_style {dihedral_style}
{improper_line}pair_style {pair_style}
pair_modify mix arithmetic tail {tail_correction}
special_bonds {special_bonds}
read_data system.data
kspace_style pppm 1.0e-5
thermo 100
thermo_style custom step temp pe ebond eangle edihed evdwl ecoul elong{improper_thermo}
timestep 1.0
# This file is an initial-energy validation input, not a dynamics protocol.
run 0
# Explicit termination avoids an exit-time fault in some packaged LAMMPS builds.
quit 0
""",
        encoding="utf-8",
    )
    metadata.update(
        {
            "data": files.data.name,
            "input": files.input.name,
            "dispersion_correction": tail_correction == "yes",
            "input_scope": "single-point run 0 validation; not a dynamics protocol",
        }
    )
    return files, metadata
