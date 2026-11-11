"""Compare chitin CHARMM36 energies term by term across GROMACS, OpenMM and LAMMPS.

For alpha and beta chitin (1 0 0) slabs with finite and periodic chains, the
same generated GROMACS topology is evaluated at the built coordinates by:

- GROMACS (reference);
- OpenMM, through the CelluloseForge OpenMM reader (``GromacsTopFile``);
- LAMMPS, through the CelluloseForge LAMMPS converter.

Compared categories: bonds, angles, Urey-Bradley, proper dihedrals, the two
acetyl impropers, 1-4 Lennard-Jones, 1-4 electrostatics, and total nonbonded
energy, plus the total potential.

GROMACS reports angles and Urey-Bradley together as ``U-B``. The split is taken
from GROMACS's own resolved parameters (``gmx dump``) evaluated at the GRO
coordinates; their sum must reproduce the GROMACS ``U-B`` energy. LAMMPS
reports angle + Urey-Bradley as ``eangle`` and tallies CHARMM 1-4 terms into
``evdwl``/``ecoul``; a second LAMMPS run with Urey-Bradley constants and
dihedral 1-4 weights set to zero separates them.

OpenMM is compared with GROMACS using a plain 1.0 nm cutoff; LAMMPS uses
``lj/charmmfsw/coul/long`` (10-12 angstrom force switch), so GROMACS is rerun
with the matching force-switch settings for that comparison.

Usage:
    python validation/validate_chitin_charmm36_engines.py --work-directory DIR \\
        --openmm-python PYTHON [--lammps lmp_serial] [--output FILE]
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import numpy.typing as npt

sys.path.insert(0, str(Path(__file__).resolve().parent))

from validate_chitin_charmm36_gromacs import _build, _energies, _gmx

from cellulose_forge.charmm36 import (
    FiniteCharmm36Topology,
    PeriodicCharmm36Topology,
    map_finite_charmm36,
    map_periodic_charmm36,
    native_finite_charmm36_topology,
    native_periodic_charmm36_topology,
    write_charmm36_system,
)
from cellulose_forge.forcefields import discover_charmm36_directory
from cellulose_forge.lammps import write_lammps_system

_KCAL = 4.184
# LAMMPS switches to the CHARMM Coulomb constant (332.0716 kcal A/mol/e^2);
# GROMACS uses 138.935458 kJ nm/mol/e^2 = 332.0637 kcal A/mol/e^2.
_LAMMPS_CHARMM_COULOMB = 332.0716
_GROMACS_COULOMB = 138.935458 * 10.0 / _KCAL
_COULOMB_RATIO = _GROMACS_COULOMB / _LAMMPS_CHARMM_COULOMB
_LAMMPS_SWITCH_INNER = 10.0
_LAMMPS_SWITCH_OUTER = 12.0
# Bonded terms use identical parameters and coordinates: agreement is limited by
# single-precision GROMACS output. Nonbonded totals differ by PME/PPPM accuracy.
_BONDED_RELATIVE_TOLERANCE = 2.0e-5
_BONDED_ABSOLUTE_TOLERANCE = 0.05
# After removing documented LAMMPS energy conventions, nonbonded totals must
# agree to this fraction of the summed nonbonded component magnitudes.
_NONBONDED_RELATIVE_TOLERANCE = 5.0e-5
_DYNAMICS_STEPS = 2000
_EXPECTED_TERMS = frozenset(
    {
        "bond",
        "angle",
        "urey_bradley",
        "proper_dihedral",
        "improper",
        "lj_14",
        "coulomb_14",
        "nonbonded",
        "potential",
    }
)
_GMX_TERMS = (
    "Bond",
    "U-B",
    "Proper-Dih.",
    "Improper-Dih.",
    "LJ-14",
    "Coulomb-14",
    "LJ-(SR)",
    "Coulomb-(SR)",
    "Coul.-recip.",
    "Potential",
)
_PLAIN_MDP = """\
integrator = md
nsteps = 0
constraints = none
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 1.0
vdwtype = cutoff
rvdw = 1.0
vdw-modifier = None
coulomb-modifier = None
ewald-rtol = 1e-6
fourierspacing = 0.08
DispCorr = no
pbc = xyz
nstenergy = 1
nstcalcenergy = 1
periodic-molecules = {periodic}
"""
_SWITCHED_MDP = """\
integrator = md
nsteps = 0
constraints = none
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 1.2
vdwtype = cutoff
rvdw = 1.2
vdw-modifier = force-switch
rvdw-switch = 1.0
coulomb-modifier = None
ewald-rtol = 1e-6
fourierspacing = 0.08
DispCorr = no
pbc = xyz
nstenergy = 1
nstcalcenergy = 1
periodic-molecules = {periodic}
"""
_OPENMM_EVALUATOR = r"""
import json
import sys

import openmm
from openmm import unit
from openmm.app import PME, GromacsGroFile, GromacsTopFile

request = json.loads(sys.argv[1])
gro = GromacsGroFile(request["coordinates"])
box = gro.getPeriodicBoxVectors()
top = GromacsTopFile(request["topology"], periodicBoxVectors=box,
                     includeDir=request["include_directory"])
system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*unit.nanometer,
                          constraints=None, ewaldErrorTolerance=1e-6)
covalent = {tuple(sorted((a.index, b.index))) for a, b in top.topology.bonds()}
periodic = request["periodic"]
groups = {}
for force in list(system.getForces()):
    if isinstance(force, openmm.HarmonicBondForce):
        bond = openmm.HarmonicBondForce()
        urey = openmm.HarmonicBondForce()
        for i in range(force.getNumBonds()):
            a, b, length, k = force.getBondParameters(i)
            target = bond if tuple(sorted((a, b))) in covalent else urey
            target.addBond(a, b, length, k)
        for created, name in ((bond, "bond"), (urey, "urey_bradley")):
            created.setUsesPeriodicBoundaryConditions(periodic)
            created.setForceGroup(len(groups) + 1)
            groups[name] = created.getForceGroup()
            system.addForce(created)
        force.setForceGroup(31)
        force.setName("replaced")
for index, force in enumerate(system.getForces()):
    if force.getName() == "replaced":
        system.removeForce(index)
        break
named = {
    openmm.HarmonicAngleForce: "angle",
    openmm.PeriodicTorsionForce: "proper_dihedral",
    openmm.CustomTorsionForce: "improper",
    openmm.NonbondedForce: "nonbonded",
}
for force in system.getForces():
    for kind, name in named.items():
        if isinstance(force, kind):
            force.setForceGroup(len(groups) + 1)
            groups[name] = force.getForceGroup()
            if isinstance(force, openmm.NonbondedForce):
                force.setUseDispersionCorrection(False)
                force.setExceptionsUsePeriodicBoundaryConditions(periodic)
            elif periodic:
                force.setUsesPeriodicBoundaryConditions(True)
context = openmm.Context(system, openmm.VerletIntegrator(0.001),
                         openmm.Platform.getPlatformByName("Reference"))
context.setPeriodicBoxVectors(*box)
context.setPositions(gro.positions)

def energy(group):
    state = context.getState(getEnergy=True, groups={group})
    return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

result = {name: energy(group) for name, group in groups.items()}
result["potential"] = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
    unit.kilojoule_per_mole
)
nonbonded = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
particles = [nonbonded.getParticleParameters(i) for i in range(nonbonded.getNumParticles())]
exceptions = [nonbonded.getExceptionParameters(i) for i in range(nonbonded.getNumExceptions())]
for i in range(nonbonded.getNumParticles()):
    nonbonded.setParticleParameters(i, 0.0, 1.0, 0.0)
nonbonded.updateParametersInContext(context)
one_four = energy(groups["nonbonded"])
for i, (a, b, q, sigma, epsilon) in enumerate(exceptions):
    nonbonded.setExceptionParameters(i, a, b, 0.0, sigma, epsilon)
nonbonded.updateParametersInContext(context)
result["lj_14"] = energy(groups["nonbonded"])
result["coulomb_14"] = one_four - result["lj_14"]
result["improper_count"] = sum(
    f.getNumTorsions() for f in system.getForces() if isinstance(f, openmm.CustomTorsionForce)
)
print(json.dumps(result))
"""


def _gromacs(directory: Path, template: str, name: str, *, periodic: bool) -> dict[str, float]:
    mdp = directory / f"{name}.mdp"
    mdp.write_text(template.format(periodic="yes" if periodic else "no"), encoding="utf-8")
    _gmx(
        "grompp", "-f", mdp.name, "-c", "structure.gro", "-p", "system.top",
        "-o", f"{name}.tpr", "-po", f"{name}.mdout.mdp", "-maxwarn", "0",
        cwd=directory,
    )  # fmt: skip
    _gmx(
        "mdrun", "-s", f"{name}.tpr", "-deffnm", name, "-ntmpi", "1", "-ntomp", "1",
        cwd=directory,
    )  # fmt: skip
    values = _energies(directory, f"{name}.edr")
    return {term: values[term] for term in _GMX_TERMS if term in values}


def _gromacs_angle_split(directory: Path, tpr: str) -> dict[str, float]:
    """Evaluate GROMACS's resolved angle and Urey-Bradley terms separately."""
    dump = _gmx("dump", "-s", tpr, cwd=directory)
    parameters = {
        int(index): (float(theta), float(ktheta), float(r13), float(kub))
        for index, theta, ktheta, r13, kub in re.findall(
            r"functype\[(\d+)\]=UREY_BRADLEY, thetaA=\s*(\S+), kthetaA=\s*(\S+), "
            r"r13A=\s*(\S+), kUBA=\s*(\S+),",
            dump,
        )
    }
    rest = dump.split("\n      U-B:", 1)[1]
    end = re.search(r"\n {6}\S", rest)
    entries = re.findall(
        r"^\s+\d+ type=(\d+) \(UREY_BRADLEY\)\s+(\d+)\s+(\d+)\s+(\d+)",
        rest[: end.start()] if end else rest,
        flags=re.MULTILINE,
    )
    molecules = sum(int(value) for value in re.findall(r"#molecules\s+=\s+(\d+)", dump))
    lines = (directory / "structure.gro").read_text(encoding="utf-8").splitlines()
    count = int(lines[1])
    xyz = np.asarray(
        [
            [float(line[20:28]), float(line[28:36]), float(line[36:44])]
            for line in lines[2 : 2 + count]
        ]
    )
    box_values = [float(value) for value in lines[2 + count].split()] + [0.0] * 6
    box = np.asarray(
        [
            [box_values[0], box_values[3], box_values[4]],
            [box_values[5], box_values[1], box_values[6]],
            [box_values[7], box_values[8], box_values[2]],
        ]
    )
    inverse = np.linalg.inv(box)
    per_molecule = count // molecules

    def image(delta: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        fractional = delta @ inverse
        result: npt.NDArray[np.float64] = (fractional - np.rint(fractional)) @ box
        return result

    angle = 0.0
    urey = 0.0
    for molecule in range(molecules):
        offset = molecule * per_molecule
        for kind, first, centre, last in entries:
            theta0, ktheta, r13, kub = parameters[int(kind)]
            a, b, c = (xyz[offset + int(i)] for i in (first, centre, last))
            u = image(a - b)
            v = image(c - b)
            cosine = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
            theta = math.degrees(math.acos(max(-1.0, min(1.0, cosine))))
            angle += 0.5 * ktheta * math.radians(theta - theta0) ** 2
            urey += 0.5 * kub * (float(np.linalg.norm(image(c - a))) - r13) ** 2
    return {"angle": angle, "urey_bradley": urey}


def _lammps_run(directory: Path, executable: str, input_name: str) -> dict[str, float]:
    completed = subprocess.run(
        [executable, "-in", input_name, "-log", f"{input_name}.log"],
        cwd=directory,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        msg = f"LAMMPS failed:\n{completed.stdout}\n{completed.stderr}"
        raise RuntimeError(msg)
    log = (directory / f"{input_name}.log").read_text(encoding="utf-8")
    header = re.search(r"^\s*Step\s+Temp\s+PotEng.*$", log, flags=re.MULTILINE)
    if header is None:
        msg = "LAMMPS thermo header not found"
        raise RuntimeError(msg)
    names = header.group(0).split()
    row = log[header.end() :].strip().splitlines()[0].split()
    return {name: float(value) * _KCAL for name, value in zip(names, row, strict=True)}


def _lammps(directory: Path, executable: str) -> dict[str, float]:
    full = _lammps_run(directory, executable, "in.chitin")
    data = (directory / "system.data").read_text(encoding="utf-8")

    def replace_section(text: str, title: str, keep: int, zeros: list[str]) -> str:
        head, body = text.split(f"\n{title}\n\n", 1)
        rows, tail = body.split("\n\n", 1)
        edited = "\n".join(" ".join([*row.split()[:keep], *zeros]) for row in rows.splitlines())
        return f"{head}\n{title}\n\n{edited}\n\n{tail}"

    # Urey-Bradley K_ub, r_ub -> 0; CHARMM dihedral 1-4 weight -> 0.
    variant = replace_section(data, "Angle Coeffs", 3, ["0", "0"])
    variant = replace_section(variant, "Dihedral Coeffs", 4, ["0.0"])
    (directory / "variant.data").write_text(variant, encoding="utf-8")
    (directory / "in.variant").write_text(
        (directory / "in.chitin")
        .read_text(encoding="utf-8")
        .replace("read_data system.data", "read_data variant.data"),
        encoding="utf-8",
    )
    reduced = _lammps_run(directory, executable, "in.variant")
    dynamics = _lammps_dynamics(directory, executable)
    offset = _charmmfsw_one_four_offset(data)
    lj_14 = full["E_vdwl"] - reduced["E_vdwl"] - offset
    coulomb_14 = (full["E_coul"] - reduced["E_coul"]) * _COULOMB_RATIO
    switch_offset = _charmmfsw_short_range_offset(data)
    nonbonded = (
        full["E_vdwl"] - offset - switch_offset + (full["E_coul"] + full["E_long"]) * _COULOMB_RATIO
    )
    return {
        "bond": full["E_bond"],
        "angle": reduced["E_angle"],
        "urey_bradley": full["E_angle"] - reduced["E_angle"],
        "proper_dihedral": full["E_dihed"],
        "improper": full.get("E_impro", 0.0),
        "lj_14": lj_14,
        "coulomb_14": coulomb_14,
        "nonbonded": nonbonded,
        "potential": full["PotEng"] - full["E_vdwl"] - full["E_coul"] - full["E_long"] + nonbonded,
        **dynamics,
        "raw_lj_14": full["E_vdwl"] - reduced["E_vdwl"],
        "raw_coulomb_14": full["E_coul"] - reduced["E_coul"],
        "raw_potential": full["PotEng"],
        "charmmfsw_one_four_lj_offset": offset,
        "charmmfsw_minus_gromacs_force_switch_lj": switch_offset,
    }


def _lammps_sections(data: str, title: str) -> list[str]:
    return data.split(f"\n{title}\n\n", 1)[1].split("\n\n", 1)[0].splitlines()


def _charmmfsw_short_range_offset(data: str) -> float:
    """Return CHARMM vfswitch minus GROMACS force-switch short-range LJ energy, kJ/mol.

    Both switch the force identically between 10 and 12 angstrom but shift the
    potential by different constants. Each convention is evaluated over every
    non-excluded pair (exclusions to three bonds) in the LAMMPS data file.
    """
    header = data.split("\nMasses\n", 1)[0].splitlines()

    def bounds(tag: str) -> float:
        low, high = (float(v) for v in next(r for r in header if r.endswith(tag)).split()[:2])
        return high - low

    tilt = next((r for r in header if r.endswith("xy xz yz")), "0 0 0 xy xz yz")
    xy, xz, yz = (float(v) for v in tilt.split()[:3])
    box = np.asarray(
        [[bounds("xlo xhi"), 0, 0], [xy, bounds("ylo yhi"), 0], [xz, yz, bounds("zlo zhi")]]
    )
    rows = [row.split() for row in _lammps_sections(data, "Atoms # full")]
    count = len(rows)
    types = np.zeros(count, dtype=int)
    xyz = np.zeros((count, 3))
    for fields in rows:
        index = int(fields[0]) - 1
        types[index] = int(fields[2])
        xyz[index] = [float(v) for v in fields[4:7]]
    size = int(types.max()) + 1
    epsilon = np.zeros((size, size))
    sigma = np.zeros((size, size))
    for row in _lammps_sections(data, "PairIJ Coeffs"):
        first, second, e, s, *_ = row.split()
        for a, b in ((int(first), int(second)), (int(second), int(first))):
            epsilon[a, b] = float(e)
            sigma[a, b] = float(s)
    neighbours: list[set[int]] = [set() for _ in range(count)]
    for row in _lammps_sections(data, "Bonds"):
        _, _, bonded_first, bonded_second = (int(v) for v in row.split())
        neighbours[bonded_first - 1].add(bonded_second - 1)
        neighbours[bonded_second - 1].add(bonded_first - 1)
    inner, outer = _LAMMPS_SWITCH_INNER, _LAMMPS_SWITCH_OUTER

    def gromacs(r: npt.NDArray[np.float64], power: int) -> npt.NDArray[np.float64]:
        a = (
            -power
            * ((power + 4) * outer - (power + 1) * inner)
            / (outer ** (power + 2) * (outer - inner) ** 2)
        )
        b = (
            power
            * ((power + 3) * outer - (power + 1) * inner)
            / (outer ** (power + 2) * (outer - inner) ** 3)
        )
        c = outer**-power - a / 3 * (outer - inner) ** 3 - b / 4 * (outer - inner) ** 4
        tail = np.clip(r - inner, 0.0, None)
        return np.asarray(r**-power - a / 3 * tail**3 - b / 4 * tail**4 - c)

    def charmm(r: npt.NDArray[np.float64], power: int) -> npt.NDArray[np.float64]:
        half = power // 2
        k = outer**half / (outer**half - inner**half)
        return np.asarray(
            np.where(
                r <= inner,
                r**-power - 1 / (inner**half * outer**half),
                k * (r**-half - outer**-half) ** 2,
            )
        )

    inverse = np.linalg.inv(box)
    total = 0.0
    for index in range(count - 1):
        excluded = {index}
        frontier = {index}
        for _ in range(3):
            frontier = {k for j in frontier for k in neighbours[j]} - excluded
            excluded |= frontier
        delta = xyz[index + 1 :] - xyz[index]
        fractional = delta @ inverse
        distance = np.linalg.norm((fractional - np.rint(fractional)) @ box, axis=1)
        partners = np.arange(index + 1, count)
        mask = (distance < outer) & ~np.isin(partners, list(excluded))
        r = distance[mask]
        e = epsilon[types[index], types[partners[mask]]]
        s = sigma[types[index], types[partners[mask]]]
        repulsion = 4 * e * s**12
        dispersion = 4 * e * s**6
        total += float(
            np.sum(
                repulsion * (charmm(r, 12) - gromacs(r, 12))
                - dispersion * (charmm(r, 6) - gromacs(r, 6))
            )
        )
    return total * _KCAL


def _lammps_dynamics(directory: Path, executable: str) -> dict[str, float]:
    """Run 1 ps NVT at 300 K (0.5 fs, no constraints) from the built coordinates.

    This exercises improper forces, which the single-point comparison does not.
    LAMMPS ``minimize`` is not used: the packaged LAMMPS 22 Jul 2025 build used
    here exits with a segmentation fault after any minimization completes,
    independently of this system's force field.
    """
    text = (directory / "in.chitin").read_text(encoding="utf-8")
    controls = (
        "timestep 0.5\n"
        "velocity all create 300.0 4928459 mom yes rot no dist gaussian\n"
        "fix thermostat all nvt temp 300.0 300.0 50.0\n"
        f"thermo {_DYNAMICS_STEPS}\n"
        f"run {_DYNAMICS_STEPS}\n"
    )
    (directory / "in.dynamics").write_text(text.replace("run 0\n", controls), encoding="utf-8")
    completed = subprocess.run(
        [executable, "-in", "in.dynamics", "-log", "dynamics.log"],
        cwd=directory,
        capture_output=True,
        text=True,
        check=False,
    )
    log = (directory / "dynamics.log").read_text(encoding="utf-8")
    if completed.returncode != 0 or re.search(r"ERROR|\bnan\b|\binf\b", log, re.IGNORECASE):
        msg = f"LAMMPS dynamics smoke test failed:\n{log[-2000:]}"
        raise RuntimeError(msg)
    headers = list(re.finditer(r"^\s*Step\s+Temp\s+PotEng.*$", log, flags=re.MULTILINE))
    names = headers[-1].group(0).split()
    rows = [
        row.split()
        for row in log[headers[-1].end() :].splitlines()
        if row.strip() and row.split()[0].isdigit()
    ]
    final = dict(zip(names, (float(v) for v in rows[-1]), strict=True))
    return {
        "dynamics_final_step": final["Step"],
        "dynamics_final_temperature_k": final["Temp"],
        "dynamics_final_potential": final["PotEng"] * _KCAL,
        "dynamics_final_improper": final.get("E_impro", math.nan) * _KCAL,
    }


def _charmmfsw_one_four_offset(data: str) -> float:
    """Return LAMMPS charmmfsw's constant (force-free) 1-4 LJ energy offset, kJ/mol.

    For each weighted 1-4 pair, dihedral_style charmmfsw adds
    -4 eps sigma^12 / (r_in^6 r_c^6) + 4 eps sigma^6 / (r_in^3 r_c^3).
    """

    def section(title: str) -> list[str]:
        return data.split(f"\n{title}\n\n", 1)[1].split("\n\n", 1)[0].splitlines()

    atom_type = {int(row.split()[0]): int(row.split()[2]) for row in section("Atoms # full")}
    one_four = {
        (int(fields[0]), int(fields[1])): (float(fields[4]), float(fields[5]))
        for fields in (row.split() for row in section("PairIJ Coeffs"))
    }
    weight = {int(row.split()[0]): float(row.split()[4]) for row in section("Dihedral Coeffs")}
    inner, outer = _LAMMPS_SWITCH_INNER, _LAMMPS_SWITCH_OUTER
    total = 0.0
    for row in section("Dihedrals"):
        fields = [int(value) for value in row.split()]
        if weight[fields[1]] != 1.0:
            continue
        first, second = sorted((atom_type[fields[2]], atom_type[fields[5]]))
        epsilon, sigma = one_four[(first, second)]
        total += -4 * epsilon * sigma**12 / (inner**6 * outer**6) + 4 * epsilon * sigma**6 / (
            inner**3 * outer**3
        )
    return total * _KCAL


def _agree(value: float, reference: float, relative: float) -> bool:
    return abs(value - reference) <= max(relative * abs(reference), _BONDED_ABSOLUTE_TOLERANCE)


def _compare(
    engine: dict[str, float], reference: dict[str, float], scale: float
) -> dict[str, dict[str, float | bool]]:
    """Compare terms; nonbonded and potential use a tolerance on ``scale``.

    ``scale`` is the summed magnitude of the GROMACS nonbonded components, so
    PME/PPPM accuracy is judged against the size of the energies that cancel.
    """
    rows: dict[str, dict[str, float | bool]] = {}
    for term, value in engine.items():
        if term not in reference:
            continue
        difference = value - reference[term]
        if term in {"nonbonded", "potential"}:
            match = abs(difference) <= _NONBONDED_RELATIVE_TOLERANCE * scale
        else:
            match = _agree(value, reference[term], _BONDED_RELATIVE_TOLERANCE)
        rows[term] = {
            "engine_kj_mol": round(value, 4),
            "gromacs_kj_mol": round(reference[term], 4),
            "difference_kj_mol": round(difference, 4),
            "match": match,
        }
    return rows


def _nonbonded_scale(gromacs: dict[str, float]) -> float:
    return sum(
        abs(gromacs[term])
        for term in ("LJ-14", "Coulomb-14", "LJ-(SR)", "Coulomb-(SR)", "Coul.-recip.")
    )


def _reference(gromacs: dict[str, float], split: dict[str, float]) -> dict[str, float]:
    return {
        "bond": gromacs["Bond"],
        "angle": split["angle"],
        "urey_bradley": split["urey_bradley"],
        "proper_dihedral": gromacs["Proper-Dih."],
        "improper": gromacs["Improper-Dih."],
        "lj_14": gromacs["LJ-14"],
        "coulomb_14": gromacs["Coulomb-14"],
        "nonbonded": gromacs["LJ-14"]
        + gromacs["Coulomb-14"]
        + gromacs["LJ-(SR)"]
        + gromacs["Coulomb-(SR)"]
        + gromacs["Coul.-recip."],
        "potential": gromacs["Potential"],
    }


def _case(  # noqa: PLR0913, PLR0917
    allomorph: str, boundary: str, work: Path, ff: Path, openmm_python: str, lammps: str
) -> dict[str, object]:
    periodic = boundary == "periodic"
    crystal, geometry = _build(allomorph, "slab", boundary, 2.4)
    dp = len(
        {
            residue
            for chain, residue in zip(crystal.chain_numbers, crystal.residue_numbers, strict=True)
            if chain == crystal.chain_numbers[0]
        }
    )
    topology: PeriodicCharmm36Topology | FiniteCharmm36Topology
    if periodic:
        periodic_topology = native_periodic_charmm36_topology(
            ff, degree_polymerization=dp, residue="BGLCNA"
        )
        topology = periodic_topology
        system = map_periodic_charmm36(crystal, periodic_topology)
    else:
        finite_topology = native_finite_charmm36_topology(
            ff, degree_polymerization=dp, residue="BGLCNA"
        )
        topology = finite_topology
        system = map_finite_charmm36(crystal, finite_topology)
    directory = work / f"{allomorph}_chitin_slab_{boundary}"
    if directory.exists():
        shutil.rmtree(directory)
    _, _, top_path, _ = write_charmm36_system(system, topology, ff, directory, polymer="chitin")

    plain = _gromacs(directory, _PLAIN_MDP, "plain", periodic=periodic)
    switched = _gromacs(directory, _SWITCHED_MDP, "switched", periodic=periodic)
    split = _gromacs_angle_split(directory, "plain.tpr")
    split_consistent = _agree(
        split["angle"] + split["urey_bradley"], plain["U-B"], _BONDED_RELATIVE_TOLERANCE
    )

    request = {
        "coordinates": str(directory / "structure.gro"),
        "topology": str(top_path),
        "include_directory": str(directory),
        "periodic": periodic,
    }
    completed = subprocess.run(
        [openmm_python, "-c", _OPENMM_EVALUATOR, json.dumps(request)],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        msg = f"OpenMM evaluation failed:\n{completed.stderr}"
        raise RuntimeError(msg)
    openmm_energies = json.loads(completed.stdout.strip().splitlines()[-1])
    improper_count = int(openmm_energies.pop("improper_count"))

    lammps_directory = directory / "lammps"
    _, lammps_metadata = write_lammps_system(
        system,
        top_path,
        directory / "structure.gro",
        lammps_directory,
        force_field="charmm36",
        python_executable=openmm_python,
        polymer="chitin",
    )
    lammps_energies = _lammps(lammps_directory, lammps)

    openmm_comparison = _compare(openmm_energies, _reference(plain, split), _nonbonded_scale(plain))
    lammps_comparison = _compare(
        lammps_energies, _reference(switched, split), _nonbonded_scale(switched)
    )
    expected_impropers = 2 * dp * system.chain_count
    checks = {
        "gromacs_angle_split_reproduces_u_b": split_consistent,
        "openmm_improper_count": improper_count == expected_impropers,
        "lammps_improper_count": lammps_metadata.get("improper_count") == expected_impropers,
        "openmm_all_expected_terms_present": set(openmm_comparison) == _EXPECTED_TERMS,
        "lammps_all_expected_terms_present": set(lammps_comparison) == _EXPECTED_TERMS,
        "lammps_dynamics_smoke_test_completed": lammps_energies.get("dynamics_final_step")
        == _DYNAMICS_STEPS
        and all(
            math.isfinite(lammps_energies[key])
            for key in (
                "dynamics_final_temperature_k",
                "dynamics_final_potential",
                "dynamics_final_improper",
            )
        ),
        "openmm_terms_match": all(bool(row["match"]) for row in openmm_comparison.values()),
        "lammps_terms_match": all(bool(row["match"]) for row in lammps_comparison.values()),
    }
    return {
        "case": directory.name,
        **geometry,
        "atoms": len(system.atom_names),
        "impropers": expected_impropers,
        "gromacs_plain_cutoff": plain,
        "gromacs_force_switch": switched,
        "gromacs_angle_split": split,
        "openmm_vs_gromacs_plain_cutoff": openmm_comparison,
        "lammps_vs_gromacs_force_switch": lammps_comparison,
        "lammps_raw": {
            key: round(value, 4)
            for key, value in lammps_energies.items()
            if key.startswith(("raw_", "charmmfsw", "dynamics_"))
        },
        "lammps_conversion": lammps_metadata,
        "checks": checks,
        "passed": all(checks.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--openmm-python", required=True)
    parser.add_argument("--lammps", default=shutil.which("lmp_serial") or "lmp")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    ff = discover_charmm36_directory()
    args.work_directory.mkdir(parents=True, exist_ok=True)
    results = [
        _case(allomorph, boundary, args.work_directory, ff, args.openmm_python, args.lammps)
        for allomorph in ("alpha", "beta")
        for boundary in ("periodic", "finite")
    ]
    report = {
        "test": "chitin CHARMM36 term-by-term energies: GROMACS vs OpenMM and LAMMPS",
        "tolerances": {
            "bonded_relative": _BONDED_RELATIVE_TOLERANCE,
            "bonded_absolute_kj_mol": _BONDED_ABSOLUTE_TOLERANCE,
            "nonbonded_and_potential_relative": _NONBONDED_RELATIVE_TOLERANCE,
        },
        "lammps_improper_conversion": "K_LAMMPS = k_GROMACS / 2, then kJ -> kcal",
        "lammps_conventions_removed_before_comparison": {
            "coulomb_constant": f"LAMMPS CHARMM {_LAMMPS_CHARMM_COULOMB} vs GROMACS "
            f"{_GROMACS_COULOMB:.4f} kcal A/mol/e^2; LAMMPS electrostatics scaled by "
            f"{_COULOMB_RATIO:.8f}",
            "charmmfsw_one_four_lj_offset": "constant per-pair energy offset of dihedral_style "
            "charmmfsw (10-12 angstrom), computed analytically and subtracted",
            "short_range_lj_switch": "CHARMM vfswitch (LAMMPS lj/charmmfsw) and GROMACS "
            "force-switch differ by a constant potential shift; the difference is evaluated "
            "over every non-excluded pair and subtracted",
        },
        "angle_split": "GROMACS reports angle + Urey-Bradley as U-B; the split is evaluated "
        "from GROMACS's resolved parameters (gmx dump) and must sum to U-B",
        "nonbonded_tolerance": "relative to the summed magnitude of GROMACS nonbonded components",
        "results": results,
        "passed": all(result["passed"] for result in results),
    }
    text = json.dumps(report, indent=2)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
