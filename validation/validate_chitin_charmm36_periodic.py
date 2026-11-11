"""Validate periodic alpha- and beta-chitin CHARMM36 bulk systems with GROMACS.

For each allomorph this builds a covalently periodic bulk crystal, maps it to
the native BGLCNA/14bb topology, and checks:

- every boundary-crossing and internal O4--C1 bond of every chain, using
  minimum-image distances (expected 1.4--1.5 angstrom);
- 27 atoms and zero charge per periodic residue;
- two function-2 impropers per residue;
- ``gmx grompp`` without warnings, and finite zero-step energies;
- steepest-descent minimization to Fmax <= 1000 kJ/mol/nm and a 1 ps NVT smoke
  test, with all O4--C1 bonds rechecked in the final NVT coordinates.

It calls the library directly because the CLI still refuses chitin CHARMM36
output until this validation passes.

Usage:
    python validation/validate_chitin_charmm36_periodic.py --work-directory DIR \
        [--minimum-box-nm 2.4] [--output chitin_charmm36_periodic.json]
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np

from cellulose_forge.charmm36 import (
    map_periodic_charmm36,
    native_periodic_charmm36_topology,
    write_charmm36_system,
)
from cellulose_forge.crystallography import load_allomorph
from cellulose_forge.forcefields import discover_charmm36_directory
from cellulose_forge.geometry import build_bulk_crystal

_ROOT = Path(__file__).resolve().parents[1]
_ATOMS_PER_RESIDUE = 27
_IMPROPERS_PER_RESIDUE = 2
_CHARGE_TOLERANCE = 1e-6
_BOND_RANGE = (1.35, 1.55)
# After 1 ps at 300 K, thermal fluctuation widens the accepted O4--C1 range.
_DYNAMICS_BOND_RANGE = (1.25, 1.70)
_MAXIMUM_INITIAL_BOND_PER_ATOM = 50.0
_MAXIMUM_MINIMIZED_FORCE = 1000.0
_NVT_FINAL_TIME_PS = 1.0
# mdrun echoes input parameters as 'name = value' (e.g. epsilon-rf = inf).
_PARAMETER_ECHO = re.compile(r"^\s+[\w-]+(?:\[[^]]*\])?\s+=\s")
_LOG_PROBLEMS = re.compile(r"WARNING|LINCS WARNING|\bnan\b|\binf\b|not finite", re.IGNORECASE)
_TERMS = (
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
_SINGLE_POINT_MDP = """\
integrator          = md
nsteps              = 0
dt                  = 0.001
continuation        = yes
constraints         = none
cutoff-scheme       = Verlet
nstlist             = 10
rlist               = 1.0
coulombtype         = PME
rcoulomb            = 1.0
vdwtype             = Cut-off
rvdw                = 1.0
pbc                 = xyz
periodic-molecules  = yes
nstenergy           = 1
nstcalcenergy       = 1
"""


def _gmx(*arguments: str, cwd: Path, stdin: str | None = None, check: bool = True) -> str:
    completed = subprocess.run(
        [shutil.which("gmx") or "gmx", *arguments],
        cwd=cwd,
        input=stdin,
        capture_output=True,
        text=True,
        check=False,
    )
    if check and completed.returncode != 0:
        msg = f"gmx {' '.join(arguments)} failed:\n{completed.stdout}\n{completed.stderr}"
        raise RuntimeError(msg)
    return completed.stdout + completed.stderr


def _grompp(directory: Path, mdp: Path, tpr: str, coordinates: str = "structure.gro") -> str:
    output = _gmx(
        "grompp",
        "-f",
        str(mdp),
        "-c",
        coordinates,
        "-p",
        "system.top",
        "-o",
        tpr,
        "-po",
        f"{tpr}.mdout.mdp",
        "-maxwarn",
        "0",
        cwd=directory,
    )
    return output


def _energies(directory: Path, edr: str) -> dict[str, float]:
    """Return the selected energy terms from the last frame of an EDR file."""
    listing = _gmx("energy", "-f", edr, "-o", f"{edr}.xvg", cwd=directory, stdin="0\n", check=False)
    numbers = {
        name.replace(" ", "-"): number
        for number, name in re.findall(
            r"(\d+)\s+(\S+)", listing.split("End your selection")[1].split("Program:")[0]
        )
    }
    wanted = [numbers[term] for term in _TERMS if term in numbers]
    _gmx("energy", "-f", edr, "-o", f"{edr}.xvg", cwd=directory, stdin="\n".join(wanted) + "\n0\n")
    xvg = (directory / f"{edr}.xvg").read_text(encoding="utf-8").splitlines()
    rows = [line.split() for line in xvg if line and line[0] not in "#@"]
    legends = [
        match.group(1).replace(" ", "-")
        for line in xvg
        if (match := re.match(r'@ s\d+ legend "(.*)"', line))
    ]
    return {name: float(value) for name, value in zip(legends, rows[-1][1:], strict=True)}


def _temperature(directory: Path, edr: str) -> float:
    output = _gmx(
        "energy", "-f", edr, "-o", "temperature.xvg", cwd=directory, stdin="Temperature\n0\n"
    )
    match = re.search(r"^Temperature\s+([-+0-9.eE]+)", output, flags=re.MULTILINE)
    if match is None:
        msg = "could not read the NVT temperature"
        raise RuntimeError(msg)
    return float(match.group(1))


def _read_gro(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return GRO coordinates and box vectors in angstrom."""
    lines = path.read_text(encoding="utf-8").splitlines()
    count = int(lines[1])
    coordinates = np.asarray(
        [
            [float(line[20:28]), float(line[28:36]), float(line[36:44])]
            for line in lines[2 : 2 + count]
        ]
    )
    values = [float(value) for value in lines[2 + count].split()] + [0.0] * 6
    box = np.asarray(
        [
            [values[0], values[3], values[4]],
            [values[5], values[1], values[6]],
            [values[7], values[8], values[2]],
        ]
    )
    return 10.0 * coordinates, 10.0 * box


def _periodic_bonds(
    system,
    coordinates: np.ndarray | None = None,
    basis: np.ndarray | None = None,
    bond_range: tuple[float, float] = _BOND_RANGE,
) -> dict[str, float | int | bool]:
    """Check every O4--C1 bond (internal and boundary) with minimum images."""
    basis = np.asarray(system.basis_vectors if basis is None else basis, dtype=float)
    inverse = np.linalg.inv(basis)
    coordinates = np.asarray(
        system.cartesian_coordinates if coordinates is None else coordinates, dtype=float
    )
    lookup = {
        (chain, residue, name): index
        for index, (name, residue, chain) in enumerate(
            zip(system.atom_names, system.residue_numbers, system.chain_numbers, strict=True)
        )
    }
    distances: list[float] = []
    crossing: list[float] = []
    chains = tuple(dict.fromkeys(system.chain_numbers))
    for chain in chains:
        residues = tuple(
            dict.fromkeys(
                residue
                for residue, owner in zip(system.residue_numbers, system.chain_numbers, strict=True)
                if owner == chain
            )
        )
        for residue, following in zip(residues, residues[1:] + residues[:1], strict=True):
            raw = (
                coordinates[lookup[(chain, following, "C1")]]
                - coordinates[lookup[(chain, residue, "O4")]]
            )
            fractional = raw @ inverse
            image = fractional - np.rint(fractional)
            distance = float(np.linalg.norm(image @ basis))
            distances.append(distance)
            if following == residues[0]:
                crossing.append(distance)
    return {
        "within_range": len(crossing) == len(chains)
        and all(bond_range[0] <= d <= bond_range[1] for d in distances),
        "range_angstrom": list(bond_range),
        "chains": len(chains),
        "bonds_checked": len(distances),
        "boundary_crossing_bonds": len(crossing),
        "minimum_angstrom": round(min(distances), 4),
        "maximum_angstrom": round(max(distances), 4),
        "boundary_minimum_angstrom": round(min(crossing), 4),
        "boundary_maximum_angstrom": round(max(crossing), 4),
    }


def _validate(allomorph: str, minimum_box_nm: float, work: Path, ff: Path) -> dict:
    definition = load_allomorph(allomorph, polymer="chitin")
    cell = definition.unit_cell
    replication = tuple(
        math.ceil(10.0 * minimum_box_nm / length) for length in (cell.a, cell.b, cell.c)
    )
    crystal = build_bulk_crystal(definition, replication)
    residues_per_chain = len(
        {
            residue
            for chain, residue in zip(crystal.chain_numbers, crystal.residue_numbers, strict=True)
            if chain == crystal.chain_numbers[0]
        }
    )
    topology = native_periodic_charmm36_topology(
        ff, degree_polymerization=residues_per_chain, residue="BGLCNA"
    )
    system = map_periodic_charmm36(crystal, topology)
    bonds = _periodic_bonds(system)

    residue_atoms: dict[int, list] = defaultdict(list)
    for atom in topology.topology_atoms:
        residue_atoms[atom.residue_number].append(atom)
    if not all(
        len(atoms) == _ATOMS_PER_RESIDUE and abs(sum(a.charge for a in atoms)) < _CHARGE_TOLERANCE
        for atoms in residue_atoms.values()
    ):
        msg = "periodic residues are not 27-atom neutral repeats"
        raise ValueError(msg)

    directory = work / f"{allomorph}_chitin"
    if directory.exists():
        shutil.rmtree(directory)
    write_charmm36_system(system, topology, ff, directory)
    itp = (directory / "cellulose_charmm36.itp").read_text(encoding="utf-8")
    impropers = sum(
        1
        for line in itp.splitlines()
        if len(line.split()) == 5 and line.split()[-1] == "2"  # noqa: PLR2004
    )
    if impropers != _IMPROPERS_PER_RESIDUE * residues_per_chain:
        msg = f"expected {2 * residues_per_chain} impropers per chain, found {impropers}"
        raise ValueError(msg)

    (directory / "single_point.mdp").write_text(_SINGLE_POINT_MDP, encoding="utf-8")
    logs = [_grompp(directory, directory / "single_point.mdp", "single_point.tpr")]
    logs.append(_mdrun(directory, "single_point"))
    initial = _energies(directory, "single_point.edr")

    logs.append(_grompp(directory, _ROOT / "validation/periodic_pme_minim.mdp", "minim.tpr"))
    minim_log = _mdrun(directory, "minim")
    logs.append(minim_log)
    minimized = _energies(directory, "minim.edr")
    fmax_match = re.search(r"Maximum force\s*=\s*([-+0-9.eE]+)", minim_log)
    fmax = float(fmax_match.group(1)) if fmax_match else math.nan

    logs.append(
        _grompp(directory, _ROOT / "validation/periodic_pme_nvt_1ps.mdp", "nvt.tpr", "minim.gro")
    )
    logs.append(_mdrun(directory, "nvt"))
    final = _energies(directory, "nvt.edr")
    temperature = _temperature(directory, "nvt.edr")
    last_frame = re.findall(
        r"Last frame\s+(\d+)\s+time\s+([-+0-9.eE]+)",
        _gmx("check", "-f", "nvt.xtc", cwd=directory),
    )
    final_time = float(last_frame[-1][1]) if last_frame else math.nan
    final_coordinates, final_box = _read_gro(directory / "nvt.gro")
    final_bonds = _periodic_bonds(
        system, final_coordinates, final_box, bond_range=_DYNAMICS_BOND_RANGE
    )

    atom_count = len(system.atom_names)
    bond_per_atom = initial.get("Bond", math.nan) / atom_count
    problems = sorted(
        {
            match.group(0)
            for log in logs
            for line in log.splitlines()
            if not _PARAMETER_ECHO.match(line)
            for match in _LOG_PROBLEMS.finditer(line)
        }
    )
    checks = {
        "initial_bonds_within_range": bool(bonds["within_range"]),
        "all_energy_terms_present_and_finite": all(
            set(_TERMS) <= set(energies) and all(map(math.isfinite, energies.values()))
            for energies in (initial, minimized, final)
        ),
        "initial_bond_term_per_atom_below_50": bond_per_atom < _MAXIMUM_INITIAL_BOND_PER_ATOM,
        "minimization_maximum_force_at_most_1000": fmax <= _MAXIMUM_MINIMIZED_FORCE,
        "no_warning_lincs_or_nonfinite_messages": not problems,
        "nvt_final_frame_present": math.isclose(final_time, _NVT_FINAL_TIME_PS, abs_tol=1e-6),
        "final_nvt_bonds_within_range": bool(final_bonds["within_range"]),
    }
    return {
        "allomorph": f"chitin/{allomorph}",
        "replication": list(replication),
        "atoms": atom_count,
        "residues_per_chain": residues_per_chain,
        "periodic_bonds_initial": bonds,
        "periodic_bonds_final_nvt": final_bonds,
        "impropers_per_chain": impropers,
        "log_problems": problems,
        "initial_energy_kj_mol": initial,
        "initial_bond_energy_per_atom_kj_mol": round(bond_per_atom, 4),
        "minimized_energy_kj_mol": minimized,
        "minimization_maximum_force_kj_mol_nm": fmax,
        "nvt_final_frame_time_ps": final_time,
        "nvt_final_energy_kj_mol": final,
        "nvt_smoke_test_mean_temperature_k": temperature,
        "checks": checks,
        "passed": all(checks.values()),
    }


def _mdrun(directory: Path, name: str) -> str:
    """Run one single-rank mdrun and return its terminal output plus log file."""
    output = _gmx(
        "mdrun", "-s", f"{name}.tpr", "-deffnm", name, "-ntmpi", "1", "-ntomp", "1", cwd=directory
    )
    return output + (directory / f"{name}.log").read_text(encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument(
        "--minimum-box-nm", type=float, default=2.4, help="minimum supercell edge (nm)"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    ff = discover_charmm36_directory()
    args.work_directory.mkdir(parents=True, exist_ok=True)
    results = [
        _validate(allomorph, args.minimum_box_nm, args.work_directory, ff)
        for allomorph in ("alpha", "beta")
    ]
    report = {
        "test": "periodic chitin CHARMM36 bulk (BGLCNA + 14bb) with GROMACS",
        "note": (
            "the 1 ps NVT is a smoke test, not equilibration or temperature validation; "
            "the Bond term is GROMACS's bond-stretch term only"
        ),
        "gromacs_force_field": ff.name,
        "gromacs_version": _gmx("--version", cwd=args.work_directory)
        .split("GROMACS version:")[1]
        .split()[0],
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
