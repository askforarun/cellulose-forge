"""Validate the experimental OPLS-AA and OPLS-AA-SEI cellulose outputs.

Checks, for GROMACS (mixed and double precision), OpenMM (Reference platform) and LAMMPS:

1. Transcription: Kony et al. (2002) Table 3 "OPLS-AA" column, converted to
   Ryckaert-Bellemans form, and the explicit bond and angle parameters must
   equal the corresponding types in the installed GROMACS ``oplsaa.ff`` (an
   independent copy of the Damm et al. 1997 parameters).
2. Independent recomputation: torsion, 1-4 and SEI-pair energies recomputed
   here in NumPy from the written GROMACS topology must match GROMACS, and
   the OPLS-AA-SEI minus OPLS-AA nonbonded energy must equal the predicted
   effect of scaling each SEI pair's Coulomb term.
3. Cross-engine: single-point energies of identical, unrounded coordinates
   from double-precision GROMACS, OpenMM and LAMMPS must
   agree term by term and in total.  The comparison uses unshifted cutoffs,
   tight Ewald settings and no dispersion correction in every engine.
4. Runtime: every geometry offered (finite and periodic slabs, periodic bulk,
   36-chain fibrils) for every allomorph passes ``gmx grompp -maxwarn 0``,
   minimization and a 1 ps NVT smoke test; the serialized OpenMM system
   minimizes and runs 1 ps; the generated LAMMPS input runs.

This validates topology construction and engine translation only, not the
accuracy of either force field for crystalline cellulose.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "validation"))

from _provenance import provenance  # noqa: E402

from cellulose_forge.crystallography import MillerIndex, load_allomorph  # noqa: E402
from cellulose_forge.geometry import (  # noqa: E402
    build_36_chain_fibril,
    build_bulk_crystal,
    build_chain_parallel_slab,
)
from cellulose_forge.opls import (  # noqa: E402
    _ANGLE_PARAMETERS,
    _BOND_PARAMETERS,
    _LJ,
    _TORSIONS,
    KCAL_TO_KJ,
    OplsBuild,
    build_opls_system,
    rb_coefficients,
    write_opls_system,
)
from cellulose_forge.opls_engines import write_opls_lammps, write_opls_openmm  # noqa: E402
from cellulose_forge.run_inputs import (  # noqa: E402
    dispersion_correction,
    write_gromacs_run_inputs,
)

COULOMB = 138.935458  # kJ mol-1 nm e-2, GROMACS ONE_4PI_EPS0
# Individual terms (bonded, 1-4/SEI pairs, real-space Lennard-Jones) are
# compared against double-precision GROMACS; measured agreement is 1e-7 or
# better, limited by OpenMM/LAMMPS internal arithmetic.
TERM_RELATIVE_TOLERANCE = 1.0e-6
# Totals include reciprocal-space Coulomb from three different Ewald
# implementations (GROMACS PME, OpenMM PME, LAMMPS PPPM) at tight settings;
# measured agreement is about 1e-6.
TOTAL_RELATIVE_TOLERANCE = 1.0e-5
PREDICTION_RELATIVE_TOLERANCE = 1.0e-5
# Boxes with a tilted third vector (I-alpha bulk) are evaluated on the OpenMM
# CPU platform (single-precision nonbonded kernels), because OpenMM 8.6.1
# Reference nonbonded energies disagreed with exact sums for such boxes in
# measure_openmm_tilted_box.py; see cellulose_forge.opls_engines.validation_platform.
OPENMM_CPU_RELATIVE_TOLERANCE = 1.0e-5
SLAB_FACES = {"i_alpha": (0, 1, 0), "i_beta": (1, 1, 0), "ii": (1, 1, 0), "iii_i": (1, 1, 0)}
BULK_CELLS = {"i_alpha": (3, 4, 4), "i_beta": (3, 3, 3), "ii": (3, 3, 3), "iii_i": (5, 3, 3)}
VARIANTS = ("opls-aa", "opls-aa-sei")
POINT_MDP = """integrator = md
nsteps = 0
cutoff-scheme = Verlet
coulombtype = PME
coulomb-modifier = None
rcoulomb = 1.0
fourierspacing = 0.06
ewald-rtol = 1e-6
pme-order = 8
vdwtype = Cut-off
vdw-modifier = None
rvdw = 1.0
DispCorr = no
constraints = none
pbc = xyz
periodic-molecules = {periodic}
nstenergy = 1
nstcalcenergy = 1
"""
GROMACS_TERMS = ("Bond", "Angle", "Ryckaert-Bell.", "LJ-14", "Coulomb-14", "LJ-(SR)",
                 "Coulomb-(SR)", "Coul.-recip.", "Potential")
# GROMACS oplsaa.ff dihedral types and the Kony class they must reproduce.
TRANSCRIPTION = {
    "O-C-O-C": ("CT", "OS", "CO", "OS"),
    "O-C-O-H": ("HO", "OH", "CO", "OS"),
    "OH-C-C-OH": ("OH", "CT", "CT", "OH"),
    "OS-C-C-OH": ("OH", "CT", "CT", "OS"),
    "OS-C-C-OS": ("OS", "CT", "CT", "OS"),
}


def _gmx_version(gmx: str) -> str:
    output = _run([gmx, "--version"], ROOT)
    return output.split("GROMACS version:")[1].split()[0] + (
        " (double)" if "Precision:           double" in output else " (mixed)"
    )


def _gmx_share() -> Path:
    result = subprocess.run(["gmx", "--version"], capture_output=True, text=True, check=True)
    match = re.search(r"Data prefix:\s+(\S+)", result.stdout + result.stderr)
    if match is None:
        raise RuntimeError("could not locate the GROMACS data directory")
    return Path(match.group(1)) / "share" / "gromacs" / "top" / "oplsaa.ff"


def _ffbonded() -> dict[str, dict[tuple[str, ...], tuple[float, ...]]]:
    tables: dict[str, dict[tuple[str, ...], tuple[float, ...]]] = {}
    section = ""
    width = {"bondtypes": 2, "angletypes": 3, "dihedraltypes": 4}
    for raw in (_gmx_share() / "ffbonded.itp").read_text(encoding="utf-8").splitlines():
        line = raw.split(";")[0].strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("["):
            section = line.strip("[] ")
            continue
        if section in width:
            fields = line.split()
            n = width[section]
            tables.setdefault(section, {})[tuple(fields[:n])] = tuple(
                float(value) for value in fields[n + 1:]
            )
    return tables


def transcription_check() -> dict[str, object]:
    """Compare explicit parameters with GROMACS oplsaa.ff."""
    tables = _ffbonded()

    def lookup(section: str, key: tuple[str, ...]) -> tuple[float, ...]:
        found = tables[section].get(key) or tables[section].get(key[::-1])
        if found is None:
            raise RuntimeError(f"oplsaa.ff has no {section} entry {key}")
        return found

    torsions = {}
    for name, types in TRANSCRIPTION.items():
        gromacs = lookup("dihedraltypes", types)
        ours = rb_coefficients(_TORSIONS["opls-aa"][name])
        difference = max(abs(a - b) for a, b in zip(ours, gromacs, strict=True))
        # Kony prints three decimals in kcal/mol; allow that rounding.
        torsions[name] = {
            "gromacs_types": "-".join(types),
            "max_rb_difference_kj_mol": difference,
            "passed": difference <= 0.002 * KCAL_TO_KJ * 2,
        }
    bonds = {
        "-".join(key): lookup("bondtypes", key)[:2] == value
        for key, value in _BOND_PARAMETERS.items()
    }
    angles = {
        "-".join(key): lookup("angletypes", key)[:2] == value
        for key, value in _ANGLE_PARAMETERS.items()
    }
    return {
        "torsion_classes": torsions,
        "bond_types_equal": bonds,
        "angle_types_equal": angles,
        "passed": all(item["passed"] for item in torsions.values())
        and all(bonds.values())
        and all(angles.values()),
    }


def _crystal(kind: str, allomorph: str, boundary: str) -> object:
    definition = load_allomorph(allomorph)
    padding = 12.0 if boundary == "finite" else None
    if kind == "slab":
        return build_chain_parallel_slab(
            definition, MillerIndex(*SLAB_FACES[allomorph]), layers=3, chains_per_layer=4,
            degree_polymerization=6, vacuum=30.0, padding_z=padding or 0.0,
            chain_boundary=boundary,
        )
    if kind == "bulk":
        return build_bulk_crystal(definition, BULK_CELLS[allomorph])
    return build_36_chain_fibril(
        definition, degree_polymerization=6, vacuum=20.0, padding_z=padding,
        chain_boundary=boundary,
    )


def _write_g96(build: OplsBuild, path: Path) -> None:
    """Write full-precision coordinates and box for the GROMACS energy reference.

    A GRO file rounds positions to 0.001 nm, which places many pairs exactly
    on the 1.0 nm cutoff, where engines legitimately differ; every engine
    therefore receives the same unrounded geometry.
    """
    mapped = build.mapped
    xyz = np.asarray(mapped.cartesian_coordinates, dtype=np.float64) * 0.1
    a, b, c = np.asarray(mapped.basis_vectors, dtype=np.float64) * 0.1
    lines = ["TITLE", "CelluloseForge OPLS validation", "END", "POSITION"]
    for index, (name, residue, position) in enumerate(
        zip(mapped.atom_names, mapped.residue_numbers, xyz, strict=True), start=1
    ):
        lines.append(
            f"{residue % 100000:5d} BGC   {name:<5s}{index % 10000000:7d}"
            f"{position[0]:15.9f}{position[1]:15.9f}{position[2]:15.9f}"
        )
    box = (a[0], b[1], c[2], a[1], a[2], b[0], b[2], c[0], c[1])
    lines += ["END", "BOX", "".join(f"{value:15.9f}" for value in box), "END"]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _read_g96(path: Path) -> tuple[np.ndarray, np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
    start, box_line = lines.index("POSITION") + 1, lines[lines.index("BOX") + 1]
    stop = lines.index("END", start)
    xyz = np.array([[float(value) for value in line[24:].split()] for line in lines[start:stop]])
    fields = [float(value) for value in box_line.split()]
    box = np.zeros((3, 3))
    box[0, 0], box[1, 1], box[2, 2] = fields[:3]
    box[0, 1], box[0, 2], box[1, 0], box[1, 2], box[2, 0], box[2, 1] = fields[3:]
    return xyz, box


def _run(command: list[str], directory: Path, stdin: str | None = None) -> str:
    result = subprocess.run(
        command, cwd=directory, input=stdin, capture_output=True, text=True, check=False
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"{' '.join(command)} failed in {directory}:\n{result.stdout}\n{result.stderr}"
        )
    return result.stdout + result.stderr


def _gromacs_terms(directory: Path, periodic: bool, gmx: str) -> dict[str, float]:
    (directory / "point.mdp").write_text(
        POINT_MDP.format(periodic="yes" if periodic else "no"), encoding="utf-8"
    )
    _run([gmx, "grompp", "-f", "point.mdp", "-c", "structure.g96", "-p", "system.top",
          "-o", "point.tpr", "-maxwarn", "0"], directory)
    _run([gmx, "mdrun", "-deffnm", "point", "-ntomp", "1"], directory)
    _run([gmx, "energy", "-f", "point.edr", "-o", "point.xvg"], directory,
         stdin="\n".join(GROMACS_TERMS) + "\n\n")
    legends, values = [], []
    for line in (directory / "point.xvg").read_text(encoding="utf-8").splitlines():
        match = re.match(r'@ s\d+ legend "(.+)"', line)
        if match:
            # gmx energy selects "LJ-(SR)" but labels it "LJ (SR)".
            legends.append(match.group(1).replace(" ", "-"))
        elif line and line[0] not in "#@":
            values = [float(value) for value in line.split()[1:]]
    return dict(zip(legends, values, strict=True))


def _parse_itp(path: Path) -> dict[str, list[list[str]]]:
    sections: dict[str, list[list[str]]] = {}
    current = ""
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split(";")[0].strip()
        if not line:
            continue
        if line.startswith("["):
            current = line.strip("[] ")
            continue
        sections.setdefault(current, []).append(line.split())
    return sections


def _minimum_image(delta: np.ndarray, box: np.ndarray, periodic: bool) -> np.ndarray:
    if not periodic:
        return delta
    fractional = delta @ np.linalg.inv(box)
    fractional -= np.round(fractional)
    return fractional @ box


NUMPY_GROMACS_TERMS = ("Ryckaert-Bell.", "LJ-14", "Coulomb-14")


def _numpy_terms(directory: Path, chains: int, periodic: bool) -> dict[str, float]:
    """Recompute torsion, 1-4 and SEI-pair energies from the written files.

    ``SEI-LJ`` and ``SEI-Coulomb`` are the parts of ``LJ-14`` and
    ``Coulomb-14`` contributed by the GROMACS function-2 (SEI) pairs.
    """
    itp = _parse_itp(directory / "cellulose_opls.itp")
    atoms = itp["atoms"]
    n = len(atoms)
    types = [row[1] for row in atoms]
    charges = np.array([float(row[6]) for row in atoms])
    xyz, box = _read_g96(directory / "structure.g96")
    torsion = lj14 = coulomb14 = sei_lj = sei_coulomb = 0.0
    for chain in range(chains):
        position = xyz[chain * n: (chain + 1) * n]
        for row in itp["dihedrals"]:
            a, b, c, d = (int(value) - 1 for value in row[:4])
            coefficients = [float(value) for value in row[5:11]]
            b1 = _minimum_image(position[b] - position[a], box, periodic)
            b2 = _minimum_image(position[c] - position[b], box, periodic)
            b3 = _minimum_image(position[d] - position[c], box, periodic)
            n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
            m1 = np.cross(n1, b2 / np.linalg.norm(b2))
            phi = np.arctan2(np.dot(m1, n2), np.dot(n1, n2))
            cosine = np.cos(phi - np.pi)
            torsion += sum(value * cosine**power for power, value in enumerate(coefficients))
        for row in itp["pairs"]:
            a, b = int(row[0]) - 1, int(row[1]) - 1
            r = np.linalg.norm(_minimum_image(position[b] - position[a], box, periodic))
            if row[2] == "1":
                sigma = np.sqrt(_LJ[types[a]][0] * _LJ[types[b]][0])
                epsilon = np.sqrt(_LJ[types[a]][1] * _LJ[types[b]][1])
                fudge_qq, fudge_lj = 0.5, 0.5
            else:
                fudge_qq, sigma, epsilon = float(row[3]), float(row[6]), float(row[7])
                fudge_lj = 1.0
            lj = fudge_lj * 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)
            coulomb = fudge_qq * COULOMB * charges[a] * charges[b] / r
            lj14 += lj
            coulomb14 += coulomb
            if row[2] == "2":
                sei_lj += lj
                sei_coulomb += coulomb
    return {
        "Ryckaert-Bell.": torsion, "LJ-14": lj14, "Coulomb-14": coulomb14,
        "SEI-LJ": sei_lj, "SEI-Coulomb": sei_coulomb,
    }


def _sei_prediction(directory: Path, chains: int, periodic: bool) -> float:
    """Energy change from replacing full Coulomb by the SEI-scaled pair term."""
    itp = _parse_itp(directory / "cellulose_opls.itp")
    charges = np.array([float(row[6]) for row in itp["atoms"]])
    n = len(charges)
    xyz, box = _read_g96(directory / "structure.g96")
    total = 0.0
    for chain in range(chains):
        position = xyz[chain * n: (chain + 1) * n]
        for row in itp["pairs"]:
            if row[2] != "2":
                continue
            a, b, fudge = int(row[0]) - 1, int(row[1]) - 1, float(row[3])
            r = np.linalg.norm(_minimum_image(position[b] - position[a], box, periodic))
            total += (fudge - 1.0) * COULOMB * charges[a] * charges[b] / r
    return total


def _relative(first: float, second: float) -> float:
    return abs(first - second) / max(abs(first), abs(second), 1.0)


def _lammps_terms(directory: Path) -> dict[str, float]:
    """Run the generated LAMMPS input at tight PPPM accuracy (kJ/mol)."""
    source = (directory / "in.cellulose").read_text(encoding="utf-8")
    source = source.replace("kspace_style pppm 1.0e-5", "kspace_style pppm 1.0e-7")
    hybrid = "bond_style hybrid" in source
    if hybrid:
        # compute bond reports each hybrid sub-style: [1] harmonic, [2] SEI table.
        source = source.replace(
            "thermo_style custom step temp pe ebond eangle edihed evdwl ecoul elong",
            "compute bondparts all bond\n"
            "thermo_style custom step temp pe ebond eangle edihed evdwl ecoul elong "
            "c_bondparts[1] c_bondparts[2]",
        )
    (directory / "in.point").write_text(source, encoding="utf-8")
    log = _run(["lmp_serial", "-in", "in.point", "-log", "none"], directory)
    lines = log.splitlines()
    header = next(index for index, line in enumerate(lines) if line.split()[:2] == ["Step", "Temp"])
    names = lines[header].split()
    values = [float(value) for value in lines[header + 1].split()]
    row = {name: value * KCAL_TO_KJ for name, value in zip(names, values, strict=True)}
    return {
        "potential": row["PotEng"],
        "bond": row["c_bondparts[1]"] if hybrid else row["E_bond"],
        "sei_pairs": row["c_bondparts[2]"] if hybrid else 0.0,
        "angle": row["E_angle"],
        "torsion": row["E_dihed"],
        "pair_lennard_jones": row["E_vdwl"],
        "pair_coulomb": row["E_coul"] + row["E_long"],
    }


def _cross_engine(
    build: OplsBuild, directory: Path, python: str, periodic: bool, gmx: str
) -> dict[str, object]:
    """Compare single points from the three engines (no dispersion correction)."""
    write_opls_system(build, directory)
    _write_g96(build, directory / "structure.g96")
    gromacs = _gromacs_terms(directory, periodic, gmx)
    _, openmm = write_opls_openmm(
        build, directory / "openmm", dispersion_correction=False, python_executable=python,
        minimization_iterations=1, validation_ewald_error_tolerance=1.0e-6,
    )
    write_opls_lammps(build, directory, dispersion_correction=False)
    lammps = _lammps_terms(directory)
    energies = openmm["initial_energies_kj_mol"]  # type: ignore[index]
    exceptions = energies["NonbondedForce exceptions (1-4)"]
    sei = _numpy_terms(directory, build.mapped.chain_count, periodic)
    # Each category: (reference, OpenMM value, LAMMPS value, total-type?).  The
    # LAMMPS pair style sees ordinary pairs and 0.5-scaled 1-4 pairs; the SEI
    # pairs are in its tables.  Their references therefore subtract the SEI
    # pair contributions (recomputed independently and checked against
    # GROMACS LJ-14/Coulomb-14 in numpy_recomputation).
    categories = {
        "bond": (gromacs["Bond"], energies["HarmonicBondForce"], lammps["bond"], False),
        "angle": (gromacs["Angle"], energies["HarmonicAngleForce"], lammps["angle"], False),
        "torsion": (
            gromacs["Ryckaert-Bell."], energies["RBTorsionForce"], lammps["torsion"], False
        ),
        "pairs_1_4_and_sei": (gromacs["LJ-14"] + gromacs["Coulomb-14"], exceptions, None, False),
        "sei_pairs": (sei["SEI-LJ"] + sei["SEI-Coulomb"], None, lammps["sei_pairs"], False),
        "lennard_jones": (gromacs["LJ-(SR)"], energies["OPLSLennardJones"], None, False),
        "lammps_pair_lennard_jones": (
            gromacs["LJ-(SR)"] + gromacs["LJ-14"] - sei["SEI-LJ"], None,
            lammps["pair_lennard_jones"], False,
        ),
        "lammps_pair_coulomb": (
            gromacs["Coulomb-(SR)"] + gromacs["Coul.-recip."] + gromacs["Coulomb-14"]
            - sei["SEI-Coulomb"], None, lammps["pair_coulomb"], True,
        ),
        "potential": (
            gromacs["Potential"], openmm["initial_potential_energy_kj_mol"],
            lammps["potential"], True,
        ),
    }
    platform = str(openmm["validation_platform"])
    comparison: dict[str, object] = {"openmm_platform": platform}
    # Ewald errors are absolute and set by the accuracy settings, while a net
    # Coulomb energy can be a small difference of large numbers; Ewald-dependent
    # categories are therefore normalized by the magnitude of the total potential.
    scale = max(abs(gromacs["Potential"]), 1.0)
    for name, (reference, openmm_value, lammps_value, is_total) in categories.items():
        tolerance = TOTAL_RELATIVE_TOLERANCE if is_total else TERM_RELATIVE_TOLERANCE

        def score(value: float, is_total: bool = is_total, reference: float = reference) -> float:
            if is_total:
                return abs(value - reference) / max(scale, abs(reference))
            return _relative(reference, value)

        entry: dict[str, object] = {
            "reference_kj_mol": reference,
            "normalization": "max(|reference|, |potential|)" if is_total else "max(|values|, 1)",
        }
        passed = True
        if openmm_value is not None:
            openmm_tolerance = OPENMM_CPU_RELATIVE_TOLERANCE if platform == "CPU" else tolerance
            difference = score(float(openmm_value))  # type: ignore[arg-type]
            entry.update({"openmm_kj_mol": openmm_value, "openmm_relative_difference": difference,
                          "openmm_tolerance": openmm_tolerance})
            passed = difference <= openmm_tolerance
        if lammps_value is not None:
            difference = score(lammps_value)
            entry.update({"lammps_kj_mol": lammps_value, "lammps_relative_difference": difference,
                          "lammps_absolute_difference_kj_mol": abs(lammps_value - reference),
                          "lammps_tolerance": tolerance})
            passed = passed and difference <= tolerance
        entry["passed"] = passed
        comparison[name] = entry
    return {"gromacs_terms": gromacs, "terms": comparison}


def _gromacs_smoke(directory: Path) -> dict[str, object]:
    logs = _run(["gmx", "grompp", "-f", "minim.mdp", "-c", "structure.gro", "-p", "system.top",
                 "-o", "minim.tpr", "-maxwarn", "0"], directory)
    logs += _run(["gmx", "mdrun", "-deffnm", "minim", "-ntmpi", "1", "-ntomp", "2"], directory)
    logs += _run(["gmx", "grompp", "-f", "nvt.mdp", "-c", "minim.gro", "-p", "system.top",
                  "-o", "nvt.tpr", "-maxwarn", "0"], directory)
    logs += _run(["gmx", "mdrun", "-deffnm", "nvt", "-ntmpi", "1", "-ntomp", "2"], directory)
    combined = logs + (directory / "nvt.log").read_text(encoding="utf-8")
    markers = [name for name, pattern in {"LINCS": r"LINCS WARNING", "NaN": r"\bnan\b",
               "non-finite": r"non-finite"}.items()
               if re.search(pattern, combined, flags=re.IGNORECASE)]
    reached = "Statistics over 1001 steps" in combined
    return {"grompp_maxwarn": 0, "nvt_reached_final_step": reached,
            "failure_markers": markers, "passed": reached and not markers}


def _openmm_smoke(directory: Path, python: str) -> dict[str, object]:
    output = _run([python, "run_openmm.py", "--platform", "CPU", "--minimize-iterations", "100",
                   "--steps", "1000", "--timestep-fs", "1.0", "--report-interval", "500"],
                  directory)
    match = re.search(r"Final potential energy: (\S+) kJ/mol", output)
    energy = float(match.group(1)) if match else float("nan")
    return {"steps": 1000, "timestep_fs": 1.0, "final_potential_energy_kj_mol": energy,
            "passed": bool(np.isfinite(energy))}


def _lammps_smoke(directory: Path) -> dict[str, object]:
    log = _run(["lmp_serial", "-in", "in.cellulose", "-log", "none"], directory)
    finished = "Loop time" in log
    errors = [line for line in log.splitlines() if line.startswith("ERROR")]
    return {"input": "in.cellulose (run 0)", "completed": finished, "errors": errors,
            "passed": finished and not errors}


def validate_case(
    kind: str, allomorph: str, boundary: str, python: str, gmx_double: str
) -> dict[str, object]:
    crystal = _crystal(kind, allomorph, boundary)
    periodic = boundary == "periodic"
    shape = "bulk" if kind == "bulk" else kind
    dispersion = dispersion_correction("opls-aa", shape)
    result: dict[str, object] = {"kind": kind, "allomorph": allomorph, "chain_boundary": boundary}
    potentials: dict[str, dict[str, float]] = {}
    with tempfile.TemporaryDirectory(prefix=f"opls_{kind}_{allomorph}_") as temporary:
        for variant in VARIANTS:
            build = build_opls_system(crystal, variant, boundary)  # type: ignore[arg-type]
            point = Path(temporary) / variant / "point"
            cross = _cross_engine(build, point, python, periodic, gmx_double)
            gromacs: dict[str, float] = cross["gromacs_terms"]  # type: ignore[assignment]
            potentials[variant] = gromacs
            recomputed: dict[str, dict[str, object]] = {}
            ours = _numpy_terms(point, build.mapped.chain_count, periodic)
            for term in NUMPY_GROMACS_TERMS:
                value = ours[term]
                difference = _relative(gromacs[term], value)
                recomputed[term] = {
                    "gromacs_kj_mol": gromacs[term], "numpy_kj_mol": value,
                    "relative_difference": difference,
                    "passed": difference <= TERM_RELATIVE_TOLERANCE,
                }

            # Shipped outputs, as the CLI writes them, for the runtime checks.
            run = Path(temporary) / variant / "run"
            write_opls_system(build, run)
            write_gromacs_run_inputs(run, force_field=variant, periodic_chains=periodic,
                                     dispersion=dispersion)
            write_opls_openmm(build, run, dispersion_correction=dispersion,
                              python_executable=python)
            write_opls_lammps(build, run, dispersion_correction=dispersion)
            smoke = {
                "gromacs_smoke_test": _gromacs_smoke(run),
                "openmm_smoke_test": _openmm_smoke(run, python),
                "lammps_smoke_test": _lammps_smoke(run),
            }
            entry: dict[str, object] = {
                "atom_count": build.mapped.atom_count,
                "chain_count": build.mapped.chain_count,
                "dispersion_correction_in_shipped_outputs": dispersion,
                "numpy_recomputation": recomputed,
                "cross_engine": cross["terms"],
                **smoke,
            }
            checks = [item["passed"] for item in recomputed.values()]
            checks += [
                item["passed"] for item in cross["terms"].values()  # type: ignore[attr-defined]
                if isinstance(item, dict)
            ]
            checks += [item["passed"] for item in smoke.values()]
            if variant == "opls-aa-sei":
                predicted = _sei_prediction(point, build.mapped.chain_count, periodic)
                nonbonded = ("LJ-14", "Coulomb-14", "LJ-(SR)", "Coulomb-(SR)", "Coul.-recip.")
                observed = sum(
                    potentials["opls-aa-sei"][t] - potentials["opls-aa"][t] for t in nonbonded
                )
                difference = abs(observed - predicted) / max(abs(predicted), 1.0)
                entry["sei_scaling_check"] = {
                    "predicted_nonbonded_change_kj_mol": predicted,
                    "gromacs_nonbonded_change_kj_mol": observed,
                    "relative_difference": difference,
                    "passed": difference <= PREDICTION_RELATIVE_TOLERANCE,
                }
                checks.append(difference <= PREDICTION_RELATIVE_TOLERANCE)
            entry["passed"] = bool(all(checks))
            result[variant] = entry
    result["passed"] = all(result[v]["passed"] for v in VARIANTS)  # type: ignore[index]
    return result


CASES = [
    *[("slab", a, b) for a in SLAB_FACES for b in ("finite", "periodic")],
    *[("bulk", a, "periodic") for a in BULK_CELLS],
    *[("fibril", a, b) for a in ("i_alpha", "i_beta", "ii") for b in ("finite", "periodic")],
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", help="kind:allomorph:boundary")
    parser.add_argument("--openmm-python", default=sys.executable,
                        help="Python interpreter with OpenMM")
    parser.add_argument("--gmx-double", default="gmx_d",
                        help="double-precision GROMACS used for the energy reference")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "validation/results/opls_cellulose_matrix.json")
    args = parser.parse_args()
    for tool in ("gmx", args.gmx_double, "lmp_serial"):
        if shutil.which(tool) is None:
            raise SystemExit(f"{tool} is required")
    selected = [tuple(case.split(":")) for case in args.case] if args.case else CASES
    cases = {}
    for kind, allomorph, boundary in selected:
        name = f"{kind}_{allomorph}_{boundary}"
        cases[name] = validate_case(
            kind, allomorph, boundary, args.openmm_python, args.gmx_double
        )
        print(f"{name}: passed={cases[name]['passed']}", file=sys.stderr, flush=True)
    openmm_version = _run([args.openmm_python, "-c", "import openmm; print(openmm.__version__)"],
                          ROOT).strip()
    lammps_banner = _run(["lmp_serial", "-h"], ROOT).splitlines()
    transcription = transcription_check()
    result = {
        "matrix": "experimental OPLS-AA and OPLS-AA-SEI cellulose, GROMACS/OpenMM/LAMMPS",
        "provenance": provenance(("numpy",)),
        "engines": {
            "gromacs_mixed_runtime": _gmx_version("gmx"),
            "gromacs_double_energy_reference": _gmx_version(args.gmx_double),
            "openmm": openmm_version,
            "lammps": next((line.strip() for line in lammps_banner if "Simulator" in line), ""),
        },
        "transcription_check": transcription,
        "term_relative_tolerance": TERM_RELATIVE_TOLERANCE,
        "total_relative_tolerance": TOTAL_RELATIVE_TOLERANCE,
        "cases": cases,
        "passed": bool(transcription["passed"] and all(c["passed"] for c in cases.values())),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, default=lambda value: value.item()) + "\n",
        encoding="utf-8",
    )
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
