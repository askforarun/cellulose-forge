"""Experimental CHARMM Drude polarizable cellulose output for NAMD."""

from __future__ import annotations

import hashlib
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cellulose_forge.crystallography import FloatArray
from cellulose_forge.geometry import BulkCrystal
from cellulose_forge.polymer_mapping import chain_template, map_crystal_to_chain_template

DRUDE_RELEASE = "2023a"
DRUDE_MAIN_FILE = "toppar_drude_main_protein_2023a.str"
DRUDE_CARBOHYDRATE_FILE = "toppar_drude_carbohydrate_2023a.str"
_PHYSICAL_ATOMS = (
    "C1",
    "H1",
    "C2",
    "H2",
    "C3",
    "H3",
    "C4",
    "H4",
    "C5",
    "H5",
    "C6",
    "H61",
    "H62",
    "O2",
    "O3",
    "O4",
    "O5",
    "O6",
    "HO2",
    "HO3",
    "HO6",
)


@dataclass(frozen=True)
class DrudeSystemFiles:
    """Files written for one experimental NAMD Drude system."""

    psf: Path
    pdb: Path
    configuration: Path
    build_script: Path
    main_stream: Path
    carbohydrate_stream: Path
    source_sha256: dict[str, str]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_drude_toppar_directory(directory: str | Path) -> tuple[Path, Path]:
    """Return the two official 2023a stream files required by the backend."""
    root = Path(directory).expanduser().resolve()
    main = root / DRUDE_MAIN_FILE
    carbohydrate = root / DRUDE_CARBOHYDRATE_FILE
    missing = [path.name for path in (main, carbohydrate) if not path.is_file()]
    if missing:
        raise ValueError(
            f"CHARMM Drude directory {root} is missing {', '.join(missing)}; "
            "use the official MacKerell drude_toppar 2023 release"
        )
    main_text = main.read_text(encoding="utf-8", errors="replace")
    carbohydrate_text = carbohydrate.read_text(encoding="utf-8", errors="replace")
    if "MASS    -1 DRUD" not in main_text or "MASS    -1 LPD" not in main_text:
        raise ValueError(f"{main.name} does not define the required DRUD and LPD atom types")
    if "RESI BGLC" not in carbohydrate_text or "PRES 14bb" not in carbohydrate_text:
        raise ValueError(f"{carbohydrate.name} lacks BGLC or the 14bb cellulose linkage")
    return main, carbohydrate


def _finite_template(degree_polymerization: int) -> tuple[tuple[str, int], ...]:
    atoms: list[tuple[str, int]] = []
    for residue in range(1, degree_polymerization + 1):
        names = list(_PHYSICAL_ATOMS)
        if residue == 1:
            names.extend(("O1", "HO1"))
        if residue == degree_polymerization:
            names.append("HO4")
        atoms.extend((name, residue) for name in names)
    return tuple(atoms)


def _render_coordinate_pdb(system: object) -> str:
    atom_names = system.atom_names
    residue_numbers = system.residue_numbers
    chain_numbers = system.chain_numbers
    coordinates = system.cartesian_coordinates
    lines = ["REMARK physical coordinates for CHARMM Drude psfgen construction"]
    for serial, (name, residue, chain, xyz) in enumerate(
        zip(atom_names, residue_numbers, chain_numbers, coordinates, strict=True), start=1
    ):
        x, y, z = xyz
        segment = f"C{chain:03d}"
        lines.append(
            f"ATOM  {serial:5d} {name:>4s} BGLC A{residue:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00      {segment:4s}{name[0]:>2s}"
        )
    lines.append("END")
    return "\n".join(lines) + "\n"


def _render_psfgen_script(chain_count: int, degree_polymerization: int) -> str:
    segments = " ".join(f"C{chain:03d}" for chain in range(1, chain_count + 1))
    residues = "\n".join(
        f"    residue {residue} BGLC" for residue in range(1, degree_polymerization + 1)
    )
    patches = "\n".join(
        f"  patch 14bb $seg:{residue} $seg:{residue + 1}"
        for residue in range(1, degree_polymerization)
    )
    return f"""package require psfgen
topology {DRUDE_MAIN_FILE}
topology {DRUDE_CARBOHYDRATE_FILE}
foreach seg {{{segments}}} {{
  segment $seg {{
    first NONE
    last NONE
{residues}
  }}
{patches}
  coordpdb physical_coordinates.pdb $seg
}}
guesscoord
writepsf system.psf
writepdb structure_namd.pdb
quit
"""


def render_drude_namd_configuration(vectors: FloatArray) -> str:
    """Render a conservative minimization and 1 ps dual-thermostat smoke test."""
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
            f"parameters {DRUDE_MAIN_FILE}",
            f"parameters {DRUDE_CARBOHYDRATE_FILE}",
            "paraTypeCharmm on",
            "exclude scaled1-4",
            "1-4scaling 1.0",
            "switching on",
            "switchdist 10.0",
            "cutoff 12.0",
            "pairlistdist 14.0",
            "PME yes",
            "PMEGridSpacing 1.0",
            *vector_lines,
            "wrapAll off",
            "timestep 1.0",
            "rigidBonds none",
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
            "temperature 300",
            "outputName drude",
            "minimize 1000",
            "reinitvels 300",
            "run 1000",
            "",
        ]
    )


def write_drude_namd_system(
    crystal: BulkCrystal,
    toppar_directory: str | Path,
    output_directory: str | Path,
    *,
    vmd_executable: str | Path,
) -> DrudeSystemFiles:
    """Build a finite-chain CHARMM Drude PSF/PDB and NAMD smoke-test input."""
    main_source, carbohydrate_source = validate_drude_toppar_directory(toppar_directory)
    degree_polymerization = len(set(crystal.residue_numbers))
    mapped = map_crystal_to_chain_template(
        crystal,
        chain_template(_finite_template(degree_polymerization)),
        "finite",
        "BGLC",
    )
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=True)
    main_stream = directory / DRUDE_MAIN_FILE
    carbohydrate_stream = directory / DRUDE_CARBOHYDRATE_FILE
    shutil.copy2(main_source, main_stream)
    shutil.copy2(carbohydrate_source, carbohydrate_stream)
    physical = directory / "physical_coordinates.pdb"
    physical.write_text(_render_coordinate_pdb(mapped), encoding="utf-8")
    build_script = directory / "build_drude_psf.tcl"
    build_script.write_text(
        _render_psfgen_script(mapped.chain_count, degree_polymerization), encoding="utf-8"
    )
    executable = Path(vmd_executable).expanduser()
    if not executable.is_file():
        raise ValueError(f"VMD executable not found: {executable}")
    result = subprocess.run(
        [str(executable), "-dispdev", "text", "-e", build_script.name],
        cwd=directory,
        capture_output=True,
        text=True,
        check=False,
    )
    (directory / "psfgen.log").write_text(result.stdout + result.stderr, encoding="utf-8")
    psf = directory / "system.psf"
    pdb = directory / "structure_namd.pdb"
    if result.returncode or not psf.is_file() or not pdb.is_file():
        raise ValueError("VMD/psfgen failed to construct the Drude system; see psfgen.log")
    psf_text = psf.read_text(encoding="utf-8", errors="replace")
    required = ("PSF EXT DRUDE", "!NUMLP", "!NANISO")
    missing = [marker for marker in required if marker not in psf_text]
    if missing:
        raise ValueError(
            f"psfgen output is not a complete Drude PSF (missing {', '.join(missing)})"
        )
    configuration = directory / "run_drude_namd.conf"
    configuration.write_text(
        render_drude_namd_configuration(mapped.basis_vectors), encoding="utf-8"
    )
    return DrudeSystemFiles(
        psf=psf,
        pdb=pdb,
        configuration=configuration,
        build_script=build_script,
        main_stream=main_stream,
        carbohydrate_stream=carbohydrate_stream,
        source_sha256={
            main_source.name: _sha256(main_source),
            carbohydrate_source.name: _sha256(carbohydrate_source),
        },
    )
