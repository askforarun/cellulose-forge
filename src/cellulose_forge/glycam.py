"""GLYCAM finite-cellulose-chain reference topology generation."""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

GLYCAM_LEAPRC = "leaprc.GLYCAM_06j-1"
GLYCAM_CHARGE_SCALE = 18.2223
_SCALE_TOLERANCE = 1.0e-8
_GROMACS_MOLECULE_NAME = "CELLULOSE_CHAIN"
_MINIMUM_PERIODIC_DEGREE = 2


@dataclass(frozen=True)
class GlycamChainSpec:
    """GLYCAM residue sequence for one finite beta(1->4) cellulose chain."""

    degree_polymerization: int
    leaprc: str = GLYCAM_LEAPRC

    def __post_init__(self) -> None:
        if (
            not isinstance(self.degree_polymerization, int)
            or isinstance(self.degree_polymerization, bool)
            or self.degree_polymerization < 1
        ):
            msg = "degree of polymerization must be a positive integer"
            raise ValueError(msg)

    @property
    def tleap_residues(self) -> tuple[str, ...]:
        """Return residues in tleap's reducing-to-nonreducing build order."""
        return ("ROH", *("4GB",) * (self.degree_polymerization - 1), "0GB")

    @property
    def chemical_residues(self) -> tuple[str, ...]:
        """Return residues in nonreducing-to-reducing chemical display order."""
        return tuple(reversed(self.tleap_residues))


@dataclass(frozen=True)
class GlycamPeriodicChainSpec:
    """GLYCAM residue sequence for one covalently periodic cellulose chain."""

    degree_polymerization: int
    leaprc: str = GLYCAM_LEAPRC

    def __post_init__(self) -> None:
        if (
            not isinstance(self.degree_polymerization, int)
            or isinstance(self.degree_polymerization, bool)
            or self.degree_polymerization < _MINIMUM_PERIODIC_DEGREE
        ):
            msg = "periodic degree of polymerization must be an integer of at least two"
            raise ValueError(msg)

    @property
    def tleap_residues(self) -> tuple[str, ...]:
        """Return the neutral internal residues used for a periodic chain."""
        return ("4GB",) * self.degree_polymerization


@dataclass(frozen=True)
class GlycamChainFiles:
    """Paths generated for one GLYCAM chain reference topology."""

    directory: Path
    prefix: str

    @property
    def leap(self) -> Path:
        return self.directory / f"{self.prefix}.leap"

    @property
    def pdb(self) -> Path:
        return self.directory / f"{self.prefix}.pdb"

    @property
    def prmtop(self) -> Path:
        return self.directory / f"{self.prefix}.prmtop"

    @property
    def inpcrd(self) -> Path:
        return self.directory / f"{self.prefix}.inpcrd"

    @property
    def metadata(self) -> Path:
        return self.directory / f"{self.prefix}.json"

    @property
    def gromacs_topology(self) -> Path:
        return self.directory / f"{self.prefix}.top"

    @property
    def gromacs_coordinates(self) -> Path:
        return self.directory / f"{self.prefix}.gro"


def write_tleap_chain_input(
    spec: GlycamChainSpec,
    directory: str | Path,
    *,
    prefix: str = "cellulose_chain",
) -> GlycamChainFiles:
    """Write a tleap input file for one finite GLYCAM cellulose chain."""
    output_directory = Path(directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    files = GlycamChainFiles(output_directory, prefix)
    sequence = " ".join(spec.tleap_residues)
    files.leap.write_text(
        "\n".join(
            (
                f"source {spec.leaprc}",
                f"chain = sequence {{ {sequence} }}",
                "check chain",
                f"savepdb chain {files.pdb.name}",
                f"saveamberparm chain {files.prmtop.name} {files.inpcrd.name}",
                "quit",
                "",
            )
        ),
        encoding="utf-8",
    )
    return files


def write_tleap_periodic_chain_input(
    spec: GlycamPeriodicChainSpec,
    directory: str | Path,
    *,
    prefix: str = "cellulose_periodic_chain",
) -> GlycamChainFiles:
    """Write tleap input that closes a GLYCAM chain across its boundary."""
    output_directory = Path(directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    files = GlycamChainFiles(output_directory, prefix)
    sequence = " ".join(spec.tleap_residues)
    files.leap.write_text(
        "\n".join(
            (
                f"source {spec.leaprc}",
                f"chain = sequence {{ {sequence} }}",
                f"bond chain.{spec.degree_polymerization}.O4 chain.1.C1",
                "check chain",
                f"savepdb chain {files.pdb.name}",
                f"saveamberparm chain {files.prmtop.name} {files.inpcrd.name}",
                "quit",
                "",
            )
        ),
        encoding="utf-8",
    )
    return files


def _prmtop_flag_tokens(path: str | Path, flag_name: str) -> list[str]:
    """Read whitespace-delimited values from one AMBER topology flag block."""
    values: list[str] = []
    in_target = False
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("%FLAG "):
            if in_target:
                break
            in_target = line.split(maxsplit=1)[1].strip() == flag_name
            continue
        if in_target and not line.startswith("%FORMAT"):
            values.extend(line.split())
    if not values:
        msg = f"could not read {flag_name} from {path}"
        raise ValueError(msg)
    return values


def validate_glycam_prmtop(
    path: str | Path,
    spec: GlycamChainSpec,
) -> dict[str, object]:
    """Validate core GLYCAM chain invariants in an AMBER topology."""
    pointer_values = [int(value) for value in _prmtop_flag_tokens(path, "POINTERS")]
    atom_count = pointer_values[0]
    residue_labels = _prmtop_flag_tokens(path, "RESIDUE_LABEL")
    expected_residues = list(spec.tleap_residues)
    if residue_labels != expected_residues:
        msg = f"GLYCAM residue sequence mismatch: {residue_labels} != {expected_residues}"
        raise ValueError(msg)

    charges = [float(value) / GLYCAM_CHARGE_SCALE for value in _prmtop_flag_tokens(path, "CHARGE")]
    total_charge = sum(charges)
    if abs(total_charge) > _SCALE_TOLERANCE:
        msg = f"finite GLYCAM cellulose chain is not neutral: charge={total_charge:.12f}"
        raise ValueError(msg)

    scaling: dict[str, list[float]] = {}
    for flag_name in ("SCEE_SCALE_FACTOR", "SCNB_SCALE_FACTOR"):
        values = sorted({float(value) for value in _prmtop_flag_tokens(path, flag_name)})
        if any(abs(value - 1.0) > _SCALE_TOLERANCE for value in values):
            msg = f"unexpected GLYCAM {flag_name}: {values}"
            raise ValueError(msg)
        scaling[flag_name] = values

    return {
        "atom_count": atom_count,
        "topology_residue_count": len(residue_labels),
        "glucose_residue_count": spec.degree_polymerization,
        "cap_residue_count": 1,
        "total_charge": total_charge,
        "scaling": scaling,
    }


def validate_periodic_glycam_prmtop(
    path: str | Path,
    spec: GlycamPeriodicChainSpec,
) -> dict[str, object]:
    """Validate atom counts, residue types, charge, and scaling of a periodic chain."""
    pointer_values = [int(value) for value in _prmtop_flag_tokens(path, "POINTERS")]
    atom_count = pointer_values[0]
    expected_atom_count = 21 * spec.degree_polymerization
    if atom_count != expected_atom_count:
        msg = f"periodic GLYCAM atom count mismatch: {atom_count} != {expected_atom_count}"
        raise ValueError(msg)
    residue_labels = _prmtop_flag_tokens(path, "RESIDUE_LABEL")
    if residue_labels != list(spec.tleap_residues):
        msg = "periodic GLYCAM chain must contain only internal 4GB residues"
        raise ValueError(msg)
    charges = [float(value) / GLYCAM_CHARGE_SCALE for value in _prmtop_flag_tokens(path, "CHARGE")]
    total_charge = sum(charges)
    if abs(total_charge) > _SCALE_TOLERANCE:
        msg = f"periodic GLYCAM cellulose chain is not neutral: charge={total_charge:.12f}"
        raise ValueError(msg)
    scaling: dict[str, list[float]] = {}
    for flag_name in ("SCEE_SCALE_FACTOR", "SCNB_SCALE_FACTOR"):
        values = sorted({float(value) for value in _prmtop_flag_tokens(path, flag_name)})
        if any(abs(value - 1.0) > _SCALE_TOLERANCE for value in values):
            msg = f"unexpected periodic GLYCAM {flag_name}: {values}"
            raise ValueError(msg)
        scaling[flag_name] = values
    return {
        "atom_count": atom_count,
        "residue_count": len(residue_labels),
        "total_charge": total_charge,
        "scaling": scaling,
    }


def run_tleap_periodic_chain_reference(
    spec: GlycamPeriodicChainSpec,
    directory: str | Path,
    *,
    prefix: str = "cellulose_periodic_chain",
    executable: str | Path = "tleap",
) -> dict[str, object]:
    """Run tleap and validate one covalently periodic GLYCAM reference chain."""
    executable_text = str(executable)
    resolved_executable = shutil.which(executable_text)
    if resolved_executable is None:
        executable_path = Path(executable_text)
        if not executable_path.is_file():
            msg = f"tleap executable was not found: {executable}"
            raise FileNotFoundError(msg)
        resolved_executable = str(executable_path)
    files = write_tleap_periodic_chain_input(spec, directory, prefix=prefix)
    result = subprocess.run(
        [resolved_executable, "-f", files.leap.name],
        cwd=files.directory,
        check=False,
        capture_output=True,
        text=True,
    )
    combined_output = f"{result.stdout}\n{result.stderr}"
    if result.returncode != 0 or "Errors = 1" in combined_output:
        msg = f"tleap periodic-chain generation failed:\n{combined_output}"
        raise RuntimeError(msg)
    validation = validate_periodic_glycam_prmtop(files.prmtop, spec)
    metadata: dict[str, object] = {
        "force_field": "GLYCAM",
        "chain_boundary": "periodic",
        "leaprc": spec.leaprc,
        "degree_polymerization": spec.degree_polymerization,
        "tleap_residues": list(spec.tleap_residues),
        "closure_bond": {
            "from": f"{spec.degree_polymerization}:O4",
            "to": "1:C1",
        },
        "tleap_executable": resolved_executable,
        "validation": validation,
    }
    files.metadata.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return metadata


def run_tleap_chain_reference(
    spec: GlycamChainSpec,
    directory: str | Path,
    *,
    prefix: str = "cellulose_chain",
    executable: str | Path = "tleap",
) -> dict[str, object]:
    """Run tleap and validate one finite-chain GLYCAM reference topology."""
    executable_text = str(executable)
    resolved_executable = shutil.which(executable_text)
    if resolved_executable is None:
        executable_path = Path(executable_text)
        if not executable_path.is_file():
            msg = f"tleap executable was not found: {executable}"
            raise FileNotFoundError(msg)
        resolved_executable = str(executable_path)

    files = write_tleap_chain_input(spec, directory, prefix=prefix)
    result = subprocess.run(
        [resolved_executable, "-f", files.leap.name],
        cwd=files.directory,
        check=False,
        capture_output=True,
        text=True,
    )
    combined_output = f"{result.stdout}\n{result.stderr}"
    failure_markers = ("FATAL", "Errors = 1", "Errors = 2", "Parameter file was not saved")
    if result.returncode != 0 or any(marker in combined_output for marker in failure_markers):
        msg = f"tleap failed with exit code {result.returncode}:\n{combined_output}"
        raise RuntimeError(msg)

    expected_outputs = (files.pdb, files.prmtop, files.inpcrd)
    missing = [str(path) for path in expected_outputs if not path.is_file()]
    if missing:
        msg = f"tleap did not create expected outputs: {missing}"
        raise RuntimeError(msg)

    validation = validate_glycam_prmtop(files.prmtop, spec)
    metadata: dict[str, object] = {
        "force_field": "GLYCAM",
        "leaprc": spec.leaprc,
        "degree_polymerization": spec.degree_polymerization,
        "tleap_residues": list(spec.tleap_residues),
        "chemical_residues": list(spec.chemical_residues),
        "tleap_executable": resolved_executable,
        "outputs": {
            "leap": files.leap.name,
            "pdb": files.pdb.name,
            "prmtop": files.prmtop.name,
            "inpcrd": files.inpcrd.name,
        },
        "validation": validation,
    }
    files.metadata.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return metadata


def convert_chain_reference_to_gromacs(
    files: GlycamChainFiles,
    *,
    python_executable: str | Path,
) -> dict[str, object]:
    """Convert a validated AMBER chain reference to standalone GROMACS files."""
    executable_text = str(python_executable)
    resolved_executable = shutil.which(executable_text)
    if resolved_executable is None:
        executable_path = Path(executable_text)
        if not executable_path.is_file():
            msg = f"Python executable for ParmEd was not found: {python_executable}"
            raise FileNotFoundError(msg)
        resolved_executable = str(executable_path)

    for required_path in (files.prmtop, files.inpcrd):
        if not required_path.is_file():
            msg = f"required AMBER chain file was not found: {required_path}"
            raise FileNotFoundError(msg)

    conversion_code = (
        "import parmed as pmd, sys; "
        "structure = pmd.load_file(sys.argv[1], xyz=sys.argv[2]); "
        "structure.save(sys.argv[3], overwrite=True); "
        "structure.save(sys.argv[4], overwrite=True)"
    )
    result = subprocess.run(
        [
            resolved_executable,
            "-c",
            conversion_code,
            files.prmtop.name,
            files.inpcrd.name,
            files.gromacs_topology.name,
            files.gromacs_coordinates.name,
        ],
        cwd=files.directory,
        check=False,
        capture_output=True,
        text=True,
    )
    combined_output = f"{result.stdout}\n{result.stderr}"
    if result.returncode != 0 or "Unhandled exception" in combined_output:
        msg = f"ParmEd failed with exit code {result.returncode}:\n{combined_output}"
        raise RuntimeError(msg)

    expected_outputs = (files.gromacs_topology, files.gromacs_coordinates)
    missing = [str(path) for path in expected_outputs if not path.is_file()]
    if missing:
        msg = f"ParmEd did not create expected GROMACS outputs: {missing}"
        raise RuntimeError(msg)

    topology_text = files.gromacs_topology.read_text(encoding="utf-8")
    topology_text = _replace_gromacs_molecule_name(
        topology_text,
        molecule_name=_GROMACS_MOLECULE_NAME,
    )
    files.gromacs_topology.write_text(topology_text, encoding="utf-8")
    defaults = _gromacs_section_rows(topology_text, "defaults")
    molecules = _gromacs_section_rows(topology_text, "molecules")
    if not defaults or defaults[0].split()[-2:] != ["1", "1"]:
        msg = "converted GROMACS topology does not preserve GLYCAM 1-4 scaling"
        raise ValueError(msg)
    if not molecules or molecules[-1].split()[-1] != "1":
        msg = "converted reference topology must contain exactly one chain"
        raise ValueError(msg)

    return {
        "gromacs_topology": files.gromacs_topology.name,
        "gromacs_coordinates": files.gromacs_coordinates.name,
        "molecule_count": 1,
        "molecule_name": _GROMACS_MOLECULE_NAME,
        "fudge_lj": 1.0,
        "fudge_qq": 1.0,
        "parmed_python_executable": resolved_executable,
    }


def _gromacs_section_rows(text: str, section_name: str) -> list[str]:
    """Return data rows from one named GROMACS topology section."""
    rows: list[str] = []
    in_section = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            in_section = stripped[1:-1].strip() == section_name
            continue
        if in_section and stripped and not stripped.startswith(";"):
            rows.append(stripped)
    return rows


def _replace_gromacs_molecule_name(text: str, *, molecule_name: str) -> str:
    """Give the generated one-chain molecule a stable, descriptive name."""
    lines = text.splitlines()
    for section_name in ("moleculetype", "molecules"):
        in_section = False
        for index, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith("[") and stripped.endswith("]"):
                in_section = stripped[1:-1].strip() == section_name
                continue
            if in_section and stripped and not stripped.startswith(";"):
                fields = stripped.split()
                fields[0] = molecule_name
                lines[index] = "  ".join(fields)
                break
    return "\n".join(lines) + "\n"


def metadata_as_json(metadata: dict[str, object]) -> str:
    """Serialize GLYCAM metadata consistently for reporting and tests."""
    return json.dumps(metadata, indent=2, sort_keys=True)


def spec_as_dict(spec: GlycamChainSpec) -> dict[str, object]:
    """Return a serializable chain specification."""
    return asdict(spec)
