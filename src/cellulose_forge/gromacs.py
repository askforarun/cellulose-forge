"""GROMACS output for force-field-mapped cellulose systems."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np

from cellulose_forge.crystallography import FloatArray
from cellulose_forge.glycam_mapping import GlycamMappedSystem
from cellulose_forge.restraints import position_restraints_section

_ANGSTROM_TO_NANOMETER = 0.1
_GRO_INDEX_MODULUS = 100_000
_ZERO_BOX_COMPONENT_TOLERANCE = 1.0e-12


def _gro_box_line(basis_vectors: FloatArray) -> str:
    vectors = basis_vectors * _ANGSTROM_TO_NANOMETER
    a, b, c = vectors
    values: Sequence[np.float64] = (
        a[0], b[1], c[2], a[1], a[2], b[0], b[2], c[0], c[1]
    )
    if all(abs(value) < _ZERO_BOX_COMPONENT_TOLERANCE for value in values[3:]):
        values = values[:3]
    return "".join(f"{value:10.5f}" for value in values)


def write_gro(system: GlycamMappedSystem, path: str | Path) -> None:
    """Write mapped coordinates in GROMACS GRO atom order."""
    lines = ["CelluloseForge force-field-mapped cellulose", f"{system.atom_count:5d}"]
    for atom_index, (atom_name, residue_name, residue_number, coordinates) in enumerate(
        zip(
            system.atom_names,
            system.residue_names,
            system.residue_numbers,
            system.cartesian_coordinates,
            strict=True,
        ),
        start=1,
    ):
        x, y, z = coordinates * _ANGSTROM_TO_NANOMETER
        gro_atom_index = atom_index % _GRO_INDEX_MODULUS
        gro_residue_number = residue_number % _GRO_INDEX_MODULUS
        lines.append(
            f"{gro_residue_number:5d}{residue_name:<5.5s}{atom_name:>5.5s}"
            f"{gro_atom_index:5d}{x:8.3f}{y:8.3f}{z:8.3f}"
        )
    lines.append(_gro_box_line(system.basis_vectors))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_repeated_chain_topology(  # noqa: PLR0913
    template_topology: str | Path,
    output_path: str | Path,
    *,
    chain_count: int,
    restrained_chain_count: int = 0,
    restraint_force: tuple[float, float, float] | None = None,
    restraint_define: str = "POSRES_CELLULOSE",
) -> None:
    """Write a standalone topology with the requested repeated-chain count."""
    if chain_count < 1:
        msg = "GROMACS topology chain count must be positive"
        raise ValueError(msg)
    lines = Path(template_topology).read_text(encoding="utf-8").splitlines()
    free_chain_count = chain_count - restrained_chain_count
    if restrained_chain_count < 0 or free_chain_count < 0:
        msg = "restrained chain count must be between zero and the total chain count"
        raise ValueError(msg)
    if restrained_chain_count and restraint_force is None:
        msg = "restraint force constants are required for restrained chains"
        raise ValueError(msg)

    template_molecule_name = _template_molecule_name(lines)
    restrained_molecule_name = f"{template_molecule_name}_RESTRAINED"
    if restrained_chain_count:
        assert restraint_force is not None
        _insert_restrained_glycam_molecule(
            lines,
            restraint_force,
            restraint_define,
            template_molecule_name,
            restrained_molecule_name,
        )
    in_molecules = False
    replaced = False
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            in_molecules = stripped[1:-1].strip() == "molecules"
            continue
        if in_molecules and stripped and not stripped.startswith(";"):
            fields = stripped.split()
            if fields[0] != template_molecule_name:
                msg = f"unexpected GLYCAM molecule name in template: {fields[0]}"
                raise ValueError(msg)
            molecule_rows: list[str] = []
            if free_chain_count:
                molecule_rows.append(f"{template_molecule_name}  {free_chain_count}")
            if restrained_chain_count:
                molecule_rows.append(
                    f"{restrained_molecule_name}  {restrained_chain_count}"
                )
            lines[index : index + 1] = molecule_rows
            replaced = True
            break
    if not replaced:
        msg = "template topology has no [ molecules ] data row"
        raise ValueError(msg)
    Path(output_path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def _template_molecule_name(lines: list[str]) -> str:
    section = ""
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped[1:-1].strip().lower()
            continue
        if section == "molecules" and stripped and not stripped.startswith(";"):
            return stripped.split()[0]
    raise ValueError("template topology has no [ molecules ] data row")


def _insert_restrained_glycam_molecule(
    lines: list[str],
    restraint_force: tuple[float, float, float],
    restraint_define: str,
    free_molecule_name: str,
    restrained_molecule_name: str,
) -> None:
    """Duplicate the one-chain molecule definition and attach position restraints."""
    molecule_start = next(
        index for index, line in enumerate(lines) if line.strip().lower() == "[ moleculetype ]"
    )
    system_index = next(
        index for index, line in enumerate(lines) if line.strip().lower() == "[ system ]"
    )
    restrained_block = list(lines[molecule_start:system_index])
    current_section = ""
    atom_names: list[str] = []
    molecule_name_replaced = False
    for index, line in enumerate(restrained_block):
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            current_section = stripped[1:-1].strip().lower()
            continue
        if not stripped or stripped.startswith(";"):
            continue
        if current_section == "moleculetype" and not molecule_name_replaced:
            fields = stripped.split()
            if fields[0] != free_molecule_name:
                msg = f"unexpected GLYCAM molecule name in template: {fields[0]}"
                raise ValueError(msg)
            fields[0] = restrained_molecule_name
            restrained_block[index] = "  ".join(fields)
            molecule_name_replaced = True
        elif current_section == "atoms" and not stripped.startswith("#"):
            fields = stripped.split()
            try:
                atom_name = fields[4]
            except IndexError as error:
                msg = "malformed GROMACS [ atoms ] row in GLYCAM topology"
                raise ValueError(msg) from error
            atom_names.append(atom_name)
    restrained_block.extend(
        position_restraints_section(
            atom_names,
            restraint_force,
            define=restraint_define,
        ).splitlines()
    )
    lines[system_index:system_index] = restrained_block


def write_glycam_gromacs_system(  # noqa: PLR0913, PLR0917
    system: GlycamMappedSystem,
    template_topology: str | Path,
    output_directory: str | Path,
    restrained_chain_count: int = 0,
    restraint_force: tuple[float, float, float] | None = None,
    restraint_define: str = "POSRES_CELLULOSE",
) -> tuple[Path, Path]:
    """Write a mapped multi-chain GLYCAM system as GRO and standalone TOP."""
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=True)
    gro_path = directory / "structure.gro"
    topology_path = directory / "system.top"
    write_gro(system, gro_path)
    write_repeated_chain_topology(
        template_topology,
        topology_path,
        chain_count=system.chain_count,
        restrained_chain_count=restrained_chain_count,
        restraint_force=restraint_force,
        restraint_define=restraint_define,
    )
    return gro_path, topology_path
