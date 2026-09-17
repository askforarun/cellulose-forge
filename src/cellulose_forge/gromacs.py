"""GROMACS output for force-field-mapped cellulose systems."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np

from cellulose_forge.crystallography import FloatArray
from cellulose_forge.glycam_mapping import GlycamMappedSystem

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
    lines = ["CelluloseForge GLYCAM-mapped cellulose", f"{system.atom_count:5d}"]
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


def write_repeated_chain_topology(
    template_topology: str | Path,
    output_path: str | Path,
    *,
    chain_count: int,
) -> None:
    """Write a standalone topology with the requested repeated-chain count."""
    if chain_count < 1:
        msg = "GROMACS topology chain count must be positive"
        raise ValueError(msg)
    lines = Path(template_topology).read_text(encoding="utf-8").splitlines()
    in_molecules = False
    replaced = False
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            in_molecules = stripped[1:-1].strip() == "molecules"
            continue
        if in_molecules and stripped and not stripped.startswith(";"):
            fields = stripped.split()
            if fields[0] != "CELLULOSE_CHAIN":
                msg = f"unexpected GLYCAM molecule name in template: {fields[0]}"
                raise ValueError(msg)
            lines[index] = f"CELLULOSE_CHAIN  {chain_count}"
            replaced = True
            break
    if not replaced:
        msg = "template topology has no [ molecules ] data row"
        raise ValueError(msg)
    Path(output_path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_glycam_gromacs_system(
    system: GlycamMappedSystem,
    template_topology: str | Path,
    output_directory: str | Path,
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
    )
    return gro_path, topology_path
