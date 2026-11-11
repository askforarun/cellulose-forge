"""Native Amber output for complete GLYCAM cellulose systems."""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from cellulose_forge.crystallography import FloatArray
from cellulose_forge.glycam import find_parmed_python
from cellulose_forge.glycam_mapping import GlycamMappedSystem


def _amber_box(basis_vectors: FloatArray) -> list[float]:
    """Return Amber lengths and angles for three row-wise box vectors."""
    vectors = np.asarray(basis_vectors, dtype=np.float64)
    lengths = np.linalg.norm(vectors, axis=1)
    if np.any(lengths <= 0.0):
        raise ValueError("Amber output requires three nonzero box vectors")

    def angle(first: int, second: int) -> float:
        cosine = float(
            np.dot(vectors[first], vectors[second])
            / (lengths[first] * lengths[second])
        )
        return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))

    alpha = angle(1, 2)
    beta = angle(0, 2)
    gamma = angle(0, 1)
    return [*map(float, lengths), alpha, beta, gamma]


def _resolve_python(executable: str | Path | None) -> str:
    requested = str(executable or find_parmed_python())
    resolved = shutil.which(requested)
    if resolved is not None:
        return resolved
    path = Path(requested)
    if path.is_file():
        return str(path)
    raise FileNotFoundError(f"Python executable for ParmEd was not found: {requested}")


def write_glycam_amber_system(
    system: GlycamMappedSystem,
    template_prmtop: str | Path,
    template_inpcrd: str | Path,
    output_directory: str | Path,
    *,
    python_executable: str | Path | None = None,
) -> tuple[Path, Path, dict[str, object]]:
    """Write a complete multi-chain GLYCAM system as Amber topology/restart files."""
    template_topology = Path(template_prmtop).resolve()
    template_coordinates = Path(template_inpcrd).resolve()
    for required in (template_topology, template_coordinates):
        if not required.is_file():
            raise FileNotFoundError(f"required Amber chain file was not found: {required}")

    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=True)
    topology_path = (directory / "structure.prmtop").resolve()
    coordinate_path = (directory / "structure.inpcrd").resolve()
    payload = {
        "coordinates": np.asarray(system.cartesian_coordinates, dtype=float).tolist(),
        "box": _amber_box(np.asarray(system.basis_vectors)),
        "chain_count": system.chain_count,
        "atom_count": system.atom_count,
    }
    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".json",
        prefix="cellulose_forge_amber_",
        dir=directory,
        encoding="utf-8",
        delete=False,
    ) as handle:
        json.dump(payload, handle)
        payload_path = Path(handle.name)

    conversion_code = """
import copy
import json
import sys

import numpy as np
import parmed as pmd

loaded = pmd.load_file(sys.argv[1], xyz=sys.argv[2])
template = loaded.copy(pmd.Structure)
with open(sys.argv[3], encoding="utf-8") as handle:
    payload = json.load(handle)
combined = copy.deepcopy(template)
for _ in range(1, int(payload["chain_count"])):
    combined += copy.deepcopy(template)
combined.coordinates = np.asarray(payload["coordinates"], dtype=float)
combined.box = np.asarray(payload["box"], dtype=float)
if len(combined.atoms) != int(payload["atom_count"]):
    raise ValueError(
        f"assembled Amber atom count {len(combined.atoms)} does not match "
        f"mapped count {payload['atom_count']}"
    )
combined.save(sys.argv[4], overwrite=True)
combined.save(sys.argv[5], format="rst7", overwrite=True)
written = pmd.load_file(sys.argv[4], xyz=sys.argv[5])
print(json.dumps({
    "atom_count": len(written.atoms),
    "residue_count": len(written.residues),
    "total_charge_e": float(sum(atom.charge for atom in written.atoms)),
    "box": [float(value) for value in written.box],
}))
"""
    executable = _resolve_python(python_executable)
    try:
        result = subprocess.run(
            [
                executable,
                "-c",
                conversion_code,
                str(template_topology),
                str(template_coordinates),
                str(payload_path.resolve()),
                str(topology_path),
                str(coordinate_path),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
    finally:
        payload_path.unlink(missing_ok=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"ParmEd failed to assemble the complete Amber system "
            f"(exit code {result.returncode}):\n{result.stdout}\n{result.stderr}"
        )
    if not topology_path.is_file() or not coordinate_path.is_file():
        raise RuntimeError("ParmEd did not create the complete Amber topology and coordinates")
    try:
        metadata = json.loads(result.stdout.strip().splitlines()[-1])
    except (IndexError, json.JSONDecodeError) as error:
        raise RuntimeError(
            f"could not parse ParmEd Amber validation output: {result.stdout}"
        ) from error
    if metadata["atom_count"] != system.atom_count:
        raise ValueError("written Amber topology atom count does not match mapped coordinates")
    metadata.update(
        {
            "topology": topology_path.name,
            "coordinates": coordinate_path.name,
            "chain_count": system.chain_count,
            "parmed_python_executable": executable,
        }
    )
    return topology_path, coordinate_path, metadata
