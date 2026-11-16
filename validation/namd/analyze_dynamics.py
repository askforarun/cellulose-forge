"""Validate a completed NAMD minimization and 1 ps NVT smoke test."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

_EXPECTED_FINAL_STEP = 1500
_DYNAMICS_FIRST_STEP = 1001
_EXPECTED_DYNAMICS_ROWS = 500
_MIN_BOND_ANGSTROM = 0.85
_MAX_BOND_ANGSTROM = 1.75
_MAX_HEAVY_RMSD_ANGSTROM = 2.0
_MIN_MEAN_TEMPERATURE_K = 150.0
_MAX_MEAN_TEMPERATURE_K = 450.0
_XSC_MINIMUM_FIELD_COUNT = 10
_HYDROGEN_MASS_THRESHOLD = 2.0


def _parse_energy_rows(path: Path) -> list[dict[str, float]]:
    fields: list[str] | None = None
    rows: list[dict[str, float]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ETITLE:"):
            fields = line.split()[1:]
        elif line.startswith("ENERGY:"):
            if fields is None:
                raise ValueError("NAMD dynamics log contains ENERGY before ETITLE")
            values = [float(value) for value in line.split()[1:]]
            if len(values) != len(fields):
                raise ValueError("NAMD dynamics ENERGY field count differs from ETITLE")
            rows.append(dict(zip(fields, values, strict=True)))
    if not rows:
        raise ValueError("NAMD dynamics log contains no ENERGY rows")
    return rows


def _parse_psf(path: Path) -> tuple[NDArray[np.float64], list[tuple[int, int]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    atom_row = next(index for index, line in enumerate(lines) if "!NATOM" in line)
    atom_count = int(lines[atom_row].split()[0])
    masses = np.asarray(
        [float(line.split()[7]) for line in lines[atom_row + 1 : atom_row + atom_count + 1]]
    )
    bond_row = next(index for index, line in enumerate(lines) if "!NBOND" in line)
    bond_count = int(lines[bond_row].split()[0])
    values: list[int] = []
    row = bond_row + 1
    while len(values) < 2 * bond_count:
        values.extend(int(value) - 1 for value in lines[row].split())
        row += 1
    bonds = list(zip(values[::2], values[1::2], strict=True))
    return masses, bonds


def _parse_pdb_coordinates(path: Path) -> NDArray[np.float64]:
    coordinates = [
        [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]
    if not coordinates:
        raise ValueError("NAMD final coordinate file contains no atoms")
    return np.asarray(coordinates, dtype=float)


def _parse_gro_coordinates(path: Path) -> NDArray[np.float64]:
    lines = path.read_text(encoding="utf-8").splitlines()
    atom_count = int(lines[1])
    return 10.0 * np.asarray(
        [
            [float(line[20:28]), float(line[28:36]), float(line[36:44])]
            for line in lines[2 : 2 + atom_count]
        ],
        dtype=float,
    )


def _parse_xsc_box(path: Path) -> NDArray[np.float64]:
    rows = [
        line.split()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    if not rows:
        raise ValueError("NAMD extended-system file contains no box row")
    values = rows[-1]
    if len(values) < _XSC_MINIMUM_FIELD_COUNT:
        raise ValueError("NAMD extended-system box row is incomplete")
    return np.asarray(
        [
            [float(value) for value in values[1:4]],
            [float(value) for value in values[4:7]],
            [float(value) for value in values[7:10]],
        ]
    )


def _minimum_image(
    displacements: NDArray[np.float64], vectors: NDArray[np.float64]
) -> NDArray[np.float64]:
    fractional = displacements @ np.linalg.inv(vectors)
    fractional -= np.rint(fractional)
    return fractional @ vectors


def analyze(build: Path) -> dict[str, object]:
    """Return pass/fail diagnostics for one generated NAMD production run."""
    log_path = build / "production_dynamics.log"
    log_text = log_path.read_text(encoding="utf-8", errors="replace")
    rows = _parse_energy_rows(log_path)
    dynamics = [row for row in rows if int(row["TS"]) >= _DYNAMICS_FIRST_STEP]
    masses, bonds = _parse_psf(build / "system.psf")
    initial = _parse_gro_coordinates(build / "single_point_coordinates.gro")
    final = _parse_pdb_coordinates(build / "namd.coor")
    vectors = _parse_xsc_box(build / "namd.xsc")
    if len(initial) != len(final) or len(final) != len(masses):
        raise ValueError("initial coordinates, final coordinates, and PSF atom counts differ")

    bond_vectors = _minimum_image(
        np.asarray([final[second] - final[first] for first, second in bonds]), vectors
    )
    bond_lengths = np.linalg.norm(bond_vectors, axis=1)
    displacements = _minimum_image(final - initial, vectors)
    heavy = masses > _HYDROGEN_MASS_THRESHOLD
    displacements -= displacements[heavy].mean(axis=0)
    heavy_rmsd = float(np.sqrt(np.mean(np.sum(displacements[heavy] ** 2, axis=1))))
    temperatures = np.asarray([row["TEMP"] for row in dynamics], dtype=float)
    final_row = rows[-1]
    finite_energies = all(math.isfinite(value) for value in final_row.values())
    completed = "End of program" in log_text
    criteria = {
        "completed_without_runtime_failure": completed,
        "final_step_is_1500": int(final_row["TS"]) == _EXPECTED_FINAL_STEP,
        "has_500_dynamics_energy_rows": len(dynamics) == _EXPECTED_DYNAMICS_ROWS,
        "final_energy_values_are_finite": finite_energies,
        "mean_temperature_is_150_to_450_k": (
            len(temperatures) > 0
            and _MIN_MEAN_TEMPERATURE_K <= float(temperatures.mean()) <= _MAX_MEAN_TEMPERATURE_K
        ),
        "all_bonds_are_0.85_to_1.75_angstrom": (
            float(bond_lengths.min()) >= _MIN_BOND_ANGSTROM
            and float(bond_lengths.max()) <= _MAX_BOND_ANGSTROM
        ),
        "heavy_atom_rmsd_below_2_angstrom": heavy_rmsd <= _MAX_HEAVY_RMSD_ANGSTROM,
    }
    return {
        "protocol": {
            "minimization_steps": 1000,
            "nvt_steps": 500,
            "timestep_fs": 2.0,
            "target_temperature_k": 300.0,
            "lj_switch_angstrom": [10.0, 12.0],
            "pme_tolerance": 1.0e-6,
        },
        "final_step": int(final_row["TS"]),
        "dynamics_energy_rows": len(dynamics),
        "mean_temperature_k": float(temperatures.mean()),
        "temperature_range_k": [float(temperatures.min()), float(temperatures.max())],
        "bond_count": len(bonds),
        "bond_length_range_angstrom": [
            float(bond_lengths.min()),
            float(bond_lengths.max()),
        ],
        "heavy_atom_rmsd_angstrom": heavy_rmsd,
        "criteria": criteria,
        "passed": all(criteria.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = analyze(arguments.build.resolve())
    rendered = json.dumps(result, indent=2) + "\n"
    if arguments.output is None:
        print(rendered, end="")
    else:
        arguments.output.write_text(rendered, encoding="utf-8")
    if result["passed"] is not True:
        raise SystemExit("NAMD production dynamics validation failed")


if __name__ == "__main__":
    main()
