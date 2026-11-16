"""Normalize LAMMPS single-point energies into canonical energy reports.

The adapter reads the LAMMPS input that produced a ``run 0`` log, so every
convention it applies is derived from the evaluated settings rather than
assumed:

- ``units real`` energies are converted from kcal/mol to kJ/mol;
- ``evdwl`` includes ``etail``, so the tail correction is moved out of ``vdw``
  into ``dispersion_correction``;
- ``ecoul`` and ``elong`` are electrostatics; CHARMM ``dihedral_style
  charmmfsw`` 1--4 energies are already tallied into ``evdwl`` and ``ecoul``;
- energies of fixes with ``fix_modify ... energy yes`` are restraints.

``lj/charmmfsw/coul/long`` systems additionally carry the three documented
CHARMM conventions removed by ``validation/validate_charmm36_engines.py``
before comparison with GROMACS/OpenMM force-switch references.  They are
recorded as visible corrections and never replace the raw values:

1. the CHARMM Coulomb constant used by LAMMPS (332.0716 kcal A/mol/e^2) versus
   the GROMACS constant (138.935458 kJ nm/mol/e^2);
2. the constant, force-free 1--4 Lennard-Jones offset added by
   ``dihedral_style charmmfsw``;
3. the constant potential shift between CHARMM ``vfswitch``
   (``lj/charmmfsw``) and GROMACS force switching.

Two contributions cannot be separated from a LAMMPS thermo row and must be
supplied by the caller when present: tabulated GLYCAM explicit 1--4 pairs
carried as pseudo-bonds (``bond_style hybrid ... table``), and periodic
impropers written as LAMMPS dihedrals.  Both are recorded as decompositions,
which move energy between categories without changing ``potential``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

from cellulose_forge.energy_drivers.report import EnergyReport, EnergyValue

KJ_PER_KCAL = 4.184
# LAMMPS CHARMM pair styles use the CHARMM Coulomb constant; GROMACS uses
# 138.935458 kJ nm/mol/e^2 = 332.0637 kcal A/mol/e^2.
LAMMPS_CHARMM_COULOMB = 332.0716
GROMACS_COULOMB = 138.935458 * 10.0 / KJ_PER_KCAL
CHARMM_COULOMB_RATIO = GROMACS_COULOMB / LAMMPS_CHARMM_COULOMB

_CHARMM_PAIR_STYLE = "lj/charmmfsw/coul/long"
_PLAIN_PAIR_STYLE = "lj/cut/coul/long"
_SUPPORTED_PAIR_STYLES = frozenset({_CHARMM_PAIR_STYLE, _PLAIN_PAIR_STYLE})
_EXCLUDED_BOND_DEPTH = 3
_CHARMM_SWITCH_ARGUMENTS = 2

# LAMMPS thermo keyword -> thermo header label.
_THERMO_LABELS: dict[str, str] = {
    "step": "Step",
    "elapsed": "Elapsed",
    "elaplong": "Elaplong",
    "dt": "Dt",
    "time": "Time",
    "cpu": "CPU",
    "atoms": "Atoms",
    "temp": "Temp",
    "press": "Press",
    "pe": "PotEng",
    "ke": "KinEng",
    "etotal": "TotEng",
    "enthalpy": "Enthalpy",
    "evdwl": "E_vdwl",
    "ecoul": "E_coul",
    "epair": "E_pair",
    "ebond": "E_bond",
    "eangle": "E_angle",
    "edihed": "E_dihed",
    "eimp": "E_impro",
    "emol": "E_mol",
    "elong": "E_long",
    "etail": "E_tail",
    "vol": "Volume",
    "density": "Density",
    "lx": "Lx",
    "ly": "Ly",
    "lz": "Lz",
    "xy": "Xy",
    "xz": "Xz",
    "yz": "Yz",
    "pxx": "Pxx",
    "pyy": "Pyy",
    "pzz": "Pzz",
    "pxy": "Pxy",
    "pxz": "Pxz",
    "pyz": "Pyz",
}
# Columns that are not part of the potential-energy decomposition.
_NON_POTENTIAL_LABELS = frozenset(
    {
        "Step",
        "Elapsed",
        "Elaplong",
        "Dt",
        "Time",
        "CPU",
        "Atoms",
        "Temp",
        "Press",
        "KinEng",
        "TotEng",
        "Enthalpy",
        "Volume",
        "Density",
        "Lx",
        "Ly",
        "Lz",
        "Xy",
        "Xz",
        "Yz",
        "Pxx",
        "Pyy",
        "Pzz",
        "Pxy",
        "Pxz",
        "Pyz",
    }
)
# Redundant LAMMPS subtotals; their components are classified individually.
_SUBTOTAL_LABELS = frozenset({"E_pair", "E_mol"})
_COMPONENT_LABELS = frozenset(
    {"PotEng", "E_bond", "E_angle", "E_dihed", "E_impro", "E_vdwl", "E_coul", "E_long", "E_tail"}
)
_ALWAYS_REQUIRED = ("PotEng", "E_bond", "E_angle", "E_dihed", "E_vdwl", "E_coul")


@dataclass(frozen=True)
class LammpsSettings:
    """Energy-relevant settings parsed from one LAMMPS input file."""

    units: str
    pair_style: str
    pair_arguments: tuple[str, ...]
    dihedral_style: str | None
    special_bonds: str | None
    bond_style: tuple[str, ...]
    improper_style: str | None
    tail_correction: bool
    kspace: bool
    thermo_normalized: bool
    thermo_labels: tuple[str, ...]
    energy_fixes: tuple[str, ...]

    @property
    def charmm(self) -> bool:
        """Whether the CHARMM force-switch conventions apply."""
        return self.pair_style == _CHARMM_PAIR_STYLE

    @property
    def explicit_pair_table(self) -> bool:
        """Whether tabulated explicit 1--4 pairs are carried as pseudo-bonds."""
        return (
            bool(self.bond_style) and self.bond_style[0] == "hybrid" and "table" in self.bond_style
        )

    def as_dict(self) -> dict[str, Any]:
        """JSON-compatible diagnostic view."""
        return {
            "units": self.units,
            "pair_style": " ".join((self.pair_style, *self.pair_arguments)),
            "dihedral_style": self.dihedral_style,
            "special_bonds": self.special_bonds,
            "bond_style": " ".join(self.bond_style),
            "improper_style": self.improper_style,
            "tail_correction": self.tail_correction,
            "kspace": self.kspace,
            "energy_fixes": list(self.energy_fixes),
        }


@dataclass(frozen=True)
class CharmmfswConventions:
    """Data-file-dependent CHARMM convention energies, kJ/mol."""

    one_four_lj_offset_kj_mol: float
    switch_shift_kj_mol: float


def _logical_lines(text: str) -> list[list[str]]:
    lines: list[list[str]] = []
    pending: list[str] = []
    for raw in text.splitlines():
        tokens = raw.split("#", 1)[0].split()
        continued = bool(tokens) and tokens[-1] == "&"
        if continued:
            tokens = tokens[:-1]
        pending.extend(tokens)
        if not continued and pending:
            lines.append(pending)
            pending = []
    if pending:
        lines.append(pending)
    return lines


def _keyword_value(tokens: list[str], keyword: str) -> str | None:
    for index, token in enumerate(tokens[:-1]):
        if token == keyword:
            return tokens[index + 1]
    return None


def parse_lammps_input(text: str) -> LammpsSettings:  # noqa: PLR0912
    """Parse the energy-relevant commands of a single-point LAMMPS input."""
    units = "lj"  # LAMMPS default when no units command is present
    pair: list[str] = []
    dihedral_style: str | None = None
    special_bonds: str | None = None
    bond_style: tuple[str, ...] = ()
    improper_style: str | None = None
    tail = False
    kspace = False
    normalized = False
    thermo: list[str] | None = None
    energy_fixes: list[str] = []
    for tokens in _logical_lines(text):
        command, arguments = tokens[0], tokens[1:]
        if command == "units" and arguments:
            units = arguments[0]
        elif command == "pair_style":
            pair = arguments
        elif command == "pair_modify":
            value = _keyword_value(arguments, "tail")
            if value is not None:
                tail = value == "yes"
        elif command == "dihedral_style" and arguments:
            dihedral_style = arguments[0]
        elif command == "special_bonds":
            special_bonds = " ".join(arguments)
        elif command == "bond_style":
            bond_style = tuple(arguments)
        elif command == "improper_style" and arguments:
            improper_style = arguments[0]
        elif command == "kspace_style" and arguments:
            kspace = arguments[0] != "none"
        elif command == "thermo_modify":
            normalized = _keyword_value(arguments, "norm") == "yes" or normalized
        elif command == "thermo_style":
            if not arguments or arguments[0] != "custom":
                raise ValueError("LAMMPS input must use thermo_style custom")
            thermo = arguments[1:]
        elif command == "fix_modify" and len(arguments) > 1:
            if _keyword_value(arguments[1:], "energy") == "yes":
                energy_fixes.append(arguments[0])
    if not pair:
        raise ValueError("LAMMPS input has no pair_style")
    if thermo is None:
        raise ValueError("LAMMPS input has no thermo_style")
    return LammpsSettings(
        units=units,
        pair_style=pair[0],
        pair_arguments=tuple(pair[1:]),
        dihedral_style=dihedral_style,
        special_bonds=special_bonds,
        bond_style=bond_style,
        improper_style=improper_style,
        tail_correction=tail,
        kspace=kspace,
        thermo_normalized=normalized,
        thermo_labels=tuple(_THERMO_LABELS.get(keyword, keyword) for keyword in thermo),
        energy_fixes=tuple(energy_fixes),
    )


def _numeric_row(line: str) -> list[float] | None:
    try:
        return [float(value) for value in line.split()]
    except ValueError:
        return None


def parse_lammps_thermo(source: str | Path) -> dict[str, float]:
    """Parse the single thermo row of a ``run 0`` LAMMPS log, in native units.

    A single-point adapter must not choose one frame from a trajectory, so a
    log with more than one thermo block or row is rejected.
    """
    text = source.read_text(encoding="utf-8") if isinstance(source, Path) else source
    lines = text.splitlines()
    headers = [index for index, line in enumerate(lines) if line.split()[:1] == ["Step"]]
    if len(headers) != 1:
        raise ValueError(f"LAMMPS log must contain one thermo block, found {len(headers)}")
    names = lines[headers[0]].split()
    if len(set(names)) != len(names):
        raise ValueError("LAMMPS thermo header contains duplicate columns")
    rows: list[list[float]] = []
    for line in lines[headers[0] + 1 :]:
        row = _numeric_row(line)
        if row is None or not row:
            break
        rows.append(row)
    if len(rows) != 1:
        raise ValueError(f"LAMMPS single-point log must contain one thermo row, found {len(rows)}")
    if len(rows[0]) != len(names):
        raise ValueError("LAMMPS thermo row does not match its header")
    if not all(math.isfinite(value) for value in rows[0]):
        raise ValueError("LAMMPS thermo row contains a non-finite value")
    return dict(zip(names, rows[0], strict=True))


def _validate_settings(settings: LammpsSettings) -> None:
    if settings.units != "real":
        raise ValueError(f"LAMMPS units {settings.units!r} are unsupported; expected real")
    if settings.thermo_normalized:
        raise ValueError("LAMMPS thermo energies are normalized per atom (thermo_modify norm yes)")
    if settings.pair_style not in _SUPPORTED_PAIR_STYLES:
        raise ValueError(f"unsupported LAMMPS pair_style {settings.pair_style!r}")
    if settings.charmm:
        if settings.dihedral_style != "charmmfsw" or settings.special_bonds != "charmm":
            raise ValueError(
                "lj/charmmfsw/coul/long requires dihedral_style charmmfsw and special_bonds charmm"
            )
        if len(settings.pair_arguments) < _CHARMM_SWITCH_ARGUMENTS:
            raise ValueError("lj/charmmfsw/coul/long requires inner and outer switch distances")
        if settings.explicit_pair_table:
            raise ValueError("tabulated explicit pairs are unsupported with CHARMM pair styles")
    elif settings.dihedral_style == "charmmfsw":
        raise ValueError("dihedral_style charmmfsw requires pair_style lj/charmmfsw/coul/long")


def charmm_switch_distances(settings: LammpsSettings) -> tuple[float, float]:
    """Inner and outer CHARMM force-switch distances, angstrom."""
    inner, outer = (float(value) for value in settings.pair_arguments[:_CHARMM_SWITCH_ARGUMENTS])
    if not 0.0 < inner < outer:
        raise ValueError("CHARMM switch distances must satisfy 0 < inner < outer")
    return inner, outer


def _section(data: str, title: str) -> list[str]:
    marker = f"\n{title}\n\n"
    if marker not in data:
        raise ValueError(f"LAMMPS data file has no {title} section")
    return data.split(marker, 1)[1].split("\n\n", 1)[0].splitlines()


def charmmfsw_one_four_offset_kj_mol(data: str, inner: float, outer: float) -> float:
    """Return LAMMPS charmmfsw's constant (force-free) 1-4 LJ energy offset, kJ/mol.

    For each weighted 1-4 pair, dihedral_style charmmfsw adds
    -4 eps sigma^12 / (r_in^6 r_c^6) + 4 eps sigma^6 / (r_in^3 r_c^3).
    This is the convention documented in validate_charmm36_engines.py.
    """
    atom_type = {int(row.split()[0]): int(row.split()[2]) for row in _section(data, "Atoms # full")}
    one_four = {
        (int(fields[0]), int(fields[1])): (float(fields[4]), float(fields[5]))
        for fields in (row.split() for row in _section(data, "PairIJ Coeffs"))
    }
    weight = {
        int(row.split()[0]): float(row.split()[4]) for row in _section(data, "Dihedral Coeffs")
    }
    total = 0.0
    for row in _section(data, "Dihedrals"):
        fields = [int(value) for value in row.split()]
        if weight[fields[1]] != 1.0:
            continue
        first, second = sorted((atom_type[fields[2]], atom_type[fields[5]]))
        epsilon, sigma = one_four[(first, second)]
        total += -4 * epsilon * sigma**12 / (inner**6 * outer**6) + 4 * epsilon * sigma**6 / (
            inner**3 * outer**3
        )
    return total * KJ_PER_KCAL


def _gromacs_switched(
    r: npt.NDArray[np.float64], power: int, inner: float, outer: float
) -> npt.NDArray[np.float64]:
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


def _charmm_switched(
    r: npt.NDArray[np.float64], power: int, inner: float, outer: float
) -> npt.NDArray[np.float64]:
    half = power // 2
    k = outer**half / (outer**half - inner**half)
    return np.asarray(
        np.where(
            r <= inner,
            r**-power - 1 / (inner**half * outer**half),
            k * (r**-half - outer**-half) ** 2,
        )
    )


def charmmfsw_switch_shift_kj_mol(data: str, inner: float, outer: float) -> float:
    """Return CHARMM vfswitch minus GROMACS force-switch short-range LJ energy, kJ/mol.

    Both switch the force identically between ``inner`` and ``outer`` but shift
    the potential by different constants. Each convention is evaluated over
    every non-excluded pair (exclusions to three bonds) in the LAMMPS data
    file, as documented in validate_charmm36_engines.py.
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
    rows = [row.split() for row in _section(data, "Atoms # full")]
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
    for row in _section(data, "PairIJ Coeffs"):
        first, second, e, s, *_ = row.split()
        for a, b in ((int(first), int(second)), (int(second), int(first))):
            epsilon[a, b] = float(e)
            sigma[a, b] = float(s)
    neighbours: list[set[int]] = [set() for _ in range(count)]
    for row in _section(data, "Bonds"):
        _, _, bonded_first, bonded_second = (int(v) for v in row.split())
        neighbours[bonded_first - 1].add(bonded_second - 1)
        neighbours[bonded_second - 1].add(bonded_first - 1)

    inverse = np.linalg.inv(box)
    total = 0.0
    for index in range(count - 1):
        excluded = {index}
        frontier = {index}
        for _ in range(_EXCLUDED_BOND_DEPTH):
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
                repulsion
                * (_charmm_switched(r, 12, inner, outer) - _gromacs_switched(r, 12, inner, outer))
                - dispersion
                * (_charmm_switched(r, 6, inner, outer) - _gromacs_switched(r, 6, inner, outer))
            )
        )
    return total * KJ_PER_KCAL


def charmmfsw_conventions(settings: LammpsSettings, data: str) -> CharmmfswConventions:
    """Evaluate the two data-file-dependent CHARMM conventions."""
    inner, outer = charmm_switch_distances(settings)
    return CharmmfswConventions(
        one_four_lj_offset_kj_mol=charmmfsw_one_four_offset_kj_mol(data, inner, outer),
        switch_shift_kj_mol=charmmfsw_switch_shift_kj_mol(data, inner, outer),
    )


def _finite(value: object, description: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise TypeError(f"{description} is not numeric")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{description} is not finite")
    return number


def _classify_columns(
    thermo: Mapping[str, float], settings: LammpsSettings
) -> tuple[dict[str, float], dict[str, float]]:
    """Return (component columns, restraint columns) in kJ/mol; fail closed."""
    if set(thermo) != set(settings.thermo_labels):
        raise ValueError("LAMMPS thermo columns do not match the input thermo_style")
    required = list(_ALWAYS_REQUIRED)
    if settings.kspace:
        required.append("E_long")
    if settings.tail_correction:
        required.append("E_tail")
    if settings.improper_style is not None:
        required.append("E_impro")
    required.extend(f"f_{fix}" for fix in settings.energy_fixes)
    missing = [label for label in required if label not in thermo]
    if missing:
        raise ValueError(f"LAMMPS thermo lacks energy columns the input requires: {missing}")
    restraint_labels = {f"f_{fix}" for fix in settings.energy_fixes}
    components: dict[str, float] = {}
    restraints: dict[str, float] = {}
    unsupported: list[str] = []
    for label, raw in thermo.items():
        value = _finite(raw, f"LAMMPS column {label!r}")
        if label in _NON_POTENTIAL_LABELS or label in _SUBTOTAL_LABELS:
            continue
        if label in _COMPONENT_LABELS:
            components[label] = value * KJ_PER_KCAL
        elif label in restraint_labels:
            restraints[label] = value * KJ_PER_KCAL
        elif value != 0:
            unsupported.append(label)
    if unsupported:
        raise ValueError(
            f"LAMMPS output contains unsupported nonzero energy terms: {sorted(unsupported)}"
        )
    return components, restraints


def _decompose(
    grouped: dict[str, float],
    settings: LammpsSettings,
    *,
    periodic_improper_count: int,
    periodic_improper_kj_mol: float | None,
    explicit_pair_kj_mol: Mapping[str, float] | None,
) -> list[dict[str, Any]]:
    decompositions: list[dict[str, Any]] = []
    if grouped["dispersion_correction"] != 0:
        decompositions.append(
            {
                "name": "lammps_evdwl_includes_etail",
                "from": "vdw",
                "to": "dispersion_correction",
                "energy_kj_mol": grouped["dispersion_correction"],
                "reason": "LAMMPS evdwl includes the etail long-range tail correction",
            }
        )
    if isinstance(periodic_improper_count, bool) or not isinstance(periodic_improper_count, int):
        raise TypeError("periodic_improper_count must be an integer")
    if periodic_improper_count < 0:
        raise ValueError("periodic_improper_count must be non-negative")
    if periodic_improper_count == 0:
        if periodic_improper_kj_mol is not None:
            raise ValueError("periodic_improper_kj_mol was given but the topology has none")
    else:
        if periodic_improper_kj_mol is None:
            raise ValueError(
                "the topology has periodic impropers, which LAMMPS reports inside edihed; "
                "supply their separately evaluated energy"
            )
        energy = _finite(periodic_improper_kj_mol, "periodic_improper_kj_mol")
        grouped["dihedral"] -= energy
        grouped["improper"] += energy
        decompositions.append(
            {
                "name": "lammps_periodic_impropers_in_edihed",
                "from": "dihedral",
                "to": "improper",
                "energy_kj_mol": energy,
                "count": periodic_improper_count,
                "reason": "periodic impropers are written as LAMMPS dihedrals; classified "
                "from the CelluloseForge topology",
            }
        )
    if settings.explicit_pair_table:
        if explicit_pair_kj_mol is None:
            raise ValueError(
                "tabulated explicit 1-4 pairs are tallied in ebond; supply their separately "
                "evaluated vdw and electrostatics energies"
            )
        if set(explicit_pair_kj_mol) != {"vdw", "electrostatics"}:
            raise ValueError("explicit_pair_kj_mol must contain exactly vdw and electrostatics")
        vdw = _finite(explicit_pair_kj_mol["vdw"], "explicit pair vdw")
        electrostatics = _finite(
            explicit_pair_kj_mol["electrostatics"], "explicit pair electrostatics"
        )
        grouped["bond"] -= vdw + electrostatics
        grouped["vdw"] += vdw
        grouped["electrostatics"] += electrostatics
        decompositions.append(
            {
                "name": "lammps_explicit_pairs_in_ebond",
                "from": "bond",
                "to": {"vdw": vdw, "electrostatics": electrostatics},
                "energy_kj_mol": vdw + electrostatics,
                "reason": "explicit GLYCAM 1-4 pairs are tabulated pseudo-bonds tallied in ebond",
            }
        )
    elif explicit_pair_kj_mol is not None:
        raise ValueError("explicit_pair_kj_mol was given but the input has no pair table")
    return decompositions


def _charmm_corrections(
    grouped: dict[str, float], conventions: CharmmfswConventions
) -> list[dict[str, Any]]:
    """Apply the three documented CHARMM conventions; return visible records."""
    one_four = _finite(conventions.one_four_lj_offset_kj_mol, "1-4 LJ offset")
    switch = _finite(conventions.switch_shift_kj_mol, "switch shift")
    corrections: list[dict[str, Any]] = []

    raw = grouped["electrostatics"]
    comparable = raw * CHARMM_COULOMB_RATIO
    corrections.append(
        {
            "name": "lammps_charmm_coulomb_constant",
            "kind": "scale",
            "term": "electrostatics",
            "factor": CHARMM_COULOMB_RATIO,
            "raw_kj_mol": raw,
            "correction_kj_mol": comparable - raw,
            "comparable_kj_mol": comparable,
            "reason": f"LAMMPS CHARMM Coulomb constant {LAMMPS_CHARMM_COULOMB} vs GROMACS "
            f"{GROMACS_COULOMB:.4f} kcal A/mol/e^2; applied to ecoul and elong",
            "reference": "validation/validate_charmm36_engines.py",
        }
    )
    grouped["electrostatics"] = comparable

    for name, delta, reason in (
        (
            "charmmfsw_one_four_lj_offset",
            -one_four,
            "constant per-pair energy offset of dihedral_style charmmfsw, computed "
            "analytically and subtracted",
        ),
        (
            "charmmfsw_force_switch_potential_shift",
            -switch,
            "CHARMM vfswitch (lj/charmmfsw) and GROMACS force-switch differ by a constant "
            "potential shift, evaluated over every non-excluded pair and subtracted",
        ),
    ):
        raw = grouped["vdw"]
        grouped["vdw"] = raw + delta
        corrections.append(
            {
                "name": name,
                "kind": "additive",
                "term": "vdw",
                "raw_kj_mol": raw,
                "correction_kj_mol": delta,
                "comparable_kj_mol": grouped["vdw"],
                "reason": reason,
                "reference": "validation/validate_charmm36_engines.py",
            }
        )
    return corrections


def normalize_lammps_energies(  # noqa: PLR0913
    thermo: Mapping[str, float],
    *,
    settings: LammpsSettings,
    coordinate_sha256: str,
    periodic_improper_count: int,
    periodic_improper_kj_mol: float | None = None,
    explicit_pair_kj_mol: Mapping[str, float] | None = None,
    charmm_conventions: CharmmfswConventions | None = None,
) -> EnergyReport:
    """Return a complete canonical decomposition of a LAMMPS single point.

    ``thermo`` is one native-unit thermo row keyed by header label.  Raw
    thermo values, the uncorrected canonical terms, decompositions and
    corrections are all retained in the report.
    """
    _validate_settings(settings)
    if settings.charmm and charmm_conventions is None:
        raise ValueError("lj/charmmfsw/coul/long requires the CHARMM convention energies")
    if not settings.charmm and charmm_conventions is not None:
        raise ValueError("CHARMM conventions apply only to lj/charmmfsw/coul/long")
    components, restraints = _classify_columns(thermo, settings)

    def column(label: str) -> float:
        return components.get(label, 0.0)

    tail = column("E_tail")
    grouped: dict[str, float] = {
        "bond": column("E_bond"),
        "angle": column("E_angle"),  # angle_style charmm includes Urey-Bradley
        "dihedral": column("E_dihed"),
        "improper": column("E_impro"),
        "vdw": column("E_vdwl") - tail,
        "electrostatics": column("E_coul") + column("E_long"),
        "restraint": sum(restraints.values(), 0.0),
        "dispersion_correction": tail,
    }
    decompositions = _decompose(
        grouped,
        settings,
        periodic_improper_count=periodic_improper_count,
        periodic_improper_kj_mol=periodic_improper_kj_mol,
        explicit_pair_kj_mol=explicit_pair_kj_mol,
    )
    uncorrected = dict(grouped)
    raw_potential = column("PotEng")
    corrections: list[dict[str, Any]] = []
    if charmm_conventions is not None:
        corrections = _charmm_corrections(grouped, charmm_conventions)
    total_correction = sum(item["correction_kj_mol"] for item in corrections)
    terms: dict[str, EnergyValue] = {**grouped, "potential": raw_potential + total_correction}
    return EnergyReport(
        engine="lammps",
        coordinate_sha256=coordinate_sha256,
        terms=terms,
        corrections=tuple(corrections),
        extra={
            "raw_thermo_kcal_mol": {
                label: float(value)
                for label, value in thermo.items()
                if label not in _NON_POTENTIAL_LABELS
            },
            "raw_terms_kj_mol": {**uncorrected, "potential": raw_potential},
            "decompositions": decompositions,
            "lammps_settings": settings.as_dict(),
            "unit_conversion": f"LAMMPS real kcal/mol x {KJ_PER_KCAL}",
        },
    )


def lammps_energy_report(  # noqa: PLR0913
    log: str | Path,
    *,
    input_text: str,
    coordinate_sha256: str,
    periodic_improper_count: int,
    periodic_improper_kj_mol: float | None = None,
    explicit_pair_kj_mol: Mapping[str, float] | None = None,
    data_text: str | None = None,
) -> EnergyReport:
    """Parse a ``run 0`` log with its input (and data file for CHARMM) into a report."""
    settings = parse_lammps_input(input_text)
    _validate_settings(settings)
    conventions: CharmmfswConventions | None = None
    if settings.charmm:
        if data_text is None:
            raise ValueError("CHARMM conventions require the LAMMPS data file")
        conventions = charmmfsw_conventions(settings, data_text)
    return normalize_lammps_energies(
        parse_lammps_thermo(log),
        settings=settings,
        coordinate_sha256=coordinate_sha256,
        periodic_improper_count=periodic_improper_count,
        periodic_improper_kj_mol=periodic_improper_kj_mol,
        explicit_pair_kj_mol=explicit_pair_kj_mol,
        charmm_conventions=conventions,
    )
