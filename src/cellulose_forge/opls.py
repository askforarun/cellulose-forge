"""Experimental OPLS-AA and OPLS-AA-SEI topologies for crystalline cellulose.

Atom types, charges, Lennard-Jones parameters, bonds, and angles are the
standard OPLS-AA carbohydrate parameters (Damm, Frontera, Tirado-Rives and
Jorgensen, J. Comput. Chem. 1997, 18, 1955) as distributed in the GROMACS
``oplsaa.ff``, which every GROMACS installation provides.  Torsions are
written explicitly from Table 3 of Kony, Damm, Stoll and van Gunsteren,
J. Comput. Chem. 2002, 23, 1416 (doi:10.1002/jcc.10139): its OPLS-AA column
for ``opls-aa`` and its OPLS-AA-SEI column for ``opls-aa-sei``.  OPLS-AA-SEI
additionally divides 1-5 electrostatic interactions within each
hexopyranose by 1.26 and the 1-6 interactions between the hydroxymethyl and
the vicinal C4 hydroxyl by 1.22, with Lennard-Jones kept at full strength.

Interpretations required to apply the monosaccharide parameter set to
cellulose (all recorded in the build report):

* The 1-5 and 1-6 scaling is applied within each residue only, as published.
  The 1-6 term is applied only where O4 is a hydroxyl (the non-reducing end);
  in internal residues O4 is a glycosidic ether, not a hydroxyl.
* OS-C-C-OS (O5-C5-C4-O4 across a glycosidic ether) does not occur in a free
  hexopyranose and is absent from Table 3; the OPLS-AA polyether value
  (V1 = -0.55 kcal/mol, GROMACS ``oplsaa.ff``) is used in both variants.
* C-C-C-OS_ace is assigned to C3-C2-C1-O where O is the exocyclic acetal
  oxygen on C1 (the glycosidic oxygen); C3-C2-C1-O5 uses C-C-C-OS.
"""

from __future__ import annotations

import json
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from cellulose_forge.geometry import BulkCrystal
from cellulose_forge.glycam_mapping import GlycamMappedSystem
from cellulose_forge.gromacs import write_gro

KCAL_TO_KJ = 4.184
VARIANTS = ("opls-aa", "opls-aa-sei")
SEI_ONE_FIVE_SCALE = 1.26
SEI_ONE_SIX_SCALE = 1.22
_NEUTRALITY_TOLERANCE_E = 1.0e-9
VALIDATION_SCOPE = (
    "topology construction validated against the published OPLS-AA/OPLS-AA-SEI "
    "parameters; force-field accuracy for crystalline cellulose not assessed"
)

# (name, element) in the same per-residue order as the OpenFF builder.
_RESIDUE_ATOMS = (
    ("C1", "C"), ("H1", "H"), ("O5", "O"), ("C5", "C"), ("H5", "H"), ("C6", "C"),
    ("H62", "H"), ("H61", "H"), ("O6", "O"), ("HO6", "H"), ("C4", "C"), ("H4", "H"),
    ("C3", "C"), ("H3", "H"), ("O3", "O"), ("HO3", "H"), ("C2", "C"), ("H2", "H"),
    ("O2", "O"), ("HO2", "H"), ("O4", "O"),
)
_RESIDUE_BONDS = (
    ("C1", "H1"), ("C1", "O5"), ("C1", "C2"), ("O5", "C5"), ("C5", "H5"), ("C5", "C6"),
    ("C5", "C4"), ("C6", "H62"), ("C6", "H61"), ("C6", "O6"), ("O6", "HO6"),
    ("C4", "H4"), ("C4", "C3"), ("C4", "O4"), ("C3", "H3"), ("C3", "O3"), ("C3", "C2"),
    ("O3", "HO3"), ("C2", "H2"), ("C2", "O2"), ("O2", "HO2"),
)

# OPLS-AA atom types (GROMACS oplsaa.ff ffnonbonded.itp); charge in e.
_OPLS_TYPES: dict[str, tuple[str, float, float]] = {
    # OPLS type -> (bond_type, charge, mass)
    "opls_140": ("HC", 0.060, 1.008),    # H on hydroxyl-bearing / CH2OH carbon
    "opls_157": ("CT", 0.145, 12.011),   # C6 (CH2OH)
    "opls_158": ("CT", 0.205, 12.011),   # CH bearing OH (C2, C3, terminal C4)
    "opls_169": ("OH", -0.700, 15.9994), # hydroxyl O (diol)
    "opls_170": ("HO", 0.435, 1.008),    # hydroxyl H (diol)
    "opls_183": ("CT", 0.170, 12.011),   # CH bonded to ether O (C5, internal C4)
    "opls_185": ("HC", 0.030, 1.008),    # H on ether CH
    "opls_186": ("OS", -0.400, 15.9994), # acetal ether O (O5, glycosidic O4)
    "opls_187": ("OH", -0.700, 15.9994), # hemiacetal hydroxyl O (O1)
    "opls_188": ("HO", 0.435, 1.008),    # hemiacetal hydroxyl H (HO1)
    "opls_193": ("CO", 0.300, 12.011),   # acetal CH (internal C1)
    "opls_194": ("HC", 0.100, 1.008),    # H on acetal CH
    "opls_195": ("CO", 0.365, 12.011),   # hemiacetal CH (reducing-end C1)
    "opls_196": ("HC", 0.100, 1.008),    # H on hemiacetal CH
}

# Kony et al. 2002 Table 3 (and the OS-C-C-OS polyether term), kcal/mol
# Fourier coefficients V1, V2, V3 for E = V1/2(1+cos) + V2/2(1-cos2) + V3/2(1+cos3).
_TORSIONS: dict[str, dict[str, tuple[float, float, float]]] = {
    "opls-aa": {
        "C-C-C-C": (1.740, -0.157, 0.279),
        "C-C-C-H": (0.000, 0.000, 0.366),
        "H-C-O-C": (0.000, 0.000, 0.760),
        "C-O-C-C": (0.650, -0.250, 0.670),
        "H-C-C-O": (0.000, 0.000, 0.468),
        "H-C-C-H": (0.000, 0.000, 0.318),
        "H-C-O-H": (0.000, 0.000, 0.450),
        "O-C-O-C": (-0.375, -1.358, 0.004),
        "OS-C-C-OH": (4.319, 0.000, 0.000),
        "C-C-O-H": (2.674, -2.883, 1.026),
        "OH-C-C-OH": (9.066, 0.000, 0.000),
        "O-C-O-H": (-1.257, -1.806, 0.003),
        "C-C-C-OH": (-1.336, 0.000, 0.000),
        "C-C-C-OS": (-1.336, 0.000, 0.000),
        "C-C-C-OS_ace": (-1.336, 0.000, 0.000),
        "OS-C-C-OS": (-0.550, 0.000, 0.000),
    },
}
_TORSIONS["opls-aa-sei"] = {
    **_TORSIONS["opls-aa"],
    "O-C-O-C": (-3.687, 0.005, 0.537),
    "OS-C-C-OH": (1.593, 0.000, 0.000),
    "C-C-O-H": (0.205, -1.459, 0.010),
    "OH-C-C-OH": (8.519, -6.829, 4.793),
    "O-C-O-H": (-1.673, -2.200, 0.001),
    "C-C-C-OH": (-2.587, 0.000, 0.000),
    "C-C-C-OS": (-0.089, -0.729, 0.870),
    "C-C-C-OS_ace": (1.501, 0.000, 0.000),
}

# Lennard-Jones parameters (sigma nm, epsilon kJ/mol) from GROMACS oplsaa.ff.
_LJ = {
    "opls_140": (0.250, 0.125520), "opls_157": (0.350, 0.276144),
    "opls_158": (0.350, 0.276144), "opls_169": (0.307, 0.711280),
    "opls_170": (0.000, 0.000000), "opls_183": (0.350, 0.276144),
    "opls_185": (0.250, 0.125520), "opls_186": (0.290, 0.585760),
    "opls_187": (0.307, 0.711280), "opls_188": (0.000, 0.000000),
    "opls_193": (0.350, 0.276144), "opls_194": (0.250, 0.125520),
    "opls_195": (0.350, 0.276144), "opls_196": (0.250, 0.125520),
}

# Harmonic bonds (b0 nm, kb kJ mol-1 nm-2; E = kb/2 (b - b0)^2) and angles
# (theta0 degrees, k kJ mol-1 rad-2; E = k/2 (theta - theta0)^2) for the bond
# types above, copied from GROMACS oplsaa.ff ffbonded.itp.  They are written
# explicitly so that the GROMACS, OpenMM and LAMMPS outputs share one source;
# validate_opls_cellulose.py checks them against the installed oplsaa.ff.
_BOND_PARAMETERS: dict[tuple[str, str], tuple[float, float]] = {
    ("CO", "CT"): (0.1529, 224262.4), ("CO", "HC"): (0.109, 284512.0),
    ("CO", "OH"): (0.138, 267776.0), ("CO", "OS"): (0.138, 267776.0),
    ("CT", "CT"): (0.1529, 224262.4), ("CT", "HC"): (0.109, 284512.0),
    ("CT", "OH"): (0.141, 267776.0), ("CT", "OS"): (0.141, 267776.0),
    ("HO", "OH"): (0.0945, 462750.4),
}
# Keyed by (outer, centre, outer) with the outer types in sorted order.
_ANGLE_PARAMETERS: dict[tuple[str, str, str], tuple[float, float]] = {
    ("CT", "CT", "CO"): (112.7, 488.273), ("CO", "CT", "HC"): (110.7, 313.8),
    ("CO", "CT", "OH"): (109.5, 418.4), ("CO", "OH", "HO"): (108.5, 460.24),
    ("CO", "OS", "CT"): (109.5, 502.08), ("CT", "CO", "HC"): (110.7, 313.8),
    ("CT", "CO", "OH"): (109.5, 418.4), ("CT", "CO", "OS"): (109.5, 418.4),
    ("CT", "CT", "CT"): (112.7, 488.273), ("CT", "CT", "HC"): (110.7, 313.8),
    ("CT", "CT", "OH"): (109.5, 418.4), ("CT", "CT", "OS"): (109.5, 418.4),
    ("CT", "OH", "HO"): (108.5, 460.24), ("HC", "CO", "OH"): (109.5, 292.88),
    ("HC", "CO", "OS"): (109.5, 292.88), ("HC", "CT", "HC"): (107.8, 276.144),
    ("HC", "CT", "OH"): (109.5, 292.88), ("HC", "CT", "OS"): (109.5, 292.88),
    ("OH", "CO", "OS"): (111.55, 774.876), ("OS", "CO", "OS"): (111.55, 774.876),
}
# Fudge factors of the OPLS-AA [ defaults ] directive.
OPLS_FUDGE_QQ = 0.5
OPLS_FUDGE_LJ = 0.5


@dataclass(frozen=True)
class _Atom:
    name: str
    element: str
    residue: int
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class _Chain:
    """Lightweight chain graph accepted by the shared coordinate mapper."""

    atoms: tuple[_Atom, ...]
    bonds: tuple[tuple[int, int], ...]

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)


@dataclass(frozen=True)
class OplsTopology:
    """One cellulose chain parameterized with OPLS-AA or OPLS-AA-SEI."""

    variant: str
    degree_polymerization: int
    chain_boundary: str
    atoms: tuple[_Atom, ...]
    atom_types: tuple[str, ...]
    charges: tuple[float, ...]
    bonds: tuple[tuple[int, int], ...]
    angles: tuple[tuple[int, int, int], ...]
    dihedrals: tuple[tuple[int, int, int, int, str], ...]
    pairs_14: tuple[tuple[int, int], ...]
    scaled_pairs: tuple[tuple[int, int, float, str], ...] = ()

    def class_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for *_, name in self.dihedrals:
            counts[name] = counts.get(name, 0) + 1
        return dict(sorted(counts.items()))


def build_chain_graph(degree_polymerization: int, chain_boundary: str) -> _Chain:
    """Return the finite or covalently periodic beta-(1->4) glucan graph."""
    if degree_polymerization < 2:  # noqa: PLR2004
        raise ValueError("OPLS cellulose chains require DP >= 2")
    if chain_boundary not in {"finite", "periodic"}:
        raise ValueError("OPLS chain boundary must be finite or periodic")
    atoms: list[_Atom] = []
    index: dict[tuple[int, str], int] = {}

    def add(name: str, element: str, residue: int) -> None:
        index[residue, name] = len(atoms)
        atoms.append(_Atom(name, element, residue, {"residue_number": residue}))

    for residue in range(1, degree_polymerization + 1):
        for name, element in _RESIDUE_ATOMS:
            add(name, element, residue)
        if chain_boundary == "finite" and residue == 1:
            add("O1", "O", residue)
            add("HO1", "H", residue)
        if chain_boundary == "finite" and residue == degree_polymerization:
            add("HO4", "H", residue)
    bonds = [
        (index[residue, left], index[residue, right])
        for residue in range(1, degree_polymerization + 1)
        for left, right in _RESIDUE_BONDS
    ]
    bonds += [
        (index[residue, "O4"], index[residue + 1, "C1"])
        for residue in range(1, degree_polymerization)
    ]
    if chain_boundary == "periodic":
        bonds.append((index[degree_polymerization, "O4"], index[1, "C1"]))
    else:
        bonds += [
            (index[1, "C1"], index[1, "O1"]),
            (index[1, "O1"], index[1, "HO1"]),
            (index[degree_polymerization, "O4"], index[degree_polymerization, "HO4"]),
        ]
    return _Chain(tuple(atoms), tuple(bonds))


def _neighbours(chain: _Chain) -> list[list[int]]:
    neighbours: list[list[int]] = [[] for _ in chain.atoms]
    for left, right in chain.bonds:
        neighbours[left].append(right)
        neighbours[right].append(left)
    return neighbours


_FIXED_TYPES = {
    "O5": "opls_186", "C5": "opls_183", "H5": "opls_185", "C6": "opls_157",
    "H61": "opls_140", "H62": "opls_140", "H2": "opls_140", "H3": "opls_140",
    "C2": "opls_158", "C3": "opls_158", "O2": "opls_169", "O3": "opls_169",
    "O6": "opls_169", "HO2": "opls_170", "HO3": "opls_170", "HO6": "opls_170",
    "HO4": "opls_170", "O1": "opls_187", "HO1": "opls_188",
}
# (type in a reducing-end or hydroxyl-terminated residue, type otherwise)
_ROLE_TYPES = {
    "C1": ("opls_195", "opls_193"), "H1": ("opls_196", "opls_194"),
    "C4": ("opls_158", "opls_183"), "H4": ("opls_140", "opls_185"),
    "O4": ("opls_169", "opls_186"),
}


def _assign_types(chain: _Chain) -> tuple[str, ...]:
    """Assign Damm 1997 OPLS-AA carbohydrate atom types."""
    types: list[str] = []
    for atom in chain.atoms:
        if atom.name in _FIXED_TYPES:
            types.append(_FIXED_TYPES[atom.name])
            continue
        if atom.name not in _ROLE_TYPES:
            raise ValueError(f"no OPLS-AA type for atom {atom.name}")
        # C1/H1 change at the reducing end (hemiacetal); C4/H4/O4 change at the
        # non-reducing end, where O4 is a hydroxyl instead of a glycosidic ether.
        marker = "O1" if atom.name in {"C1", "H1"} else "HO4"
        terminal, internal = _ROLE_TYPES[atom.name]
        types.append(terminal if _residue_has(chain, atom.residue, marker) else internal)
    return tuple(types)


def _residue_has(chain: _Chain, residue: int, name: str) -> bool:
    return any(atom.residue == residue and atom.name == name for atom in chain.atoms)


def _distances(neighbours: list[list[int]], source: int, limit: int) -> dict[int, int]:
    seen = {source: 0}
    queue = deque([source])
    while queue:
        node = queue.popleft()
        if seen[node] == limit:
            continue
        for other in neighbours[node]:
            if other not in seen:
                seen[other] = seen[node] + 1
                queue.append(other)
    return seen


def _oxygen_class(atom_index: int, neighbours: list[list[int]], chain: _Chain) -> str:
    """Return OH for a hydroxyl oxygen and OS for an ether or acetal oxygen."""
    return (
        "OH"
        if any(chain.atoms[other].element == "H" for other in neighbours[atom_index])
        else "OS"
    )


def _classify(
    quad: tuple[int, int, int, int], neighbours: list[list[int]], chain: _Chain
) -> str:
    """Return the Kony Table 3 torsion class of one dihedral."""
    atoms = [chain.atoms[index] for index in quad]
    detailed = [
        _oxygen_class(index, neighbours, chain) if atom.element == "O" else atom.element
        for index, atom in zip(quad, atoms, strict=True)
    ]
    generic = ["O" if label in {"OH", "OS"} else label for label in detailed]

    def matches(pattern: str, labels: list[str]) -> bool:
        parts = pattern.split("-")
        return parts == labels or parts == labels[::-1]

    # O-C-C-O classes distinguish hydroxyl and ether oxygens.
    if generic == ["O", "C", "C", "O"]:
        ends = sorted((detailed[0], detailed[3]))
        return {
            ("OH", "OH"): "OH-C-C-OH",
            ("OH", "OS"): "OS-C-C-OH",
            ("OS", "OS"): "OS-C-C-OS",
        }[(ends[0], ends[1])]
    # C-C-C-O classes: hydroxyl, ether, or the exocyclic acetal oxygen on C1.
    for labels, order in ((generic, quad), (generic[::-1], quad[::-1])):
        if labels == ["C", "C", "C", "O"]:
            oxygen, carbon = order[3], order[2]
            if _oxygen_class(oxygen, neighbours, chain) == "OH":
                return "C-C-C-OH"
            anomeric = chain.atoms[carbon].name == "C1"
            exocyclic = chain.atoms[oxygen].name != "O5"
            return "C-C-C-OS_ace" if anomeric and exocyclic else "C-C-C-OS"
    for pattern in (
        "O-C-O-C", "O-C-O-H", "C-C-O-H", "H-C-O-H", "H-C-O-C", "C-O-C-C",
        "C-C-C-C", "C-C-C-H", "H-C-C-O", "H-C-C-H",
    ):
        if matches(pattern, generic):
            return pattern
    raise ValueError(f"no OPLS-AA torsion class for {'-'.join(detailed)}")


def build_opls_topology(  # noqa: PLR0912 - explicit topology enumeration
    variant: str, degree_polymerization: int, chain_boundary: str
) -> OplsTopology:
    """Construct one OPLS-AA or OPLS-AA-SEI cellulose chain topology."""
    if variant not in VARIANTS:
        raise ValueError(f"unknown OPLS variant {variant!r}")
    chain = build_chain_graph(degree_polymerization, chain_boundary)
    neighbours = _neighbours(chain)
    types = _assign_types(chain)
    charges = tuple(_OPLS_TYPES[atom_type][1] for atom_type in types)
    if abs(sum(charges)) > _NEUTRALITY_TOLERANCE_E:
        raise ValueError(f"OPLS cellulose chain is not neutral: {sum(charges):.6f} e")

    angles = sorted(
        {
            (min(a, c), b, max(a, c))
            for b, bonded in enumerate(neighbours)
            for a in bonded
            for c in bonded
            if a < c
        },
        key=lambda item: (item[1], item[0], item[2]),
    )
    dihedrals: set[tuple[int, int, int, int]] = set()
    for b, c in chain.bonds:
        for a in neighbours[b]:
            for d in neighbours[c]:
                if len({a, b, c, d}) == 4:  # noqa: PLR2004
                    quad = (a, b, c, d) if a < d else (d, c, b, a)
                    dihedrals.add(quad)
    classified = tuple(
        (*quad, _classify(quad, neighbours, chain)) for quad in sorted(dihedrals)
    )

    pairs_14: list[tuple[int, int]] = []
    scaled: list[tuple[int, int, float, str]] = []
    for source in range(len(chain.atoms)):
        distances = _distances(neighbours, source, 5)
        for target, distance in distances.items():
            if target <= source:
                continue
            if distance == 3:  # noqa: PLR2004
                pairs_14.append((source, target))
            same_residue = chain.atoms[source].residue == chain.atoms[target].residue
            if variant != "opls-aa-sei" or not same_residue:
                continue
            if distance == 4:  # noqa: PLR2004
                scaled.append((source, target, 1.0 / SEI_ONE_FIVE_SCALE, "1-5"))
            elif distance == 5 and _is_hydroxymethyl_vicinal_pair(  # noqa: PLR2004
                source, target, chain
            ):
                scaled.append((source, target, 1.0 / SEI_ONE_SIX_SCALE, "1-6"))
    return OplsTopology(
        variant=variant,
        degree_polymerization=degree_polymerization,
        chain_boundary=chain_boundary,
        atoms=chain.atoms,
        atom_types=types,
        charges=charges,
        bonds=tuple(sorted(tuple(sorted(bond)) for bond in chain.bonds)),  # type: ignore[misc]
        angles=tuple(angles),
        dihedrals=classified,
        pairs_14=tuple(pairs_14),
        scaled_pairs=tuple(scaled),
    )


def _is_hydroxymethyl_vicinal_pair(first: int, second: int, chain: _Chain) -> bool:
    """1-6 pairs between the C6 hydroxyl and a C4 hydroxyl (Kony Fig. 3)."""
    names = {chain.atoms[first].name, chain.atoms[second].name}
    return names in ({"HO6", "O4"}, {"O6", "HO4"}) and _residue_has(
        chain, chain.atoms[first].residue, "HO4"
    )


def _bond_type(topology: OplsTopology, index: int) -> str:
    return _OPLS_TYPES[topology.atom_types[index]][0]


def bond_parameters(topology: OplsTopology, first: int, second: int) -> tuple[float, float]:
    """Return (b0 nm, kb kJ mol-1 nm-2) for one bond."""
    key = (_bond_type(topology, first), _bond_type(topology, second))
    parameters = _BOND_PARAMETERS.get(key) or _BOND_PARAMETERS.get(key[::-1])
    if parameters is None:
        raise ValueError(f"no OPLS-AA bond parameters for {'-'.join(key)}")
    return parameters


def angle_parameters(
    topology: OplsTopology, first: int, centre: int, last: int
) -> tuple[float, float]:
    """Return (theta0 degrees, k kJ mol-1 rad-2) for one angle."""
    outer = sorted((_bond_type(topology, first), _bond_type(topology, last)))
    key = (outer[0], _bond_type(topology, centre), outer[1])
    for candidate in (key, key[::-1]):
        if candidate in _ANGLE_PARAMETERS:
            return _ANGLE_PARAMETERS[candidate]
    raise ValueError(f"no OPLS-AA angle parameters for {'-'.join(key)}")


def lj_pair(topology: OplsTopology, first: int, second: int) -> tuple[float, float]:
    """Geometric-mean (sigma nm, epsilon kJ/mol), OPLS combination rule 3."""
    left = _LJ[topology.atom_types[first]]
    right = _LJ[topology.atom_types[second]]
    return float(np.sqrt(left[0] * right[0])), float(np.sqrt(left[1] * right[1]))


def lj_parameters(atom_type: str) -> tuple[float, float]:
    """Return (sigma nm, epsilon kJ/mol) of one OPLS-AA atom type."""
    return _LJ[atom_type]


def atom_mass(atom_type: str) -> float:
    return _OPLS_TYPES[atom_type][2]


def nonbonded_exceptions(
    topology: OplsTopology,
) -> tuple[tuple[int, int, float, float, str], ...]:
    """Every intramolecular pair whose nonbonded term differs from the default.

    Each entry is (i, j, Coulomb scale, Lennard-Jones scale, kind): 1-2 and 1-3
    pairs are excluded (0, 0), 1-4 pairs use the OPLS fudge factors, and
    OPLS-AA-SEI 1-5/1-6 pairs scale Coulomb only.
    """
    excluded = {tuple(sorted(bond)) for bond in topology.bonds}
    excluded |= {(min(a, c), max(a, c)) for a, _, c in topology.angles}
    entries: list[tuple[int, int, float, float, str]] = [
        (a, b, 0.0, 0.0, "excluded") for a, b in sorted(excluded)
    ]
    entries += [(a, b, OPLS_FUDGE_QQ, OPLS_FUDGE_LJ, "1-4") for a, b in topology.pairs_14]
    entries += [(a, b, scale, 1.0, kind) for a, b, scale, kind in topology.scaled_pairs]
    seen = [(a, b) for a, b, *_ in entries]
    if len(seen) != len(set(seen)):
        raise ValueError("an OPLS pair is assigned more than one nonbonded treatment")
    return tuple(entries)


def rb_coefficients(fourier_kcal: tuple[float, float, float]) -> tuple[float, ...]:
    """Convert OPLS Fourier V1-V3 (kcal/mol) to GROMACS RB C0-C5 (kJ/mol)."""
    v1, v2, v3 = (value * KCAL_TO_KJ for value in fourier_kcal)
    return (v2 + (v1 + v3) / 2.0, (-v1 + 3.0 * v3) / 2.0, -v2, -2.0 * v3, 0.0, 0.0)


def _chain_stand_in(topology: OplsTopology) -> _Chain:
    return _Chain(topology.atoms, topology.bonds)


@dataclass(frozen=True)
class OplsBuild:
    topology: OplsTopology
    mapped: GlycamMappedSystem


def build_opls_system(crystal: BulkCrystal, variant: str, chain_boundary: str) -> OplsBuild:
    """Map a crystal onto repeated OPLS chains with the shared cellulose mapper."""
    from cellulose_forge.openff_sage import _mapped_coordinates  # noqa: PLC0415

    first_chain = next(iter(crystal.chain_numbers))
    degree = len(
        {
            residue
            for chain, residue in zip(crystal.chain_numbers, crystal.residue_numbers, strict=True)
            if chain == first_chain
        }
    )
    topology = build_opls_topology(variant, degree, chain_boundary)
    mapped = _mapped_coordinates(crystal, _chain_stand_in(topology), chain_boundary)
    return OplsBuild(topology, mapped)


def _format_topology(topology: OplsTopology, molecule_name: str) -> str:
    lines = [
        f"; CelluloseForge {topology.variant} cellulose chain "
        f"(DP {topology.degree_polymerization}, {topology.chain_boundary})",
        "[ moleculetype ]",
        "; name  nrexcl",
        f"{molecule_name}  3",
        "",
        "[ atoms ]",
        "; nr type resnr res atom cgnr charge mass",
    ]
    for number, (atom, atom_type, charge) in enumerate(
        zip(topology.atoms, topology.atom_types, topology.charges, strict=True), start=1
    ):
        mass = _OPLS_TYPES[atom_type][2]
        lines.append(
            f"{number:6d} {atom_type:>9s} {atom.residue:5d}  BGC {atom.name:>5s} "
            f"{number:6d} {charge:9.4f} {mass:9.4f}"
        )
    lines += ["", "[ bonds ]", "; explicit OPLS-AA parameters (b0 nm, kb kJ mol-1 nm-2)"]
    for a, b in topology.bonds:
        length, constant = bond_parameters(topology, a, b)
        lines.append(f"{a + 1:6d} {b + 1:6d} 1 {length:.4f} {constant:.1f}")
    lines += ["", "[ pairs ]", "; 1-4 pairs (OPLS fudgeQQ = fudgeLJ = 0.5)"]
    lines += [f"{a + 1:6d} {b + 1:6d} 1" for a, b in topology.pairs_14]
    if topology.scaled_pairs:
        lines.append("; OPLS-AA-SEI scaled 1-5/1-6 electrostatics, full Lennard-Jones")
        for a, b, fudge, kind in topology.scaled_pairs:
            sigma = float(np.sqrt(_LJ[topology.atom_types[a]][0] * _LJ[topology.atom_types[b]][0]))
            epsilon = float(
                np.sqrt(_LJ[topology.atom_types[a]][1] * _LJ[topology.atom_types[b]][1])
            )
            lines.append(
                f"{a + 1:6d} {b + 1:6d} 2 {fudge:.8f} {topology.charges[a]:8.4f} "
                f"{topology.charges[b]:8.4f} {sigma:.6f} {epsilon:.6f} ; {kind}"
            )
    lines += ["", "[ angles ]", "; explicit OPLS-AA parameters (theta0 deg, k kJ mol-1 rad-2)"]
    for a, b, c in topology.angles:
        theta, constant = angle_parameters(topology, a, b, c)
        lines.append(f"{a + 1:6d} {b + 1:6d} {c + 1:6d} 1 {theta:.3f} {constant:.3f}")
    lines += ["", "[ dihedrals ]", "; Ryckaert-Bellemans from Kony et al. 2002 Table 3"]
    table = _TORSIONS[topology.variant]
    for a, b, c, d, name in topology.dihedrals:
        coefficients = " ".join(f"{value:.10f}" for value in rb_coefficients(table[name]))
        lines.append(f"{a + 1:6d} {b + 1:6d} {c + 1:6d} {d + 1:6d} 3 {coefficients} ; {name}")
    if topology.scaled_pairs:
        lines += ["", "[ exclusions ]", "; SEI pairs replace the normal nonbonded term"]
        lines += [f"{a + 1:6d} {b + 1:6d}" for a, b, *_ in topology.scaled_pairs]
    return "\n".join(lines) + "\n"


def write_opls_system(build: OplsBuild, directory: Path) -> dict[str, str]:
    """Write ``structure.gro``, ``cellulose_opls.itp`` and ``system.top``."""
    directory.mkdir(parents=True, exist_ok=True)
    gro = directory / "structure.gro"
    itp = directory / "cellulose_opls.itp"
    top = directory / "system.top"
    write_gro(build.mapped, gro)
    itp.write_text(_format_topology(build.topology, "CELLULOSE"), encoding="utf-8")
    top.write_text(
        "; CelluloseForge experimental OPLS cellulose system\n"
        '#include "oplsaa.ff/forcefield.itp"\n'
        '#include "cellulose_opls.itp"\n\n'
        "[ system ]\nCelluloseForge cellulose\n\n"
        f"[ molecules ]\nCELLULOSE  {build.mapped.chain_count}\n",
        encoding="utf-8",
    )
    return {
        "gromacs_coordinates": gro.name,
        "gromacs_molecule_topology": itp.name,
        "gromacs_topology": top.name,
    }


def topology_report(topology: OplsTopology) -> dict[str, object]:
    """Machine-readable record of the parameter sources and choices."""
    one_five = sum(1 for *_, kind in topology.scaled_pairs if kind == "1-5")
    one_six = sum(1 for *_, kind in topology.scaled_pairs if kind == "1-6")
    return {
        "variant": topology.variant,
        "nonbonded_bonds_angles_source": (
            "OPLS-AA carbohydrate types (Damm et al. 1997) from GROMACS oplsaa.ff"
        ),
        "torsion_source": (
            "Kony et al. 2002 Table 3, "
            + ("OPLS-AA-SEI" if topology.variant == "opls-aa-sei" else "OPLS-AA")
            + " column; OS-C-C-OS from the OPLS-AA polyether term"
        ),
        "torsion_classes_per_chain": topology.class_counts(),
        "sei_scaling": (
            {
                "one_five_divisor": SEI_ONE_FIVE_SCALE,
                "one_six_divisor": SEI_ONE_SIX_SCALE,
                "scope": "within each residue, as published",
                "one_five_pairs_per_chain": one_five,
                "one_six_pairs_per_chain": one_six,
                "one_six_applies_to": "HO6/O4 and O6/HO4 where O4 is a hydroxyl",
            }
            if topology.variant == "opls-aa-sei"
            else None
        ),
        "chain_charge_e": float(sum(topology.charges)),
        "validation_scope": VALIDATION_SCOPE,
    }


def write_opls_report(path: Path, report: dict[str, object]) -> None:
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
