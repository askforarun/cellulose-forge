"""Fail-closed native topology for CHARMM Drude polarizable cellulose.

Only the audited ``RESI BGLC`` and ``PRES 14bb`` definitions of the official
CHARMM Drude 2023a streams are read.  Every directive inside those two
definitions is either implemented or refused; nothing is silently skipped
except records that carry no topology (``GROUP``, ``IC``, ``DONOR``,
``ACCEPTOR`` and ``PATCH FIRST NONE LAST NONE``).

The chain builder applies ``14bb`` along a finite or cyclic chain and then
performs the operations CHARMM runs after patching (``AUTOGENERATE ANGLES
DIHEDRALS DRUDE``): it inserts one Drude particle after every polarizable
atom, moves the Drude charge and mass off the parent, and generates angles and
proper dihedrals over the real atoms.  Lone pairs and Drude particles never
enter an angle or dihedral.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

DRUDE_MAIN_FILE = "toppar_drude_main_protein_2023a.str"
DRUDE_CARBOHYDRATE_FILE = "toppar_drude_carbohydrate_2023a.str"

# SHA-256 of the stream files whose BGLC and 14bb definitions were audited.
# The audit holds for these bytes only, so any other file is refused.
AUDITED_STREAM_SHA256 = {
    DRUDE_MAIN_FILE: "083fa30a7bb00fdf0737a0bcfceb687c4939da0fb53047ab23054ca1016836ab",
    DRUDE_CARBOHYDRATE_FILE: "bfd45af59bce51742aa3e6b5075089e526e63befc86942471168804763d69cfb",
}

# CHARMM's Coulomb constant (kcal A / mol / e^2).  CHARMM derives the Drude
# charge from the polarizability with this value, and the acceptance reference
# for this backend is a CHARMM-generated PSF.  psfgen uses NAMD's 332.0636,
# which changes each Drude charge by about 1.2e-5 relative.
CHARMM_COULOMB_CONSTANT = 332.0716
NAMD_COULOMB_CONSTANT = 332.0636
DRUDE_MASS = 0.4
DRUDE_TYPE = "DRUD"
LONE_PAIR_TYPE = "LPD"
LONE_PAIR_MODES = frozenset({"relative", "bisector"})

_BLOCK_HEADER_FIELDS = 3
_ATOM_FIELDS = 4
_LONE_PAIR_FIELDS = 12
_ANISOTROPY_FIELDS = 9
_DELETE_ATOM_FIELDS = 3
_MASS_FIELDS = 4
_BOND_PARAMETER_FIELDS = 4
_ANGLE_PARAMETER_FIELDS = 5
_DIHEDRAL_PARAMETER_FIELDS = 7
_MINIMUM_PATCH_REFERENCE_LENGTH = 2
_MINIMUM_CHAIN_DP = 2
_CHARGE_TOLERANCE = 1.0e-8
_ISOTROPIC_TRACE = 3.0
_ONE_THREE_GRAPH_DISTANCE = 2

# Records inside a selected definition that define no topology here.
_NON_TOPOLOGY_DIRECTIVES = frozenset({"GROUP", "IC", "DONOR", "ACCEPTOR"})
_PARAMETER_SECTIONS = {
    "BOND": "bonds",
    "ANGL": "angles",
    "THET": "angles",
    "DIHE": "dihedrals",
    "PHI": "dihedrals",
}
_OTHER_SECTIONS = ("IMPR", "IMPH", "CMAP", "NONB", "NBON", "NBFI", "HBON", "THOL", "END", "READ")

AtomReference = tuple[int, str]


@dataclass(frozen=True)
class DrudeAtomDefinition:
    """One atom row from a Drude residue or patch."""

    name: str
    atom_type: str
    charge: float
    alpha: float | None = None
    thole: float | None = None


@dataclass(frozen=True)
class LonePairDefinition:
    """One CHARMM lone-pair construction record, as named in the stream."""

    mode: str
    site: str
    hosts: tuple[str, str, str]
    distance: float
    angle: float
    dihedral: float


@dataclass(frozen=True)
class AnisotropyDefinition:
    """One CHARMM anisotropic-polarizability record, as named in the stream."""

    atoms: tuple[str, str, str, str]
    a11: float
    a22: float


@dataclass(frozen=True)
class DrudeTopologyBlock:
    """Parsed topology operations for one selected RESI or PRES block."""

    name: str
    charge: float
    atoms: tuple[DrudeAtomDefinition, ...]
    bonds: tuple[tuple[str, str], ...]
    deleted_atoms: tuple[str, ...]
    lone_pairs: tuple[LonePairDefinition, ...]
    anisotropies: tuple[AnisotropyDefinition, ...]


@dataclass(frozen=True)
class DrudeParameterCoverage:
    """Bonded parameter keys defined by the streams, for coverage checks only."""

    bonds: frozenset[tuple[str, str]]
    angles: frozenset[tuple[str, str, str]]
    dihedrals: frozenset[tuple[str, str, str, str]]
    drude_force_constant: float


@dataclass(frozen=True)
class DrudeCelluloseDefinitions:
    """Masses, selected residue/linkage definitions, and parameter coverage."""

    masses: dict[str, float]
    residue: DrudeTopologyBlock
    patch: DrudeTopologyBlock
    parameters: DrudeParameterCoverage


@dataclass(frozen=True)
class DrudeSite:
    """One particle of a patched chain, in PSF order."""

    residue: int
    name: str
    atom_type: str
    charge: float
    mass: float
    kind: str  # "atom", "drude" or "lonepair"
    alpha: float = 0.0
    thole: float = 0.0


@dataclass(frozen=True)
class DrudeLonePair:
    """One lone-pair site and its three hosts, as chain site indices."""

    mode: str
    site: int
    hosts: tuple[int, int, int]
    distance: float
    angle: float
    dihedral: float


@dataclass(frozen=True)
class DrudeAnisotropy:
    """One anisotropic Drude response, as chain site indices and CHARMM constants."""

    atoms: tuple[int, int, int, int]
    a11: float
    a22: float
    k11: float
    k22: float
    k33: float


@dataclass(frozen=True)
class DrudeChainTopology:
    """Complete Drude topology of one patched cellulose chain.

    All index tuples refer to ``sites``.  Bonds include the host--lone-pair
    bonds and, after every other bond, the parent--Drude bonds.
    """

    degree_polymerization: int
    periodic: bool
    sites: tuple[DrudeSite, ...]
    bonds: tuple[tuple[int, int], ...]
    angles: tuple[tuple[int, int, int], ...]
    dihedrals: tuple[tuple[int, int, int, int], ...]
    lone_pairs: tuple[DrudeLonePair, ...]
    anisotropies: tuple[DrudeAnisotropy, ...]
    drude_parents: tuple[tuple[int, int], ...]
    coulomb_constant: float
    drude_force_constant: float

    @property
    def total_charge(self) -> float:
        """Return the net charge of every particle."""
        return math.fsum(site.charge for site in self.sites)

    @property
    def explicit_site_count(self) -> int:
        """Return the number of real atoms and lone pairs (no Drude particles)."""
        return sum(site.kind != "drude" for site in self.sites)

    @property
    def explicit_bond_count(self) -> int:
        """Return the number of bonds that do not involve a Drude particle."""
        return sum(
            self.sites[first].kind != "drude" and self.sites[second].kind != "drude"
            for first, second in self.bonds
        )

    def index(self, residue: int, name: str) -> int:
        """Return the site index of one named particle."""
        for index, site in enumerate(self.sites):
            if site.residue == residue and site.name == name:
                return index
        raise KeyError((residue, name))


def sha256_file(path: Path) -> str:
    """Return the SHA-256 of one file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _content(raw: str) -> str:
    return raw.split("!", maxsplit=1)[0].strip()


def _selected_block(lines: list[str], kind: str, name: str) -> tuple[float, list[str]]:
    starts = [
        index
        for index, raw in enumerate(lines)
        if len(fields := _content(raw).split()) >= _BLOCK_HEADER_FIELDS
        and fields[0].upper() == kind
        and fields[1].upper() == name
    ]
    if not starts:
        raise ValueError(f"Drude carbohydrate stream does not define {kind} {name}")
    if len(starts) > 1:
        raise ValueError(f"Drude carbohydrate stream defines {kind} {name} more than once")
    charge = float(_content(lines[starts[0]]).split()[2])
    selected: list[str] = []
    for raw in lines[starts[0] + 1 :]:
        fields = _content(raw).split()
        if fields and fields[0].upper() in {"RESI", "PRES"}:
            break
        selected.append(raw)
    return charge, selected


def _atom(fields: list[str]) -> DrudeAtomDefinition:
    if len(fields) < _ATOM_FIELDS:
        raise ValueError("malformed Drude ATOM directive")
    alpha: float | None = None
    thole: float | None = None
    tail = fields[4:]
    if len(tail) % 2:
        raise ValueError(f"malformed Drude ATOM options for {fields[1]}")
    for key, value in zip(tail[::2], tail[1::2], strict=True):
        if key.upper() == "ALPHA":
            alpha = float(value)
        elif key.upper() == "THOLE":
            thole = float(value)
        else:
            raise ValueError(f"unsupported Drude ATOM option {key!r} for {fields[1]}")
    if (alpha is None) != (thole is None):
        raise ValueError(f"Drude atom {fields[1]} must define ALPHA and THOLE together")
    return DrudeAtomDefinition(fields[1], fields[2], float(fields[3]), alpha, thole)


def _lone_pair(fields: list[str]) -> LonePairDefinition:
    if (
        len(fields) != _LONE_PAIR_FIELDS
        or fields[6].lower() != "distance"
        or fields[8].lower() != "angle"
        or fields[10].lower() != "dihe"
    ):
        raise ValueError("malformed Drude LONEPAIR directive")
    mode = fields[1].lower()
    if mode not in LONE_PAIR_MODES:
        raise ValueError(f"unsupported Drude LONEPAIR mode {fields[1]!r}")
    return LonePairDefinition(
        mode,
        fields[2],
        (fields[3], fields[4], fields[5]),
        float(fields[7]),
        float(fields[9]),
        float(fields[11]),
    )


def _anisotropy(fields: list[str]) -> AnisotropyDefinition:
    if (
        len(fields) != _ANISOTROPY_FIELDS
        or fields[5].upper() != "A11"
        or fields[7].upper() != "A22"
    ):
        raise ValueError("malformed Drude ANISOTROPY directive")
    return AnisotropyDefinition(
        (fields[1], fields[2], fields[3], fields[4]), float(fields[6]), float(fields[8])
    )


def _parse_block(  # noqa: PLR0912 - one explicit branch per supported directive
    kind: str, name: str, charge: float, lines: list[str]
) -> DrudeTopologyBlock:
    atoms: list[DrudeAtomDefinition] = []
    bonds: list[tuple[str, str]] = []
    deleted_atoms: list[str] = []
    lone_pairs: list[LonePairDefinition] = []
    anisotropies: list[AnisotropyDefinition] = []
    for raw in lines:
        content = _content(raw)
        if not content:
            continue
        if content.endswith("-"):
            raise ValueError(f"continuation lines are not supported in {kind} {name}: {content}")
        fields = content.split()
        directive = fields[0].upper()
        if directive in _NON_TOPOLOGY_DIRECTIVES:
            continue
        if directive == "PATCH":
            if [field.upper() for field in fields[1:]] != ["FIRST", "NONE", "LAST", "NONE"]:
                raise ValueError(f"unsupported terminal patch in {kind} {name}: {content}")
        elif directive == "ATOM":
            atoms.append(_atom(fields))
        elif directive == "BOND":
            if len(fields) < _BLOCK_HEADER_FIELDS or len(fields[1:]) % 2:
                raise ValueError(f"malformed BOND directive in {kind} {name}")
            bonds.extend(zip(fields[1::2], fields[2::2], strict=True))
        elif directive in {"DELE", "DELETE"}:
            if len(fields) == _DELETE_ATOM_FIELDS and fields[1].upper() == "ATOM":
                deleted_atoms.append(fields[2])
            else:
                raise ValueError(f"unsupported deletion in {kind} {name}: {content}")
        elif directive == "LONEPAIR":
            lone_pairs.append(_lone_pair(fields))
        elif directive == "ANISOTROPY":
            anisotropies.append(_anisotropy(fields))
        else:
            raise ValueError(f"unsupported directive {directive!r} in {kind} {name}")
    names = [atom.name for atom in atoms]
    if len(names) != len(set(names)):
        raise ValueError(f"{kind} {name} defines an atom more than once")
    return DrudeTopologyBlock(
        name,
        charge,
        tuple(atoms),
        tuple(bonds),
        tuple(deleted_atoms),
        tuple(lone_pairs),
        tuple(anisotropies),
    )


def _is_number(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True


def _parse_parameter_coverage(texts: tuple[str, ...]) -> DrudeParameterCoverage:
    bonds: set[tuple[str, str]] = set()
    angles: set[tuple[str, str, str]] = set()
    dihedrals: set[tuple[str, str, str, str]] = set()
    drude_constants: set[float] = set()
    for text in texts:
        section: str | None = None
        for raw in text.splitlines():
            fields = _content(raw).split()
            if not fields or fields[0].startswith("*"):
                continue
            head = fields[0].upper()
            parameter_section = next(
                (value for key, value in _PARAMETER_SECTIONS.items() if head.startswith(key)),
                None,
            )
            if parameter_section is not None and (len(fields) == 1 or not _is_number(fields[-1])):
                section = parameter_section
                continue
            if any(head.startswith(key) for key in _OTHER_SECTIONS):
                section = None
                continue
            if (
                section == "bonds"
                and len(fields) >= _BOND_PARAMETER_FIELDS
                and _is_number(fields[2])
            ):
                first, second = fields[0].upper(), fields[1].upper()
                bonds.add((first, second))
                bonds.add((second, first))
                if {first, second} == {"X", DRUDE_TYPE}:
                    drude_constants.add(float(fields[2]))
            elif (
                section == "angles"
                and len(fields) >= _ANGLE_PARAMETER_FIELDS
                and _is_number(fields[3])
            ):
                key = (fields[0].upper(), fields[1].upper(), fields[2].upper())
                angles.add(key)
                angles.add(key[::-1])
            elif (
                section == "dihedrals"
                and len(fields) >= _DIHEDRAL_PARAMETER_FIELDS
                and _is_number(fields[4])
            ):
                quad = (
                    fields[0].upper(),
                    fields[1].upper(),
                    fields[2].upper(),
                    fields[3].upper(),
                )
                dihedrals.add(quad)
                dihedrals.add(quad[::-1])
    if len(drude_constants) != 1:
        raise ValueError(
            "Drude streams must define exactly one 'X DRUD' bond force constant; "
            f"found {sorted(drude_constants)}"
        )
    return DrudeParameterCoverage(
        frozenset(bonds), frozenset(angles), frozenset(dihedrals), drude_constants.pop()
    )


def parse_drude_cellulose_definitions(
    main_text: str, carbohydrate_text: str
) -> DrudeCelluloseDefinitions:
    """Parse the BGLC/14bb subset from stream text, without the audit pin.

    Use :func:`read_drude_cellulose_definitions` for real stream files; this
    entry point exists so the parser can be tested on small fixtures.
    """
    masses: dict[str, float] = {}
    for text in (main_text, carbohydrate_text):
        for raw in text.splitlines():
            fields = _content(raw).split()
            if fields and fields[0].upper() == "MASS" and len(fields) >= _MASS_FIELDS:
                masses[fields[2].upper()] = float(fields[3])
    carbohydrate_lines = carbohydrate_text.splitlines()
    residue_charge, residue_lines = _selected_block(carbohydrate_lines, "RESI", "BGLC")
    patch_charge, patch_lines = _selected_block(carbohydrate_lines, "PRES", "14BB")
    residue = _parse_block("RESI", "BGLC", residue_charge, residue_lines)
    patch = _parse_block("PRES", "14BB", patch_charge, patch_lines)
    used_types = {atom.atom_type.upper() for atom in (*residue.atoms, *patch.atoms)}
    missing = sorted((used_types | {DRUDE_TYPE, LONE_PAIR_TYPE}).difference(masses))
    if missing:
        raise ValueError("Drude stream lacks masses for atom types: " + ", ".join(missing))
    if masses[DRUDE_TYPE] != 0.0 or masses[LONE_PAIR_TYPE] != 0.0:
        raise ValueError("Drude streams must define massless DRUD and LPD atom types")
    parameters = _parse_parameter_coverage((main_text, carbohydrate_text))
    return DrudeCelluloseDefinitions(masses, residue, patch, parameters)


def verify_audited_streams(
    main_stream: str | Path, carbohydrate_stream: str | Path
) -> dict[str, str]:
    """Return the stream checksums, refusing any file that was not audited."""
    checksums: dict[str, str] = {}
    for path, name in (
        (Path(main_stream), DRUDE_MAIN_FILE),
        (Path(carbohydrate_stream), DRUDE_CARBOHYDRATE_FILE),
    ):
        if not path.is_file():
            raise ValueError(f"CHARMM Drude stream not found: {path}")
        checksum = sha256_file(path)
        if checksum != AUDITED_STREAM_SHA256[name]:
            raise ValueError(
                f"{path.name} (SHA-256 {checksum}) is not the audited CHARMM Drude 2023a "
                f"{name} (SHA-256 {AUDITED_STREAM_SHA256[name]}); the BGLC and 14bb audit "
                "applies to that file only"
            )
        checksums[name] = checksum
    return checksums


def read_drude_cellulose_definitions(
    main_stream: str | Path, carbohydrate_stream: str | Path
) -> DrudeCelluloseDefinitions:
    """Read the audited BGLC/14bb subset of the official CHARMM Drude streams."""
    verify_audited_streams(main_stream, carbohydrate_stream)
    return parse_drude_cellulose_definitions(
        Path(main_stream).read_text(encoding="utf-8", errors="replace"),
        Path(carbohydrate_stream).read_text(encoding="utf-8", errors="replace"),
    )


def _patch_reference(value: str, first: int, second: int) -> AtomReference:
    if len(value) < _MINIMUM_PATCH_REFERENCE_LENGTH or value[0] not in {"1", "2"}:
        raise ValueError(f"14bb atom reference lacks a residue selector: {value}")
    return (first if value[0] == "1" else second), value[1:]


def drude_charge(alpha: float, force_constant: float, coulomb_constant: float) -> float:
    """Return the Drude particle charge CHARMM assigns for one polarizability.

    The stream gives negative ``ALPHA`` values, for which CHARMM places a
    negative charge on the Drude particle.  Positive values are refused because
    no audited definition uses them.
    """
    if alpha >= 0.0:
        raise ValueError(f"only negative CHARMM ALPHA values are supported, not {alpha}")
    return -math.sqrt(2.0 * force_constant * abs(alpha) / coulomb_constant)


def anisotropy_constants(
    a11: float, a22: float, force_constant: float
) -> tuple[float, float, float]:
    """Return CHARMM's K11, K22 and K33 for one ANISOTROPY record.

    ``A33 = 3 - A11 - A22``; each ``K = k_D / A - k_D``, and K33 is then removed
    from K11 and K22 because CHARMM adds it isotropically.
    """
    a33 = _ISOTROPIC_TRACE - a11 - a22
    if min(a11, a22, a33) <= 0.0:
        raise ValueError(f"anisotropy components must be positive: {a11}, {a22}, {a33}")
    k33 = force_constant / a33 - force_constant
    return (
        force_constant / a11 - force_constant - k33,
        force_constant / a22 - force_constant - k33,
        k33,
    )


def _patched_chain(  # noqa: PLR0912 - patch operations stay auditable in order
    definitions: DrudeCelluloseDefinitions, degree_polymerization: int, periodic: bool
) -> tuple[
    dict[int, list[str]],
    dict[AtomReference, DrudeAtomDefinition],
    set[tuple[AtomReference, AtomReference]],
    dict[AtomReference, tuple[str, tuple[AtomReference, ...], float, float, float]],
    dict[tuple[AtomReference, ...], tuple[float, float]],
]:
    residues = range(1, degree_polymerization + 1)
    order = {residue: [atom.name for atom in definitions.residue.atoms] for residue in residues}
    atoms = {
        (residue, atom.name): atom for residue in residues for atom in definitions.residue.atoms
    }
    bonds: set[tuple[AtomReference, AtomReference]] = set()
    lone_pairs: dict[AtomReference, tuple[str, tuple[AtomReference, ...], float, float, float]] = {}
    anisotropies: dict[tuple[AtomReference, ...], tuple[float, float]] = {}
    for residue in residues:
        for left, right in definitions.residue.bonds:
            bonds.add(tuple(sorted(((residue, left), (residue, right)))))  # type: ignore[arg-type]
        for lone_pair in definitions.residue.lone_pairs:
            lone_pairs[(residue, lone_pair.site)] = (
                lone_pair.mode,
                tuple((residue, host) for host in lone_pair.hosts),
                lone_pair.distance,
                lone_pair.angle,
                lone_pair.dihedral,
            )
        for anisotropy in definitions.residue.anisotropies:
            anisotropies[tuple((residue, name) for name in anisotropy.atoms)] = (
                anisotropy.a11,
                anisotropy.a22,
            )

    residue_alpha = {atom.name: atom.alpha for atom in definitions.residue.atoms}
    links = [(residue, residue + 1) for residue in range(1, degree_polymerization)]
    if periodic:
        links.append((degree_polymerization, 1))
    for first, second in links:
        deleted: set[AtomReference] = set()
        for reference in definitions.patch.deleted_atoms:
            residue, name = _patch_reference(reference, first, second)
            if (residue, name) in atoms:
                del atoms[(residue, name)]
                order[residue].remove(name)
                deleted.add((residue, name))
                continue
            # CHARMM lists the generated Drude of a deleted polarizable atom
            # explicitly (``dele atom 2DO1``).  Drude particles are generated
            # after patching here, so that deletion is accepted only when the
            # same patch application deletes the polarizable parent.
            parent = (residue, name[1:])
            parent_deleted = any(
                _patch_reference(other, first, second) == parent
                for other in definitions.patch.deleted_atoms
            )
            if (
                not name.startswith("D")
                or residue_alpha.get(parent[1]) is None
                or not parent_deleted
            ):
                raise ValueError(f"14bb deletes absent atom {residue}:{name}")
        bonds = {bond for bond in bonds if not deleted.intersection(bond)}
        lone_pairs = {
            site: value
            for site, value in lone_pairs.items()
            if site not in deleted and not deleted.intersection(value[1])
        }
        anisotropies = {
            key: value for key, value in anisotropies.items() if not deleted.intersection(key)
        }
        for atom in definitions.patch.atoms:
            residue, name = _patch_reference(atom.name, first, second)
            if (residue, name) not in atoms:
                order[residue].append(name)
            atoms[(residue, name)] = DrudeAtomDefinition(
                name, atom.atom_type, atom.charge, atom.alpha, atom.thole
            )
        for left, right in definitions.patch.bonds:
            bonds.add(
                tuple(  # type: ignore[arg-type]
                    sorted(
                        (
                            _patch_reference(left, first, second),
                            _patch_reference(right, first, second),
                        )
                    )
                )
            )
        for lone_pair in definitions.patch.lone_pairs:
            lone_pairs[_patch_reference(lone_pair.site, first, second)] = (
                lone_pair.mode,
                tuple(_patch_reference(host, first, second) for host in lone_pair.hosts),
                lone_pair.distance,
                lone_pair.angle,
                lone_pair.dihedral,
            )
        for anisotropy in definitions.patch.anisotropies:
            key = tuple(_patch_reference(name, first, second) for name in anisotropy.atoms)
            parent = key[0]
            anisotropies = {
                other: value for other, value in anisotropies.items() if other[0] != parent
            }
            anisotropies[key] = (anisotropy.a11, anisotropy.a22)
    return order, atoms, bonds, lone_pairs, anisotropies


def _require_parameters(topology: DrudeChainTopology, parameters: DrudeParameterCoverage) -> None:
    """Refuse a topology whose real-atom bonded terms lack stream parameters."""
    types = [site.atom_type.upper() for site in topology.sites]
    kinds = [site.kind for site in topology.sites]
    missing: set[str] = set()
    for first, second in topology.bonds:
        if kinds[first] != "atom" or kinds[second] != "atom":
            continue
        if (types[first], types[second]) not in parameters.bonds:
            missing.add(f"bond {types[first]}-{types[second]}")
    for angle in topology.angles:
        key = (types[angle[0]], types[angle[1]], types[angle[2]])
        if key not in parameters.angles:
            missing.add("angle " + "-".join(key))
    for dihedral in topology.dihedrals:
        key4 = (
            types[dihedral[0]],
            types[dihedral[1]],
            types[dihedral[2]],
            types[dihedral[3]],
        )
        wildcard = ("X", key4[1], key4[2], "X")
        if key4 not in parameters.dihedrals and wildcard not in parameters.dihedrals:
            missing.add("dihedral " + "-".join(key4))
    if missing:
        raise ValueError(
            "CHARMM Drude streams define no parameters for: " + ", ".join(sorted(missing))
        )


def build_drude_cellulose_chain(  # noqa: PLR0912, PLR0915 - CHARMM's generation order is explicit
    definitions: DrudeCelluloseDefinitions,
    degree_polymerization: int,
    *,
    periodic: bool,
    coulomb_constant: float = CHARMM_COULOMB_CONSTANT,
) -> DrudeChainTopology:
    """Apply 14bb natively to one finite or cyclic chain and generate its Drude topology."""
    if degree_polymerization < _MINIMUM_CHAIN_DP:
        raise ValueError("CHARMM Drude cellulose chains require a DP of at least 2")
    force_constant = definitions.parameters.drude_force_constant
    order, atoms, named_bonds, named_lone_pairs, named_anisotropies = _patched_chain(
        definitions, degree_polymerization, periodic
    )
    for bond in named_bonds:
        absent = [reference for reference in bond if reference not in atoms]
        if absent:
            raise ValueError(f"patched Drude bond references absent atoms: {absent}")

    sites: list[DrudeSite] = []
    index_of: dict[AtomReference, int] = {}
    drude_parents: list[tuple[int, int]] = []
    for residue in range(1, degree_polymerization + 1):
        for name in order[residue]:
            atom = atoms[(residue, name)]
            mass = definitions.masses[atom.atom_type.upper()]
            is_lone_pair = (residue, name) in named_lone_pairs
            if is_lone_pair != (atom.atom_type.upper() == LONE_PAIR_TYPE):
                raise ValueError(
                    f"{residue}:{name} must be an {LONE_PAIR_TYPE} site with exactly one "
                    "LONEPAIR record"
                )
            index_of[(residue, name)] = len(sites)
            if is_lone_pair:
                if atom.alpha is not None:
                    raise ValueError(f"lone pair {residue}:{name} cannot be polarizable")
                sites.append(DrudeSite(residue, name, atom.atom_type, atom.charge, 0.0, "lonepair"))
                continue
            if atom.alpha is None or atom.thole is None:
                sites.append(DrudeSite(residue, name, atom.atom_type, atom.charge, mass, "atom"))
                continue
            if (residue, "D" + name) in atoms:
                raise ValueError(f"explicit atom {residue}:D{name} collides with a Drude name")
            charge = drude_charge(atom.alpha, force_constant, coulomb_constant)
            parent = len(sites)
            sites.append(
                DrudeSite(
                    residue,
                    name,
                    atom.atom_type,
                    atom.charge - charge,
                    mass - DRUDE_MASS,
                    "atom",
                    atom.alpha,
                    atom.thole,
                )
            )
            sites.append(DrudeSite(residue, "D" + name, DRUDE_TYPE, charge, DRUDE_MASS, "drude"))
            drude_parents.append((parent + 1, parent))

    explicit_bonds = sorted(
        tuple(sorted((index_of[first], index_of[second]))) for first, second in named_bonds
    )
    neighbours: dict[int, set[int]] = {}
    for first, second in explicit_bonds:
        if sites[first].kind == "atom" and sites[second].kind == "atom":
            neighbours.setdefault(first, set()).add(second)
            neighbours.setdefault(second, set()).add(first)
        elif "lonepair" not in {sites[first].kind, sites[second].kind}:
            raise ValueError("patched Drude bond joins unsupported particle kinds")
    angles = sorted(
        (first, centre, second)
        for centre, bonded in neighbours.items()
        for first, second in combinations(sorted(bonded), 2)
    )
    dihedrals: set[tuple[int, int, int, int]] = set()
    for second, bonded in neighbours.items():
        for third in bonded:
            if second > third:
                continue
            for first in bonded - {third}:
                for fourth in neighbours[third] - {second, first}:
                    quad = (first, second, third, fourth)
                    dihedrals.add(min(quad, quad[::-1]))

    lone_pairs: list[DrudeLonePair] = []
    for site_reference, (mode, hosts, distance, angle, dihedral) in named_lone_pairs.items():
        absent = [host for host in hosts if host not in index_of]
        if absent:
            raise ValueError(f"lone pair {site_reference} references absent hosts: {absent}")
        host_indices = tuple(index_of[host] for host in hosts)
        if any(sites[host].kind != "atom" for host in host_indices):
            raise ValueError(f"lone pair {site_reference} must be hosted by real atoms")
        lone_pairs.append(
            DrudeLonePair(
                mode,
                index_of[site_reference],
                host_indices,  # type: ignore[arg-type]
                distance,
                angle,
                dihedral,
            )
        )
    lone_pairs.sort(key=lambda item: item.site)

    anisotropies: list[DrudeAnisotropy] = []
    for references, (a11, a22) in named_anisotropies.items():
        absent = [reference for reference in references if reference not in index_of]
        if absent:
            raise ValueError(f"anisotropy references absent atoms: {absent}")
        indices = tuple(index_of[reference] for reference in references)
        if sites[indices[0]].alpha == 0.0:
            raise ValueError(f"anisotropic atom {references[0]} is not polarizable")
        anisotropies.append(
            DrudeAnisotropy(
                indices,  # type: ignore[arg-type]
                a11,
                a22,
                *anisotropy_constants(a11, a22, force_constant),
            )
        )
    anisotropies.sort(key=lambda item: item.atoms)
    parents = [item.atoms[0] for item in anisotropies]
    if len(parents) != len(set(parents)):
        raise ValueError("an atom carries more than one anisotropy record")

    # NAMD accepts a Drude PSF only when the parent--Drude bonds follow every
    # ordinary and lone-pair bond; mixed in, it tries to parameterize them
    # ("UNABLE TO FIND BOND PARAMETERS FOR CD31HA DRUD", NAMD 3.0.3).
    bonds = [*explicit_bonds, *sorted((parent, drude) for drude, parent in drude_parents)]
    topology = DrudeChainTopology(
        degree_polymerization,
        periodic,
        tuple(sites),
        tuple(bonds),  # type: ignore[arg-type]
        tuple(angles),
        tuple(sorted(dihedrals)),
        tuple(lone_pairs),
        tuple(anisotropies),
        tuple(drude_parents),
        coulomb_constant,
        force_constant,
    )
    if abs(topology.total_charge) > _CHARGE_TOLERANCE:
        raise ValueError(f"patched CHARMM Drude chain is not neutral ({topology.total_charge:.8f})")
    _require_parameters(topology, definitions.parameters)
    return topology


def thole_screened_pairs(topology: DrudeChainTopology) -> tuple[tuple[int, int, float], ...]:
    """Return the 1-2 and 1-3 polarizable parent pairs and their Thole factors.

    CHARMM screens the dipole--dipole interaction of polarizable atoms that are
    excluded from ordinary nonbonded interactions, with ``a_ij = a_i + a_j``.
    Each row is ``(parent_i, parent_j, a_ij)``.  NAMD derives these pairs from a
    Drude PSF itself; an engine that needs an explicit list (OpenMM
    ``DrudeForce.addScreenedPair``) uses this function.
    """
    neighbours: dict[int, set[int]] = {}
    for first, second in topology.bonds:
        if topology.sites[first].kind == "atom" and topology.sites[second].kind == "atom":
            neighbours.setdefault(first, set()).add(second)
            neighbours.setdefault(second, set()).add(first)
    polarizable = {parent for _drude, parent in topology.drude_parents}
    pairs: set[tuple[int, int]] = set()
    for atom in polarizable:
        frontier = {atom}
        seen = {atom}
        for _ in range(_ONE_THREE_GRAPH_DISTANCE):
            frontier = {other for item in frontier for other in neighbours.get(item, ())} - seen
            seen |= frontier
            pairs.update(
                (atom, other) for other in frontier if other in polarizable and atom < other
            )
    return tuple(
        (first, second, topology.sites[first].thole + topology.sites[second].thole)
        for first, second in sorted(pairs)
    )
