"""Map crystal coordinates onto repeated beta-(1->4) chain templates.

This module is shared by the graph-based backends (OpenFF Sage and OPLS-AA).
The mapping itself is backend-independent: a backend describes one chain as an
ordered sequence of :class:`ChainTemplateAtom` entries (atom name and residue
number, residue 1 being the reducing end), and the mapper returns coordinates
for every crystal chain in exactly that atom order. Crystal atoms are copied
unchanged from the residue selected by tracing the O4->C1 glycosidic path.

Atoms absent from the crystal are placed from fixed reference geometry, which
is not polymer-specific:

* Finite-chain caps (O1/HO1 on the reducing end, HO4 on the non-reducing end)
  use glucose-derived reference coordinates, superposed on the C1/H1/O5 and
  O4/C4/H4 anchors. Those anchor atoms are the same in GlcNAc, so the caps are
  also applied to chitin; a polymer with different anchors needs its own
  reference geometry.
* Cellulose II HO2/HO3/HO6 hydrogens, absent from the crystal definition, are
  placed from the same reference residue and resolved as one hydrogen-bond
  network. This applies only to cellulose II.

The returned :class:`GlycamMappedSystem` is the shared mapped-coordinate
container used by every backend despite its historical name.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import numpy as np

from cellulose_forge.cellulose_ii import cellulose_ii_preferred_acceptor
from cellulose_forge.crystallography import FloatArray
from cellulose_forge.geometry import BulkCrystal
from cellulose_forge.glycam_mapping import GlycamMappedSystem
from cellulose_forge.hydroxyl import HydroxylHydrogenPlacer

CHAIN_BOUNDARIES = frozenset({"finite", "periodic"})

_MAX_LINK_ANGSTROM = 2.0
_MINIMUM_CHAIN_RESIDUES = 2
_LINK_ATOMS = ("C1", "O4")
_BOX_DEGENERACY_TOLERANCE = 1.0e-12
_MISSING_II_HYDROGENS = frozenset({"HO2", "HO3", "HO6"})
_REDUCING_CAP = frozenset({"O1", "HO1"})
_NONREDUCING_CAP = frozenset({"HO4"})

_INTERNAL_REFERENCE = {
    "C1": (2.720, 11.610, 12.966),
    "H1": (2.053, 11.447, 12.268),
    "O5": (3.942, 12.030, 12.413),
    "C5": (3.939, 13.328, 11.848),
    "H5": (3.382, 13.347, 11.042),
    "C6": (5.429, 13.601, 11.483),
    "H62": (5.837, 12.780, 11.164),
    "H61": (5.907, 13.890, 12.273),
    "O6": (5.535, 14.593, 10.478),
    "HO6": (6.378, 14.431, 10.000),
    "C4": (3.454, 14.380, 12.860),
    "H4": (4.094, 14.469, 13.594),
    "C3": (2.078, 13.975, 13.392),
    "H3": (1.435, 13.946, 12.653),
    "O3": (1.602, 14.853, 14.397),
    "HO3": (1.714, 15.767, 14.053),
    "C2": (2.194, 12.595, 14.004),
    "H2": (2.841, 12.639, 14.743),
    "O2": (0.944, 12.153, 14.530),
    "HO2": (1.042, 11.190, 14.695),
    "O4": (3.368, 15.606, 12.110),
}
_REDUCING_REFERENCE = {
    "HO1": (3.163, 4.521, 12.573),
    "O1": (3.479, 5.230, 12.007),
    "C1": (3.658, 6.420, 12.839),
    "H1": (4.325, 6.257, 13.537),
    "O5": (2.436, 6.840, 13.392),
}
_NONREDUCING_REFERENCE = {
    "O4": (3.368, 36.366, 12.110),
    "C4": (3.454, 35.140, 12.860),
    "H4": (4.094, 35.229, 13.594),
    "HO4": (2.949, 37.032, 12.679),
}

Transform = tuple[FloatArray, FloatArray, FloatArray]


@dataclass(frozen=True)
class ChainTemplateAtom:
    """One atom of a backend chain graph, in the backend's own atom order."""

    name: str
    residue_number: int


def chain_template(atoms: Iterable[tuple[str, int]]) -> tuple[ChainTemplateAtom, ...]:
    """Return a validated chain template from (atom name, residue number) pairs."""
    template = tuple(ChainTemplateAtom(str(name), int(residue)) for name, residue in atoms)
    if not template:
        raise ValueError("a chain template must contain at least one atom")
    seen: set[tuple[int, str]] = set()
    for atom in template:
        key = (atom.residue_number, atom.name)
        if key in seen:
            raise ValueError(f"duplicate chain-template atom {atom.name} in residue {key[0]}")
        seen.add(key)
    return template


def reduced_periodic_box(box: FloatArray) -> FloatArray:
    """Return an equivalent triclinic lattice in OpenMM reduced form."""
    a, b, c = np.asarray(box, dtype=np.float64).copy()
    if abs(a[0]) < _BOX_DEGENERACY_TOLERANCE or abs(b[1]) < _BOX_DEGENERACY_TOLERANCE:
        raise ValueError("periodic output requires a non-degenerate triclinic box")
    c -= b * round(float(c[1] / b[1]))
    c -= a * round(float(c[0] / a[0]))
    b -= a * round(float(b[0] / a[0]))
    reduced = np.vstack((a, b, c))
    reduced.setflags(write=False)
    return reduced


def minimum_image(delta: FloatArray, box: FloatArray) -> FloatArray:
    """Return the minimum-image form of a displacement in a triclinic box."""
    fractional = delta @ np.linalg.inv(box)
    fractional -= np.round(fractional)
    return np.asarray(fractional @ box, dtype=np.float64)


def map_crystal_to_chain_template(  # noqa: PLR0912, PLR0915 - one pass keeps atom order explicit
    crystal: BulkCrystal,
    template: Sequence[ChainTemplateAtom],
    chain_boundary: str,
    residue_name: str,
) -> GlycamMappedSystem:
    """Map every crystal chain onto ``template`` and return the combined system.

    Residue 1 of the template is the reducing end. Each crystal chain is
    traversed along its O4->C1 glycosidic links, so the template residue order
    does not depend on the crystal's residue numbering. The mapped system lists
    chains in crystal order and, within each chain, atoms in template order.
    """
    if chain_boundary not in CHAIN_BOUNDARIES:
        raise ValueError("chain boundary must be finite or periodic")
    if not crystal.chain_numbers or not crystal.residue_numbers:
        raise ValueError("chain mapping requires explicit chain and residue numbers")
    periodic = chain_boundary == "periodic"
    if periodic and any(atom.name in _REDUCING_CAP | _NONREDUCING_CAP for atom in template):
        raise ValueError("periodic chain templates must not contain terminal cap atoms")
    template_residues = sorted({atom.residue_number for atom in template})
    if template_residues != list(range(1, len(template_residues) + 1)):
        raise ValueError("chain-template residues must be numbered 1..DP")
    chain_numbers = tuple(dict.fromkeys(crystal.chain_numbers))
    output: list[FloatArray] = []
    names: list[str] = []
    residues: list[int] = []
    chains: list[int] = []
    placer = (
        HydroxylHydrogenPlacer(
            crystal.atom_labels, crystal.cartesian_coordinates, crystal.basis_vectors
        )
        if crystal.allomorph == "II"
        else None
    )
    for chain in chain_numbers:
        order = _residue_order(crystal, chain, periodic)
        if len(order) != len(template_residues):
            raise ValueError(
                f"chain {chain} has {len(order)} residues but the template has "
                f"{len(template_residues)}"
            )
        targets = {
            index + 1: _lookup(crystal, chain, residue) for index, residue in enumerate(order)
        }
        reducing = (
            None
            if periodic
            else _rigid_transform(_REDUCING_REFERENCE, targets[1], ("C1", "H1", "O5"))
        )
        nonreducing = (
            None
            if periodic
            else _rigid_transform(_NONREDUCING_REFERENCE, targets[len(order)], ("O4", "C4", "H4"))
        )
        hydrogen_transforms = (
            {
                residue: _rigid_transform(_INTERNAL_REFERENCE, target, ("C1", "C3", "C5"))
                for residue, target in targets.items()
            }
            if placer is not None
            else {}
        )
        for atom in template:
            residue, name = atom.residue_number, atom.name
            target = targets[residue]
            if name in _REDUCING_CAP:
                if reducing is None or residue != 1:
                    raise ValueError(f"{name} is only valid on the finite reducing-end residue")
                coordinate = _transform(_REDUCING_REFERENCE[name], reducing)
            elif name in _NONREDUCING_CAP:
                if nonreducing is None or residue != len(order):
                    raise ValueError(f"{name} is only valid on the finite non-reducing residue")
                coordinate = _transform(_NONREDUCING_REFERENCE[name], nonreducing)
            elif placer is not None and name in _MISSING_II_HYDROGENS:
                oxygen = "O" + name[2]
                carbon = "C" + name[2]
                aligned_h = _transform(_INTERNAL_REFERENCE[name], hydrogen_transforms[residue])
                placer.request(
                    len(output),
                    target[carbon],
                    target[oxygen],
                    aligned_h,
                    donor=oxygen,
                    preferred_acceptor="O5" if oxygen == "O3" else None,
                    preferred_acceptor_coordinate=cellulose_ii_preferred_acceptor(
                        crystal,
                        chain_number=chain,
                        residue_number=order[residue - 1],
                        donor=oxygen,
                    ),
                )
                coordinate = np.zeros(3)
            elif name in target:
                coordinate = target[name]
            else:
                raise ValueError(
                    f"template atom {name} of residue {residue} is absent from crystal "
                    f"chain {chain}"
                )
            output.append(coordinate)
            names.append(name)
            residues.append(residue)
            chains.append(chain)
    if placer is not None:
        for index, coordinate in placer.resolve().items():
            output[index] = coordinate
    coordinates = np.vstack(output)
    coordinates.setflags(write=False)
    return GlycamMappedSystem(
        atom_names=tuple(names),
        residue_names=tuple(residue_name for _ in names),
        residue_numbers=tuple(residues),
        chain_numbers=tuple(chains),
        cartesian_coordinates=coordinates,
        basis_vectors=crystal.basis_vectors,
    )


def _lookup(crystal: BulkCrystal, chain: int, residue: int) -> dict[str, FloatArray]:
    return {
        name: xyz
        for name, xyz, atom_chain, atom_residue in zip(
            crystal.atom_labels,
            crystal.cartesian_coordinates,
            crystal.chain_numbers,
            crystal.residue_numbers,
            strict=True,
        )
        if atom_chain == chain and atom_residue == residue
    }


def _residue_order(  # noqa: PLR0912 - each refusal names the failing residue
    crystal: BulkCrystal, chain: int, periodic: bool
) -> tuple[int, ...]:
    residues = tuple(
        dict.fromkeys(
            residue
            for atom_chain, residue in zip(
                crystal.chain_numbers, crystal.residue_numbers, strict=True
            )
            if atom_chain == chain
        )
    )
    if len(residues) < _MINIMUM_CHAIN_RESIDUES:
        raise ValueError(f"chain {chain} has fewer than two residues; no glycosidic path exists")
    lookups = {residue: _lookup(crystal, chain, residue) for residue in residues}
    for residue, lookup in lookups.items():
        missing = [name for name in _LINK_ATOMS if name not in lookup]
        if missing:
            raise ValueError(
                f"chain {chain} residue {residue} lacks {' and '.join(missing)}; "
                "the glycosidic O4->C1 path cannot be traced"
            )
    successors: dict[int, int] = {}
    for residue in residues:
        partners = []
        for other in residues:
            if other == residue:
                continue
            delta = lookups[other]["C1"] - lookups[residue]["O4"]
            if periodic:
                delta = minimum_image(delta, crystal.basis_vectors)
            distance = float(np.linalg.norm(delta))
            if distance <= _MAX_LINK_ANGSTROM:
                partners.append((distance, other))
        if len(partners) > 1:
            listed = ", ".join(
                f"residue {other} ({distance:.2f} A)" for distance, other in partners
            )
            raise ValueError(
                f"chain {chain} residue {residue} O4 has more than one C1 within "
                f"{_MAX_LINK_ANGSTROM:.1f} A ({listed}); the glycosidic path is ambiguous"
            )
        if partners:
            successors[residue] = partners[0][1]
    targets = list(successors.values())
    repeated = sorted({other for other in targets if targets.count(other) > 1})
    if repeated:
        raise ValueError(
            f"chain {chain} residue(s) {', '.join(map(str, repeated))} receive more than one "
            "O4->C1 link; the glycosidic path is ambiguous"
        )
    if periodic:
        start = residues[0]
    else:
        predecessors = set(successors.values())
        starts = [residue for residue in residues if residue not in predecessors]
        if len(starts) != 1:
            raise ValueError(f"chain {chain} does not form one finite polymer path")
        start = starts[0]
    order = [start]
    while len(order) < len(residues) and order[-1] in successors:
        following = successors[order[-1]]
        if following in order:
            break
        order.append(following)
    if len(order) != len(residues):
        raise ValueError(f"chain {chain} glycosidic connectivity is incomplete")
    if periodic and successors.get(order[-1]) != order[0]:
        raise ValueError(f"chain {chain} glycosidic connectivity is not cyclic")
    return tuple(order)


def _rigid_transform(
    reference: dict[str, tuple[float, float, float]],
    target: dict[str, FloatArray],
    anchors: tuple[str, str, str],
) -> Transform:
    source = np.vstack([reference[name] for name in anchors])
    destination = np.vstack([target[name] for name in anchors])
    source_center = source.mean(axis=0)
    destination_center = destination.mean(axis=0)
    left, _, right = np.linalg.svd((source - source_center).T @ (destination - destination_center))
    rotation = left @ right
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right
    return rotation, source_center, destination_center


def _transform(coordinate: tuple[float, float, float], transform: Transform) -> FloatArray:
    rotation, source, destination = transform
    return np.asarray((np.asarray(coordinate) - source) @ rotation + destination, dtype=np.float64)
