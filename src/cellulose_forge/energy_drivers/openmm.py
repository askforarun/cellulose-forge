"""Normalize classified OpenMM single-point components.

OpenMM force classes do not map one-to-one onto the canonical report.  In
particular, CHARMM Urey--Bradley terms are stored in a bond force, periodic
impropers can share a torsion force with propers, and dispersion corrections
are included inside nonbonded forces.  The evaluator therefore classifies and
isolates those contributions before calling this normalization layer.
"""

from __future__ import annotations

import hashlib
import importlib
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any, cast

from cellulose_forge.energy_drivers.report import EnergyReport, EnergyValue

OPENMM_COMPONENTS = frozenset(
    {
        "bond",
        "angle",
        "urey_bradley",
        "proper_dihedral",
        "periodic_improper",
        "harmonic_improper",
        "vdw",
        "vdw_14",
        "electrostatics",
        "electrostatics_14",
        "restraint",
        "dispersion_correction",
        "potential",
    }
)

_EXPLICIT_FORCE_ROLES = frozenset(
    {
        "bond",
        "angle",
        "urey_bradley",
        "proper_dihedral",
        "periodic_improper",
        "harmonic_improper",
        "vdw",
        "vdw_14",
        "electrostatics",
        "electrostatics_14",
        "restraint",
    }
)
_GRO_HEADER_AND_BOX_LINES = 3
_MIN_GRO_FIELD_WIDTH = 2
_ORTHORHOMBIC_BOX_VALUES = 3
_TRICLINIC_BOX_VALUES = 9


@dataclass(frozen=True)
class OpenMMForceClassification:
    """Explicit classification for OpenMM terms that force classes cannot identify.

    Keys are zero-based force indices in the serialized ``System``.  Interaction
    indices identify Urey--Bradley bonds or periodic impropers inside mixed
    ``HarmonicBondForce`` and ``PeriodicTorsionForce`` objects.  ``force_roles``
    is required for custom forces and may also classify an otherwise unsupported
    force.  Long-range corrections are isolated by the evaluator and must not be
    listed as a role.
    """

    urey_bradley_bonds: Mapping[int, tuple[int, ...]] = field(default_factory=dict)
    periodic_impropers: Mapping[int, tuple[int, ...]] = field(default_factory=dict)
    rb_impropers: Mapping[int, tuple[int, ...]] = field(default_factory=dict)
    force_roles: Mapping[int, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for description, mapping in (
            ("Urey--Bradley", self.urey_bradley_bonds),
            ("periodic improper", self.periodic_impropers),
            ("RB improper", self.rb_impropers),
        ):
            normalized: dict[int, tuple[int, ...]] = {}
            for force_index, interaction_indices in mapping.items():
                if type(force_index) is not int or force_index < 0:
                    raise ValueError(f"{description} force indices must be non-negative integers")
                values = tuple(interaction_indices)
                if any(type(index) is not int or index < 0 for index in values):
                    raise ValueError(
                        f"{description} interaction indices must be non-negative integers"
                    )
                if len(values) != len(set(values)):
                    raise ValueError(f"duplicate {description} interaction index")
                normalized[force_index] = values
            object.__setattr__(
                self,
                {
                    "Urey--Bradley": "urey_bradley_bonds",
                    "periodic improper": "periodic_impropers",
                    "RB improper": "rb_impropers",
                }[description],
                MappingProxyType(normalized),
            )
        roles: dict[int, str] = {}
        for force_index, role in self.force_roles.items():
            if type(force_index) is not int or force_index < 0:
                raise ValueError("explicit force indices must be non-negative integers")
            if role not in _EXPLICIT_FORCE_ROLES:
                raise ValueError(f"unsupported explicit OpenMM force role: {role!r}")
            roles[force_index] = role
        object.__setattr__(self, "force_roles", MappingProxyType(roles))


def _gro_coordinate_data(openmm: Any, unit: Any, path: Path) -> tuple[Any, Any]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < _GRO_HEADER_AND_BOX_LINES:
        raise ValueError("OpenMM GRO coordinate file is truncated")
    try:
        atom_count = int(lines[1].strip())
    except ValueError as error:
        raise ValueError("OpenMM GRO atom count is invalid") from error
    if len(lines) != atom_count + 3:
        raise ValueError("OpenMM GRO atom count does not match the coordinate rows")
    positions = []
    try:
        for line in lines[2 : 2 + atom_count]:
            first_decimal = line.index(".", 20)
            second_decimal = line.index(".", first_decimal + 1)
            field_width = second_decimal - first_decimal
            if field_width < _MIN_GRO_FIELD_WIDTH or len(line) < 20 + 3 * field_width:
                raise ValueError
            x, y, z = (
                float(line[20 + axis * field_width : 20 + (axis + 1) * field_width])
                for axis in range(3)
            )
            positions.append(openmm.Vec3(x, y, z))
        values = [float(value) for value in lines[-1].split()]
    except (ValueError, IndexError) as error:
        raise ValueError("OpenMM GRO coordinates are invalid") from error
    if len(values) == _ORTHORHOMBIC_BOX_VALUES:
        box = (
            openmm.Vec3(values[0], 0, 0),
            openmm.Vec3(0, values[1], 0),
            openmm.Vec3(0, 0, values[2]),
        )
    elif len(values) == _TRICLINIC_BOX_VALUES:
        box = (
            openmm.Vec3(values[0], values[3], values[4]),
            openmm.Vec3(values[5], values[1], values[6]),
            openmm.Vec3(values[7], values[8], values[2]),
        )
    else:
        raise ValueError("OpenMM GRO box must contain three or nine values")
    return positions * unit.nanometer, tuple(vector * unit.nanometer for vector in box)


def _coordinate_data(openmm: Any, app: Any, unit: Any, path: Path) -> tuple[Any, Any | None]:
    suffix = path.suffix.lower()
    if suffix == ".gro":
        return _gro_coordinate_data(openmm, unit, path)
    if suffix == ".pdb":
        coordinates = app.PDBFile(str(path))
        return coordinates.getPositions(), coordinates.getTopology().getPeriodicBoxVectors()
    if suffix in {".cif", ".mmcif"}:
        coordinates = app.PDBxFile(str(path))
        return coordinates.getPositions(), coordinates.getTopology().getPeriodicBoxVectors()
    raise ValueError("OpenMM shared coordinates must be .gro, .pdb, .cif, or .mmcif")


def _context_energy(  # noqa: PLR0913, PLR0917 - one explicit evaluation context
    openmm: Any, unit: Any, system: Any, positions: Any, box: Any, platform: Any
) -> float:
    context = openmm.Context(system, openmm.VerletIntegrator(0.001), platform)
    if box is not None:
        context.setPeriodicBoxVectors(*box)
    context.setPositions(positions)
    context.computeVirtualSites()
    energy = context.getState(getEnergy=True).getPotentialEnergy()
    result = float(energy.value_in_unit(unit.kilojoule_per_mole))
    del context
    return result


def _clone(openmm: Any, system: Any) -> Any:
    return openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(system))


def _one_force_system(openmm: Any, system: Any, force_index: int) -> Any:
    isolated = _clone(openmm, system)
    for index in reversed(range(isolated.getNumForces())):
        if index != force_index:
            isolated.removeForce(index)
    return isolated


def _zero_harmonic_bonds(force: Any, selected: set[int]) -> None:
    for index in range(force.getNumBonds()):
        first, second, length, constant = force.getBondParameters(index)
        if index in selected:
            force.setBondParameters(index, first, second, length, 0 * constant)


def _zero_periodic_torsions(force: Any, selected: set[int]) -> None:
    for index in range(force.getNumTorsions()):
        first, second, third, fourth, periodicity, phase, constant = force.getTorsionParameters(
            index
        )
        if index in selected:
            force.setTorsionParameters(
                index, first, second, third, fourth, periodicity, phase, 0 * constant
            )


def _zero_rb_torsions(force: Any, selected: set[int]) -> None:
    for index in range(force.getNumTorsions()):
        parameters = force.getTorsionParameters(index)
        if index in selected:
            particles = parameters[:4]
            coefficients = parameters[4:]
            force.setTorsionParameters(index, *particles, *(0 * value for value in coefficients))


def _nonbonded_components(  # noqa: PLR0913, PLR0917 - one explicit evaluation context
    openmm: Any, unit: Any, system: Any, force_index: int, positions: Any, box: Any, platform: Any
) -> dict[str, float]:
    electrostatic_system = _one_force_system(openmm, system, force_index)
    electrostatic = electrostatic_system.getForce(0)
    if (
        electrostatic.getNumParticleParameterOffsets()
        or electrostatic.getNumExceptionParameterOffsets()
    ):
        raise ValueError("OpenMM NonbondedForce parameter offsets are unsupported")
    for index in range(electrostatic.getNumParticles()):
        charge, sigma, epsilon = electrostatic.getParticleParameters(index)
        electrostatic.setParticleParameters(index, charge, sigma, 0 * epsilon)
    for index in range(electrostatic.getNumExceptions()):
        first, second, charge_product, sigma, epsilon = electrostatic.getExceptionParameters(index)
        electrostatic.setExceptionParameters(
            index, first, second, charge_product, sigma, 0 * epsilon
        )
    electrostatics = _context_energy(openmm, unit, electrostatic_system, positions, box, platform)

    vdw_system = _one_force_system(openmm, system, force_index)
    vdw_force = vdw_system.getForce(0)
    for index in range(vdw_force.getNumParticles()):
        charge, sigma, epsilon = vdw_force.getParticleParameters(index)
        vdw_force.setParticleParameters(index, 0 * charge, sigma, epsilon)
    for index in range(vdw_force.getNumExceptions()):
        first, second, charge_product, sigma, epsilon = vdw_force.getExceptionParameters(index)
        vdw_force.setExceptionParameters(index, first, second, 0 * charge_product, sigma, epsilon)
    enabled = _context_energy(openmm, unit, vdw_system, positions, box, platform)
    uses_correction = bool(vdw_force.getUseDispersionCorrection())
    if uses_correction:
        vdw_force.setUseDispersionCorrection(False)
    without = _context_energy(openmm, unit, vdw_system, positions, box, platform)
    return {
        "electrostatics": electrostatics,
        "vdw": without,
        "dispersion_correction": enabled - without,
    }


def _custom_nonbonded_components(  # noqa: PLR0913, PLR0917 - one explicit evaluation context
    openmm: Any, unit: Any, system: Any, force_index: int, positions: Any, box: Any, platform: Any
) -> dict[str, float]:
    isolated = _one_force_system(openmm, system, force_index)
    force = isolated.getForce(0)
    enabled = _context_energy(openmm, unit, isolated, positions, box, platform)
    uses_correction = bool(force.getUseLongRangeCorrection())
    if uses_correction:
        force.setUseLongRangeCorrection(False)
    without = _context_energy(openmm, unit, isolated, positions, box, platform)
    return {"vdw": without, "dispersion_correction": enabled - without}


def _add(components: dict[str, float], name: str, value: float) -> None:
    components[name] = components.get(name, 0.0) + value


def run_openmm_single_point(  # noqa: PLR0912, PLR0915
    *,
    system_file: Path,
    coordinate_file: Path,
    coordinate_sha256: str,
    classification: OpenMMForceClassification,
    platform_name: str = "Reference",
) -> EnergyReport:
    """Evaluate one serialized OpenMM system on explicit shared coordinates.

    The function never infers Urey--Bradley or improper membership from force
    classes.  Custom forces require an explicit role, and every unsupported
    nonzero force fails closed.
    """
    system_file = system_file.resolve()
    coordinate_file = coordinate_file.resolve()
    for path, description in ((system_file, "system"), (coordinate_file, "coordinate")):
        if not path.is_file():
            raise FileNotFoundError(f"OpenMM {description} file does not exist: {path}")
    actual_sha256 = hashlib.sha256(coordinate_file.read_bytes()).hexdigest()
    if actual_sha256 != coordinate_sha256:
        raise ValueError("OpenMM coordinate file does not match the shared checksum")

    openmm = importlib.import_module("openmm")
    app = importlib.import_module("openmm.app")
    unit = importlib.import_module("openmm.unit")
    system = openmm.XmlSerializer.deserialize(system_file.read_text(encoding="utf-8"))
    positions, box = _coordinate_data(openmm, app, unit, coordinate_file)
    if len(positions) != system.getNumParticles():
        raise ValueError("OpenMM coordinate and system particle counts differ")
    if system.usesPeriodicBoundaryConditions() and box is None:
        raise ValueError("periodic OpenMM systems require box vectors in the coordinate file")
    try:
        platform = openmm.Platform.getPlatformByName(platform_name)
    except Exception as error:
        raise ValueError(f"OpenMM platform is unavailable: {platform_name}") from error

    force_count = system.getNumForces()
    referenced = (
        set(classification.urey_bradley_bonds)
        | set(classification.periodic_impropers)
        | set(classification.rb_impropers)
        | set(classification.force_roles)
    )
    if any(index >= force_count for index in referenced):
        raise ValueError("OpenMM classification refers to a missing force index")
    interaction_maps = (
        set(classification.urey_bradley_bonds),
        set(classification.periodic_impropers),
        set(classification.rb_impropers),
    )
    if any(set(classification.force_roles) & indices for indices in interaction_maps):
        raise ValueError("an OpenMM force cannot have both a role and an interaction list")
    expected_maps = {
        "HarmonicBondForce": set(classification.urey_bradley_bonds),
        "PeriodicTorsionForce": set(classification.periodic_impropers),
        "RBTorsionForce": set(classification.rb_impropers),
    }
    for class_name, supplied in expected_maps.items():
        actual = {
            index
            for index in range(force_count)
            if system.getForce(index).__class__.__name__ == class_name
        }
        if supplied != actual:
            raise ValueError(
                f"OpenMM classification must explicitly list every {class_name} force "
                f"(expected {sorted(actual)}, received {sorted(supplied)})"
            )

    components: dict[str, float] = {}
    force_metadata: list[dict[str, object]] = []
    allowed_custom_roles = {
        "CustomNonbondedForce": {"vdw"},
        "CustomBondForce": {"vdw_14", "electrostatics_14"},
        "CustomTorsionForce": {"harmonic_improper", "periodic_improper"},
        "CustomExternalForce": {"restraint"},
        "CustomCVForce": {"restraint"},
        "CustomCentroidBondForce": {"restraint"},
        "CustomCompoundBondForce": {"restraint"},
    }
    for force_index in range(force_count):
        force = system.getForce(force_index)
        class_name = force.__class__.__name__
        role = classification.force_roles.get(force_index)
        isolated = _one_force_system(openmm, system, force_index)
        raw_energy = _context_energy(openmm, unit, isolated, positions, box, platform)
        force_metadata.append(
            {
                "index": force_index,
                "class": class_name,
                "name": force.getName(),
                "kj_mol": raw_energy,
            }
        )

        if role is not None:
            allowed = allowed_custom_roles.get(class_name)
            if allowed is None or role not in allowed:
                raise ValueError(
                    f"OpenMM force {force_index} ({class_name}) cannot use role {role!r}"
                )
            _add(components, role, raw_energy)
        elif class_name == "HarmonicBondForce":
            selected = set(classification.urey_bradley_bonds.get(force_index, ()))
            if any(index >= force.getNumBonds() for index in selected):
                raise ValueError("OpenMM Urey--Bradley classification refers to a missing bond")
            if selected:
                bonds = _one_force_system(openmm, system, force_index)
                _zero_harmonic_bonds(bonds.getForce(0), selected)
                bond_energy = _context_energy(openmm, unit, bonds, positions, box, platform)
                _add(components, "bond", bond_energy)
                _add(components, "urey_bradley", raw_energy - bond_energy)
            else:
                _add(components, "bond", raw_energy)
        elif class_name == "HarmonicAngleForce":
            _add(components, "angle", raw_energy)
        elif class_name == "PeriodicTorsionForce":
            selected = set(classification.periodic_impropers.get(force_index, ()))
            if any(index >= force.getNumTorsions() for index in selected):
                raise ValueError("OpenMM improper classification refers to a missing torsion")
            if selected:
                propers = _one_force_system(openmm, system, force_index)
                _zero_periodic_torsions(propers.getForce(0), selected)
                proper_energy = _context_energy(openmm, unit, propers, positions, box, platform)
                _add(components, "proper_dihedral", proper_energy)
                _add(components, "periodic_improper", raw_energy - proper_energy)
            else:
                _add(components, "proper_dihedral", raw_energy)
        elif class_name == "RBTorsionForce":
            selected = set(classification.rb_impropers[force_index])
            if any(index >= force.getNumTorsions() for index in selected):
                raise ValueError("OpenMM RB improper classification refers to a missing torsion")
            if selected:
                propers = _one_force_system(openmm, system, force_index)
                _zero_rb_torsions(propers.getForce(0), selected)
                proper_energy = _context_energy(openmm, unit, propers, positions, box, platform)
                _add(components, "proper_dihedral", proper_energy)
                _add(components, "periodic_improper", raw_energy - proper_energy)
            else:
                _add(components, "proper_dihedral", raw_energy)
        elif class_name == "NonbondedForce":
            for name, value in _nonbonded_components(
                openmm, unit, system, force_index, positions, box, platform
            ).items():
                _add(components, name, value)
        elif class_name == "CustomNonbondedForce":
            raise ValueError(f"OpenMM force {force_index} ({class_name}) requires an explicit role")
        elif class_name.startswith("Custom"):
            raise ValueError(f"OpenMM force {force_index} ({class_name}) requires an explicit role")
        elif class_name in {"CMMotionRemover", "MonteCarloBarostat"}:
            if raw_energy != 0:
                raise ValueError(f"OpenMM diagnostic force {force_index} has nonzero energy")
        elif raw_energy != 0:
            raise ValueError(
                f"OpenMM force {force_index} ({class_name}) has unsupported nonzero energy"
            )

    # Explicit CustomNonbondedForce roles need long-range correction isolation.
    for force_index, role in classification.force_roles.items():
        force = system.getForce(force_index)
        if force.__class__.__name__ == "CustomNonbondedForce":
            if role != "vdw":
                raise ValueError("CustomNonbondedForce is supported only with the vdw role")
            raw = next(item["kj_mol"] for item in force_metadata if item["index"] == force_index)
            components["vdw"] -= cast(float, raw)
            for name, value in _custom_nonbonded_components(
                openmm, unit, system, force_index, positions, box, platform
            ).items():
                _add(components, name, value)

    components["potential"] = _context_energy(openmm, unit, system, positions, box, platform)
    return openmm_energy_report(
        components,
        coordinate_sha256=coordinate_sha256,
        extra={
            "platform": platform_name,
            "forces": force_metadata,
            "classified_components_kj_mol": {
                name: value for name, value in components.items() if name != "potential"
            },
        },
    )


def openmm_energy_report(
    components: Mapping[str, EnergyValue],
    *,
    coordinate_sha256: str,
    extra: Mapping[str, object] | None = None,
) -> EnergyReport:
    """Return a canonical report from already classified OpenMM components.

    All nonbonded 1--4 contributions are folded into their corresponding
    canonical aggregate.  Urey--Bradley contributes to ``angle`` and both
    periodic and harmonic impropers contribute to ``improper``.
    """
    unknown = set(components) - OPENMM_COMPONENTS
    if unknown:
        raise ValueError(f"OpenMM components contain unsupported terms: {sorted(unknown)}")
    if "potential" not in components:
        raise ValueError("OpenMM components are missing potential")

    def total(*names: str) -> EnergyValue:
        return sum((components.get(name, 0) for name in names), start=0)

    terms: dict[str, EnergyValue] = {
        "bond": total("bond"),
        "angle": total("angle", "urey_bradley"),
        "dihedral": total("proper_dihedral"),
        "improper": total("periodic_improper", "harmonic_improper"),
        "vdw": total("vdw", "vdw_14"),
        "electrostatics": total("electrostatics", "electrostatics_14"),
        "restraint": total("restraint"),
        "dispersion_correction": total("dispersion_correction"),
        "potential": components["potential"],
    }
    return EnergyReport(
        engine="openmm",
        coordinate_sha256=coordinate_sha256,
        terms=terms,
        extra={} if extra is None else extra,
    )
