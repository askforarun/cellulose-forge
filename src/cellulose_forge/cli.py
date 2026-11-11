"""Command-line interface for CelluloseForge."""

from __future__ import annotations

import argparse
import json
import tempfile
import warnings
from collections.abc import Sequence
from pathlib import Path

from cellulose_forge import __version__
from cellulose_forge.amber import write_glycam_amber_system
from cellulose_forge.charmm36 import (
    FiniteCharmm36Topology,
    PeriodicCharmm36Topology,
    derive_periodic_charmm36_topology,
    finite_charmm36_interaction_counts,
    map_finite_charmm36,
    map_periodic_charmm36,
    native_finite_charmm36_topology,
    native_periodic_charmm36_topology,
    periodic_charmm36_interaction_counts,
    write_charmm36_system,
)
from cellulose_forge.clashes import (
    validate_finite_periodic_image_contact,
    validate_periodic_glycosidic_bonds,
    validate_periodic_gromacs_topology_bonds,
)
from cellulose_forge.coordinates import write_coordinates
from cellulose_forge.crystallography import AllomorphDefinition, MillerIndex, load_allomorph
from cellulose_forge.doctor import render_json, render_text, run_diagnostics
from cellulose_forge.forcefields import (
    CHARMM36_RELEASE,
    discover_charmm36_directory,
    install_charmm36,
)
from cellulose_forge.geometry import (
    I_BETA_FIBRIL_PROFILES,
    BulkCrystal,
    build_36_chain_fibril,
    build_bulk_crystal,
    build_i_beta_profile_fibril,
    build_slab_with_sizing,
    build_unit_cell,
    validate_cellulose_chirality,
)
from cellulose_forge.glycam import GLYCAM_LEAPRC, generate_glycam_gromacs_reference
from cellulose_forge.glycam_mapping import (
    GlycamMappedSystem,
    map_crystal_to_finite_glycam,
    map_crystal_to_periodic_glycam,
)
from cellulose_forge.gromacs import write_glycam_gromacs_system
from cellulose_forge.openmm import (
    OpenMMSystemFiles,
    heavy_atom_restraints,
    write_openmm_system,
)
from cellulose_forge.restraints import (
    ChainRestraints,
    group_mapped_chains,
    resolve_chain_restraints,
)

_MINIMUM_FINITE_LONGITUDINAL_PADDING_ANGSTROM = 3.0
_RECOMMENDED_FINITE_LONGITUDINAL_PADDING_ANGSTROM = 12.0


def _parse_fibril_profile(value: str) -> tuple[int, ...]:
    """Parse a comma-separated, bottom-to-top fibril row profile."""
    try:
        rows = tuple(int(item.strip()) for item in value.split(","))
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "fibril profiles must be comma-separated positive integers"
        ) from error
    if not rows or any(row < 1 for row in rows):
        raise argparse.ArgumentTypeError(
            "fibril profiles must be comma-separated positive integers"
        )
    return rows


def build_parser() -> argparse.ArgumentParser:
    """Create the top-level command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cellulose-forge",
        description="Build crystalline cellulose structures for GROMACS simulations.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    subparsers = parser.add_subparsers(dest="command")
    build = subparsers.add_parser("build", help="build a cellulose structure")
    build.add_argument(
        "--allomorph",
        choices=["I-alpha", "I-beta", "II", "III-I"],
        required=True,
    )
    build.add_argument(
        "--shape",
        choices=["unit-cell", "bulk", "crystal", "slab", "fibril"],
        required=True,
    )
    build.add_argument(
        "--cells",
        nargs=3,
        type=int,
        metavar=("NA", "NB", "NC"),
        help="unit-cell repetitions along crystallographic a, b, and c",
    )
    build.add_argument(
        "--face",
        nargs=3,
        type=int,
        metavar=("H", "K", "L"),
        help="signed Miller index of the exposed slab face",
    )
    build.add_argument("--layers", type=int, help="number of complete molecular slab layers")
    build.add_argument("--chains-per-layer", type=int, help="cellulose chains in each slab layer")
    build.add_argument(
        "--degree-polymerization",
        type=int,
        help="anhydroglucose units (AGUs) in each cellulose chain",
    )
    build.add_argument(
        "--fibril-model",
        choices=["18-chain", "24-chain", "36-chain"],
        help="named fibril cross-section",
    )
    build.add_argument(
        "--fibril-profile",
        type=_parse_fibril_profile,
        help=(
            "custom I-beta row counts from low to high Cartesian z, "
            "for example 2,3,4,4,3,2"
        ),
    )
    build.add_argument(
        "--chain-boundary",
        choices=["finite", "periodic"],
        help="finite termini or a covalent beta(1->4) bond across the chain boundary",
    )
    build.add_argument("--slab-width", type=float, help="minimum slab width in angstrom")
    build.add_argument("--slab-thickness", type=float, help="minimum slab thickness in angstrom")
    build.add_argument(
        "--slab-length", type=float, help="minimum chain-direction length in angstrom"
    )
    build.add_argument(
        "--vacuum",
        type=float,
        default=0.0,
        help="total slab-normal or fibril-transverse padding in angstrom (default: 0)",
    )
    build.add_argument(
        "--padding-x",
        type=float,
        help="total fibril padding along Cartesian x; overrides --vacuum for x",
    )
    build.add_argument(
        "--padding-y",
        type=float,
        help=(
            "total fibril padding along the second transverse direction "
            "(Cartesian z); overrides --vacuum"
        ),
    )
    build.add_argument(
        "--padding-longitudinal",
        type=float,
        help=(
            "total chain-direction padding for finite slabs or fibrils; "
            "periodic chains require zero"
        ),
    )
    build.add_argument(
        "--padding-z",
        type=float,
        help="deprecated alias for --padding-longitudinal",
    )
    build.add_argument(
        "--force-field",
        choices=["none", "glycam", "charmm36"],
        default="none",
        help="write coordinates only or a simulation-ready GROMACS system",
    )
    build.add_argument(
        "--topology-format",
        choices=["gromacs", "amber", "both"],
        help="deprecated alias for --md-engine on GLYCAM builds",
    )
    build.add_argument(
        "--md-engine",
        choices=["gromacs", "amber", "openmm"],
        action="append",
        help=(
            "target MD engine; repeat for multiple outputs (default: gromacs for "
            "force-field builds)"
        ),
    )
    build.add_argument(
        "--charmm-template",
        type=Path,
        help="externally supplied periodic CHARMM36 cellulose-chain ITP",
    )
    build.add_argument(
        "--charmm-force-field-directory",
        type=Path,
        help="CHARMM36 GROMACS force-field directory containing forcefield.itp",
    )
    build.add_argument(
        "--tleap-executable",
        default="tleap",
        help="tleap executable used for GLYCAM reference generation",
    )
    build.add_argument(
        "--parmed-python",
        type=Path,
        help="optional Python interpreter containing ParmEd",
    )
    build.add_argument(
        "--openmm-python",
        type=Path,
        help="optional Python interpreter containing OpenMM",
    )
    restraint_selector = build.add_mutually_exclusive_group()
    restraint_selector.add_argument(
        "--restrain-chains",
        nargs="+",
        type=int,
        metavar="CHAIN_ID",
        help="sequential CelluloseForge chain IDs to restrain",
    )
    restraint_selector.add_argument(
        "--restrain-layers",
        nargs="+",
        type=int,
        metavar="LAYER_ID",
        help=("one-based transverse layers defined by the selected slab or fibril geometry"),
    )
    build.add_argument(
        "--restraint-force",
        nargs=3,
        type=float,
        metavar=("FX", "FY", "FZ"),
        help="heavy-atom position-restraint force constants in kJ mol^-1 nm^-2",
    )
    build.add_argument("--output", type=Path, required=True)

    doctor = subparsers.add_parser(
        "doctor",
        help="check the environment for coordinate-only, GLYCAM, and GROMACS readiness",
    )
    doctor.add_argument(
        "--tleap-executable",
        default="tleap",
        help="tleap executable to check for GLYCAM reference generation",
    )
    doctor.add_argument(
        "--parmed-python",
        type=Path,
        help="Python interpreter to check for an importable ParmEd",
    )
    doctor.add_argument(
        "--openmm-python",
        type=Path,
        help="Python interpreter to check for an importable OpenMM",
    )
    doctor.add_argument("--json", action="store_true", help="emit machine-readable JSON")

    forcefield = subparsers.add_parser(
        "forcefield", help="install or locate supported force-field data"
    )
    forcefield_subparsers = forcefield.add_subparsers(dest="forcefield_command", required=True)
    forcefield_install = forcefield_subparsers.add_parser(
        "install", help="download and install a supported force-field release"
    )
    forcefield_install.add_argument("name", choices=["charmm36"])
    forcefield_install.add_argument(
        "--archive",
        type=Path,
        help="install from a local official archive instead of downloading it",
    )
    forcefield_path = forcefield_subparsers.add_parser(
        "path", help="print the installed path of a supported force field"
    )
    forcefield_path.add_argument("name", choices=["charmm36"])
    return parser


def main(argv: Sequence[str] | None = None) -> int:  # noqa: PLR0912
    """Run the command-line interface."""
    args = build_parser().parse_args(argv)
    if args.command == "forcefield":
        if args.forcefield_command == "install":
            installed = install_charmm36(archive_path=args.archive)
            print(f"Installed {CHARMM36_RELEASE}: {installed}")
        else:
            print(discover_charmm36_directory())
        return 0
    if args.command == "doctor":
        checks = run_diagnostics(
            tleap_executable=args.tleap_executable,
            parmed_python=str(args.parmed_python) if args.parmed_python else None,
            openmm_python=str(args.openmm_python) if args.openmm_python else None,
        )
        print(render_json(checks) if args.json else render_text(checks))
        optional_engines = {"gmx", "openmm"}
        return 0 if all(check.ok for check in checks if check.name not in optional_engines) else 1
    if args.command == "build":
        _resolve_md_engines(args, build_parser())
        if args.padding_longitudinal is not None and args.padding_z is not None:
            build_parser().error(
                "--padding-longitudinal and deprecated --padding-z are mutually exclusive"
            )
        if args.padding_z is not None:
            warnings.warn(
                "--padding-z is deprecated; use --padding-longitudinal",
                FutureWarning,
                stacklevel=1,
            )
        args.padding_longitudinal = (
            args.padding_z
            if args.padding_longitudinal is None
            else args.padding_longitudinal
        )
        allomorph_key = {
            "I-alpha": "i_alpha",
            "I-beta": "i_beta",
            "II": "ii",
            "III-I": "iii_i",
        }[args.allomorph]
        definition = load_allomorph(allomorph_key)
        if args.shape != "fibril" and any(
            value is not None for value in (args.padding_x, args.padding_y)
        ):
            build_parser().error("--padding-x and --padding-y are valid only with --shape fibril")
        if args.shape not in {"slab", "fibril"} and args.padding_longitudinal is not None:
            build_parser().error(
                "--padding-longitudinal is valid only with --shape slab or fibril"
            )
        _validate_longitudinal_padding(args, build_parser())
        if args.shape == "unit-cell":
            if args.cells is not None:
                build_parser().error("--cells is not valid with --shape unit-cell")
            crystal = build_unit_cell(definition)
        elif args.shape in {"bulk", "crystal"}:
            if args.cells is None:
                build_parser().error("--cells NA NB NC is required for bulk crystals")
            crystal = build_bulk_crystal(definition, tuple(args.cells))
        elif args.shape == "slab":
            crystal = _build_slab(args, definition, build_parser())
        else:
            crystal = _build_fibril(args, definition, build_parser())
        if args.force_field == "glycam":
            _write_glycam_build(args, crystal, build_parser())
        elif args.force_field == "charmm36":
            _write_charmm36_build(args, crystal, build_parser())
        else:
            if any(
                value is not None
                for value in (
                    args.restrain_chains,
                    args.restrain_layers,
                    args.restraint_force,
                )
            ):
                build_parser().error("chain restraints require --force-field glycam or charmm36")
            args.output.parent.mkdir(parents=True, exist_ok=True)
            written_path = write_coordinates(crystal, args.output)
            _write_coordinate_only_report(args, crystal, written_path)
            if written_path != args.output:
                print(
                    f"Requested structure has {crystal.atom_count} atoms; "
                    f"wrote mmCIF instead of PDB: {written_path}"
                )
    return 0


def _resolve_md_engines(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> None:
    """Resolve the engine selector and its deprecated topology-format alias."""
    if args.md_engine is not None and args.topology_format is not None:
        parser.error("--md-engine and deprecated --topology-format are mutually exclusive")
    if args.topology_format is not None:
        warnings.warn(
            "--topology-format is deprecated; use one or more --md-engine options",
            FutureWarning,
            stacklevel=1,
        )
        legacy_engines = {
            "gromacs": ("gromacs",),
            "amber": ("amber",),
            "both": ("gromacs", "amber"),
        }
        engines = legacy_engines[args.topology_format]
    elif args.md_engine is not None:
        engines = tuple(dict.fromkeys(args.md_engine))
    elif args.force_field == "none":
        engines = ()
    else:
        engines = ("gromacs",)

    if args.force_field == "none" and engines:
        parser.error("--md-engine requires --force-field glycam or charmm36")
    if args.force_field == "charmm36" and not set(engines) <= {"gromacs", "openmm"}:
        parser.error("CHARMM36 supports --md-engine gromacs and openmm")
    if args.force_field == "glycam" and not engines:
        parser.error("GLYCAM builds require at least one --md-engine")
    args.md_engines = engines


def _write_coordinate_only_report(
    args: argparse.Namespace,
    crystal: BulkCrystal,
    coordinate_path: Path,
) -> Path:
    """Write a sibling build report for a coordinate-only construction."""
    requested_output = Path(args.output)
    report_path = requested_output.with_suffix(".build_report.json")
    report = {
        "cellulose_forge_version": __version__,
        "allomorph": crystal.allomorph,
        "shape": args.shape,
        "force_field": "none",
        "md_engines": list(args.md_engines),
        "chain_boundary": args.chain_boundary,
        "degree_polymerization": _crystal_degree_polymerization(crystal),
        "chain_count": len(set(crystal.chain_numbers)),
        "crystallographic_atom_count": crystal.atom_count,
        "requested_output": str(requested_output),
        "outputs": {"coordinates": coordinate_path.name},
        "geometry_annotations": list(crystal.annotations),
        "chirality_validation": validate_cellulose_chirality(crystal),
    }
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report_path


def _validate_longitudinal_padding(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
) -> None:
    """Reject overlapping finite force-field boxes and warn below the safe baseline."""
    if args.shape not in {"slab", "fibril"} or args.chain_boundary != "finite":
        return
    padding = 0.0 if args.padding_longitudinal is None else args.padding_longitudinal
    if (
        args.force_field != "none"
        and padding < _MINIMUM_FINITE_LONGITUDINAL_PADDING_ANGSTROM
    ):
        parser.error(
            "finite force-field chains require at least 3 angstrom longitudinal "
            "padding; use --padding-longitudinal 12 or greater"
        )
    if padding < _RECOMMENDED_FINITE_LONGITUDINAL_PADDING_ANGSTROM:
        qualifier = "coordinate-only " if args.force_field == "none" else ""
        warnings.warn(
            f"{qualifier}finite chains have only {padding:g} angstrom longitudinal "
            "padding; use at least 12 angstrom before periodic simulation so chain "
            "ends do not interact with their own images",
            UserWarning,
            stacklevel=1,
        )


def _validate_written_periodic_topology(
    mapped: GlycamMappedSystem,
    molecule_topology_path: Path,
    system_topology_path: Path,
    parser: argparse.ArgumentParser,
) -> dict[str, float | int]:
    """Validate the bonds actually emitted in a GROMACS molecule topology."""
    try:
        return validate_periodic_gromacs_topology_bonds(
            mapped, molecule_topology_path, system_topology_path
        )
    except ValueError as error:
        parser.error(str(error))


def _write_openmm_output(  # noqa: PLR0913
    args: argparse.Namespace,
    mapped: GlycamMappedSystem,
    restraints: ChainRestraints | None,
    *,
    gromacs_topology: Path,
    gromacs_coordinates: Path,
    parser: argparse.ArgumentParser,
) -> tuple[OpenMMSystemFiles, dict[str, object]]:
    """Serialize and validate the OpenMM system for a force-field build."""
    try:
        return write_openmm_system(
            mapped,
            gromacs_topology,
            gromacs_coordinates,
            args.output,
            periodic_chains=args.chain_boundary == "periodic",
            restraints=(
                None
                if restraints is None
                else heavy_atom_restraints(
                    mapped, restraints.chain_ids, restraints.force_constants
                )
            ),
            python_executable=args.openmm_python,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as error:
        parser.error(str(error))


def _engine_outputs(
    coordinate_path: Path,
    engine_files: dict[str, str],
    openmm_files: OpenMMSystemFiles | None,
) -> dict[str, str]:
    """Collect the written file names for the build report."""
    outputs = {"coordinates": coordinate_path.name, **engine_files}
    if openmm_files is not None:
        outputs.update(
            {
                "openmm_system": openmm_files.system.name,
                "openmm_state": openmm_files.state.name,
                "openmm_topology": openmm_files.topology.name,
                "openmm_runner": openmm_files.runner.name,
            }
        )
    return outputs


def _write_charmm36_build(  # noqa: PLR0912, PLR0915
    args: argparse.Namespace,
    crystal: BulkCrystal,
    parser: argparse.ArgumentParser,
) -> None:
    """Map cellulose to a native CHARMM36 topology."""
    if args.output.suffix:
        parser.error("--output must be a directory when --force-field charmm36 is selected")
    degree_polymerization = _crystal_degree_polymerization(crystal)
    force_field_directory = (
        args.charmm_force_field_directory
        if args.charmm_force_field_directory is not None
        else discover_charmm36_directory()
    )
    topology: PeriodicCharmm36Topology | FiniteCharmm36Topology
    periodic_bond_validation = None
    written_periodic_bond_validation = None
    if args.chain_boundary == "periodic":
        topology = (
            derive_periodic_charmm36_topology(
                args.charmm_template,
                degree_polymerization=degree_polymerization,
            )
            if args.charmm_template is not None
            else native_periodic_charmm36_topology(
                force_field_directory,
                degree_polymerization=degree_polymerization,
            )
        )
        mapped = map_periodic_charmm36(crystal, topology)
        try:
            periodic_bond_validation = validate_periodic_glycosidic_bonds(mapped)
        except ValueError as error:
            parser.error(str(error))
        interaction_counts = periodic_charmm36_interaction_counts(topology)
    else:
        if args.charmm_template is not None:
            parser.error("--charmm-template is only supported for periodic chains")
        topology = native_finite_charmm36_topology(
            force_field_directory,
            degree_polymerization=degree_polymerization,
        )
        mapped = map_finite_charmm36(crystal, topology)
        interaction_counts = finite_charmm36_interaction_counts(topology)
    periodic_image_contact = None
    if args.chain_boundary == "finite":
        try:
            periodic_image_contact = validate_finite_periodic_image_contact(mapped)
        except ValueError as error:
            parser.error(str(error))
    restraints = _resolve_restraints(args, crystal, parser)
    chain_mapping = None
    if restraints is not None:
        mapped, chain_mapping = group_mapped_chains(
            mapped,
            restraints.chain_ids,
            free_molecule_type=topology.molecule_name,
            restrained_molecule_type=f"{topology.molecule_name}_RESTRAINED",
        )
    args.output.mkdir(parents=True, exist_ok=True)
    coordinate_path = write_coordinates(crystal, args.output / "structure.pdb")
    write_gromacs = "gromacs" in args.md_engines
    openmm_files: OpenMMSystemFiles | None = None
    openmm_metadata: dict[str, object] | None = None
    with tempfile.TemporaryDirectory(prefix="cellulose_forge_charmm36_") as temporary:
        gromacs_directory = args.output if write_gromacs else Path(temporary)
        gro_path, itp_path, topology_path, force_field_path = write_charmm36_system(
            mapped,
            topology,
            force_field_directory,
            gromacs_directory,
            restrained_chain_count=len(restraints.chain_ids) if restraints else 0,
            restraint_force=restraints.force_constants if restraints else None,
        )
        if args.chain_boundary == "periodic":
            written_periodic_bond_validation = _validate_written_periodic_topology(
                mapped, itp_path, topology_path, parser
            )
        if "openmm" in args.md_engines:
            openmm_files, openmm_metadata = _write_openmm_output(
                args,
                mapped,
                restraints,
                gromacs_topology=topology_path,
                gromacs_coordinates=gro_path,
                parser=parser,
            )
    report = {
        "cellulose_forge_version": __version__,
        "allomorph": crystal.allomorph,
        "shape": args.shape,
        "force_field": "CHARMM36",
        "md_engines": list(args.md_engines),
        "status": "supported",
        "chain_boundary": args.chain_boundary,
        "degree_polymerization": topology.degree_polymerization,
        "chain_count": mapped.chain_count,
        "crystallographic_atom_count": crystal.atom_count,
        "force_field_atom_count": mapped.atom_count,
        "generated_chain_charge": topology.total_charge,
        "interaction_counts_per_chain": interaction_counts,
        "topology_source": (
            str(args.charmm_template)
            if args.charmm_template is not None
            else f"native BGLC residue with {args.chain_boundary} 14bb linkage"
        ),
        "reference_force_field": str(force_field_directory),
        "charmm36_release": CHARMM36_RELEASE,
        "outputs": _engine_outputs(
            coordinate_path,
            {
                "gromacs_coordinates": gro_path.name,
                "gromacs_molecule_topology": itp_path.name,
                "gromacs_topology": topology_path.name,
                "gromacs_force_field": force_field_path.name,
            }
            if write_gromacs
            else {},
            openmm_files,
        ),
        "geometry_annotations": list(crystal.annotations),
        "chirality_validation": validate_cellulose_chirality(crystal),
    }
    if openmm_metadata is not None:
        report["openmm_validation"] = openmm_metadata
    if periodic_image_contact is not None:
        report["closest_longitudinal_periodic_image_contact"] = (
            periodic_image_contact.as_dict()
        )
    if periodic_bond_validation is not None:
        report["periodic_glycosidic_bond_validation"] = {
            "mapped_geometry": periodic_bond_validation,
            "written_topology": written_periodic_bond_validation,
        }
    if restraints is not None:
        report["position_restraints"] = {
            "define": "POSRES_CELLULOSE",
            "selected_chain_ids": list(restraints.chain_ids),
            "selected_layer_ids": list(restraints.layer_ids),
            "layer_mapping": [
                {"layer_id": layer_id, "chain_ids": list(chain_ids)}
                for layer_id, chain_ids in restraints.layer_mapping
            ],
            "force_constants_kj_mol_nm2": list(restraints.force_constants),
            "atom_selection": "heavy atoms in each selected chain",
            "chain_mapping": chain_mapping,
        }
    if crystal.allomorph == "II":
        report["generated_force_field_coordinates"] = {
            "atoms": ["HO2", "HO3", "HO6"],
            "scope": "each anhydroglucose residue",
            "method": "I-beta reference O-H geometry aligned on C1/C3/C5",
        }
    (args.output / "build_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_glycam_build(  # noqa: PLR0912, PLR0915
    args: argparse.Namespace,
    crystal: BulkCrystal,
    parser: argparse.ArgumentParser,
) -> None:
    """Generate a GLYCAM reference, map coordinates, and write GROMACS files."""
    if args.output.suffix:
        parser.error("--output must be a directory when --force-field glycam is selected")
    chain_boundary = args.chain_boundary
    if args.shape in {"bulk", "crystal", "unit-cell"}:
        if chain_boundary is None:
            parser.error(
                "--chain-boundary periodic is required for GLYCAM bulk and unit-cell builds"
            )
        if chain_boundary != "periodic":
            parser.error("GLYCAM bulk and unit-cell builds require --chain-boundary periodic")
    elif chain_boundary is None:
        parser.error("--chain-boundary finite|periodic is required with --force-field glycam")
    assert chain_boundary is not None
    degree_polymerization = _crystal_degree_polymerization(crystal)
    output_directory = args.output
    write_gromacs = "gromacs" in args.md_engines
    write_amber = "amber" in args.md_engines
    write_openmm = "openmm" in args.md_engines
    if write_amber and chain_boundary == "periodic":
        parser.error(
            "native Amber output currently supports finite chains only; Amber does not "
            "minimum-image bonded terms across the periodic chain boundary. Use "
            "--chain-boundary finite with longitudinal padding, or use GROMACS or "
            "OpenMM output for covalently periodic chains"
        )
    gro_path: Path | None = None
    topology_path: Path | None = None
    amber_topology_path: Path | None = None
    amber_coordinate_path: Path | None = None
    amber_metadata: dict[str, object] | None = None
    openmm_files: OpenMMSystemFiles | None = None
    openmm_metadata: dict[str, object] | None = None

    with tempfile.TemporaryDirectory(prefix="cellulose_forge_glycam_") as temporary:
        files, reference_metadata = generate_glycam_gromacs_reference(
            degree_polymerization=degree_polymerization,
            chain_boundary=chain_boundary,
            directory=temporary,
            tleap_executable=args.tleap_executable,
            parmed_python=args.parmed_python,
        )
        if chain_boundary == "periodic":
            mapped = map_crystal_to_periodic_glycam(crystal, files.pdb)
        else:
            mapped = map_crystal_to_finite_glycam(crystal, files.pdb)
        periodic_image_contact = None
        periodic_bond_validation = None
        written_periodic_bond_validation = None
        if chain_boundary == "periodic":
            try:
                periodic_bond_validation = validate_periodic_glycosidic_bonds(mapped)
            except ValueError as error:
                parser.error(str(error))
        if chain_boundary == "finite":
            try:
                periodic_image_contact = validate_finite_periodic_image_contact(mapped)
            except ValueError as error:
                parser.error(str(error))
        restraints = _resolve_restraints(args, crystal, parser)
        chain_mapping = None
        if restraints is not None:
            mapped, chain_mapping = group_mapped_chains(mapped, restraints.chain_ids)
        output_directory.mkdir(parents=True, exist_ok=True)
        coordinate_path = write_coordinates(crystal, output_directory / "structure.pdb")
        if write_gromacs or write_openmm:
            gromacs_directory = (
                output_directory if write_gromacs else Path(temporary) / "gromacs"
            )
            written_gro, written_topology = write_glycam_gromacs_system(
                mapped,
                files.gromacs_topology,
                gromacs_directory,
                restrained_chain_count=len(restraints.chain_ids) if restraints else 0,
                restraint_force=restraints.force_constants if restraints else None,
            )
            if write_gromacs:
                gro_path, topology_path = written_gro, written_topology
            if write_openmm:
                openmm_files, openmm_metadata = _write_openmm_output(
                    args,
                    mapped,
                    restraints,
                    gromacs_topology=written_topology,
                    gromacs_coordinates=written_gro,
                    parser=parser,
                )
        if write_amber:
            amber_topology_path, amber_coordinate_path, amber_metadata = (
                write_glycam_amber_system(
                    mapped,
                    files.prmtop,
                    files.inpcrd,
                    output_directory,
                    python_executable=args.parmed_python,
                )
            )
        if chain_boundary == "periodic" and topology_path is not None:
            written_periodic_bond_validation = _validate_written_periodic_topology(
                mapped, topology_path, topology_path, parser
            )

    engine_files: dict[str, str] = {}
    if gro_path is not None and topology_path is not None:
        engine_files.update(
            {
                "gromacs_coordinates": gro_path.name,
                "gromacs_topology": topology_path.name,
            }
        )
    if amber_topology_path is not None and amber_coordinate_path is not None:
        engine_files.update(
            {
                "amber_topology": amber_topology_path.name,
                "amber_coordinates": amber_coordinate_path.name,
            }
        )
    outputs = _engine_outputs(coordinate_path, engine_files, openmm_files)
    report = {
        "cellulose_forge_version": __version__,
        "allomorph": crystal.allomorph,
        "shape": args.shape,
        "force_field": "GLYCAM",
        "force_field_release": GLYCAM_LEAPRC,
        "chain_boundary": chain_boundary,
        "degree_polymerization": degree_polymerization,
        "chain_count": mapped.chain_count,
        "crystallographic_atom_count": crystal.atom_count,
        "force_field_atom_count": mapped.atom_count,
        "md_engines": list(args.md_engines),
        "outputs": outputs,
        "geometry_annotations": list(crystal.annotations),
        "chirality_validation": validate_cellulose_chirality(crystal),
        "reference_generation": reference_metadata,
    }
    if amber_metadata is not None:
        report["amber_validation"] = amber_metadata
    if openmm_metadata is not None:
        report["openmm_validation"] = openmm_metadata
    if periodic_image_contact is not None:
        report["closest_longitudinal_periodic_image_contact"] = (
            periodic_image_contact.as_dict()
        )
    if periodic_bond_validation is not None:
        report["periodic_glycosidic_bond_validation"] = {
            "mapped_geometry": periodic_bond_validation,
            "written_topology": written_periodic_bond_validation,
        }
    if restraints is not None:
        report["position_restraints"] = {
            "define": "POSRES_CELLULOSE",
            "selected_chain_ids": list(restraints.chain_ids),
            "selected_layer_ids": list(restraints.layer_ids),
            "layer_mapping": [
                {"layer_id": layer_id, "chain_ids": list(chain_ids)}
                for layer_id, chain_ids in restraints.layer_mapping
            ],
            "force_constants_kj_mol_nm2": list(restraints.force_constants),
            "atom_selection": "heavy atoms in each selected chain",
            "chain_mapping": chain_mapping,
        }
    if crystal.allomorph == "II":
        report["generated_force_field_coordinates"] = {
            "atoms": ["HO2", "HO3", "HO6"],
            "scope": "each anhydroglucose residue",
            "method": (
                "GLYCAM hydroxyl direction from C1/C3/C5 rigid alignment with "
                "the source O-H bond length preserved"
            ),
            "reason": "hydroxyl-hydrogen coordinates are absent from the diffraction dataset",
        }
    (output_directory / "build_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _crystal_degree_polymerization(crystal: BulkCrystal) -> int:
    """Return the common residue count for every generated cellulose chain."""
    degrees = {
        chain_number: len(
            {
                residue_number
                for atom_chain, residue_number in zip(
                    crystal.chain_numbers,
                    crystal.residue_numbers,
                    strict=True,
                )
                if atom_chain == chain_number
            }
        )
        for chain_number in dict.fromkeys(crystal.chain_numbers)
    }
    if not degrees or len(set(degrees.values())) != 1:
        msg = "output requires equal, nonempty cellulose chain lengths"
        raise ValueError(msg)
    return next(iter(degrees.values()))


def _resolve_restraints(
    args: argparse.Namespace,
    crystal: BulkCrystal,
    parser: argparse.ArgumentParser,
) -> ChainRestraints | None:
    """Resolve CLI restraint selectors, reporting validation as CLI errors."""
    try:
        return resolve_chain_restraints(
            crystal,
            chain_ids=args.restrain_chains,
            layer_ids=args.restrain_layers,
            force_constants=args.restraint_force,
        )
    except ValueError as error:
        parser.error(str(error))


def _build_fibril(
    args: argparse.Namespace,
    definition: AllomorphDefinition,
    parser: argparse.ArgumentParser,
) -> BulkCrystal:
    """Validate named-fibril options and build its crystallographic geometry."""
    if args.cells is not None or args.face is not None:
        parser.error("--cells and --face are not valid with --shape fibril")
    if any(
        value is not None
        for value in (
            args.layers,
            args.chains_per_layer,
            args.slab_width,
            args.slab_thickness,
            args.slab_length,
        )
    ):
        parser.error("slab sizing options are not valid with --shape fibril")
    if args.fibril_model is not None and args.fibril_profile is not None:
        parser.error("--fibril-model and --fibril-profile are mutually exclusive")
    if args.fibril_model is None and args.fibril_profile is None:
        parser.error("--fibril-model or --fibril-profile is required for fibrils")
    if args.degree_polymerization is None:
        parser.error("--degree-polymerization N is required for fibrils")
    if args.chain_boundary is None:
        parser.error("--chain-boundary finite|periodic is required for fibrils")
    if args.chain_boundary == "periodic" and args.padding_longitudinal not in {None, 0}:
        parser.error(
            "--padding-longitudinal must be zero with --chain-boundary periodic"
        )
    if args.fibril_model == "36-chain":
        return build_36_chain_fibril(
            definition,
            degree_polymerization=args.degree_polymerization,
            vacuum=args.vacuum,
            padding_x=args.padding_x,
            padding_y=args.padding_y,
            padding_z=args.padding_longitudinal,
            chain_boundary=args.chain_boundary,
        )
    model_name = args.fibril_model or "custom"
    row_counts = (
        I_BETA_FIBRIL_PROFILES[model_name]
        if args.fibril_model is not None
        else args.fibril_profile
    )
    assert row_counts is not None
    try:
        return build_i_beta_profile_fibril(
            definition,
            row_counts=row_counts,
            model_name=model_name,
            degree_polymerization=args.degree_polymerization,
            vacuum=args.vacuum,
            padding_x=args.padding_x,
            padding_y=args.padding_y,
            padding_z=args.padding_longitudinal,
            chain_boundary=args.chain_boundary,
        )
    except ValueError as error:
        parser.error(str(error))


def _build_slab(
    args: argparse.Namespace,
    definition: AllomorphDefinition,
    parser: argparse.ArgumentParser,
) -> BulkCrystal:
    """Resolve independent discrete or physical sizing for each slab axis."""
    if args.cells is not None:
        parser.error("--cells is not valid with --shape slab")
    if args.face is None:
        parser.error("--face H K L is required for slabs")
    if args.chain_boundary is None:
        parser.error("--chain-boundary finite|periodic is required for slabs")
    if args.chain_boundary == "periodic" and args.padding_longitudinal not in {None, 0}:
        parser.error(
            "--padding-longitudinal must be zero with --chain-boundary periodic"
        )
    padding_z = (
        0.0 if args.padding_longitudinal is None else args.padding_longitudinal
    )
    axis_options = (
        ("--chains-per-layer", args.chains_per_layer, "--slab-width", args.slab_width),
        ("--layers", args.layers, "--slab-thickness", args.slab_thickness),
        (
            "--degree-polymerization",
            args.degree_polymerization,
            "--slab-length",
            args.slab_length,
        ),
    )
    for count_name, count, size_name, size in axis_options:
        if count is not None and size is not None:
            parser.error(f"{count_name} and {size_name} are mutually exclusive")
        if count is None and size is None:
            parser.error(f"supply one of {count_name} or {size_name}")
    try:
        return build_slab_with_sizing(
            definition,
            MillerIndex(*args.face),
            layers=args.layers,
            chains_per_layer=args.chains_per_layer,
            degree_polymerization=args.degree_polymerization,
            width=args.slab_width,
            thickness=args.slab_thickness,
            length=args.slab_length,
            vacuum=args.vacuum,
            padding_z=padding_z,
            chain_boundary=args.chain_boundary,
        )
    except ValueError as error:
        parser.error(str(error))


if __name__ == "__main__":
    raise SystemExit(main())
