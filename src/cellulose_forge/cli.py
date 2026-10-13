"""Command-line interface for CelluloseForge."""

from __future__ import annotations

import argparse
import json
import tempfile
from collections.abc import Sequence
from pathlib import Path

from cellulose_forge import __version__
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
from cellulose_forge.coordinates import write_coordinates
from cellulose_forge.crystallography import AllomorphDefinition, MillerIndex, load_allomorph
from cellulose_forge.doctor import render_json, render_text, run_diagnostics
from cellulose_forge.forcefields import (
    CHARMM36_RELEASE,
    discover_charmm36_directory,
    install_charmm36,
)
from cellulose_forge.geometry import (
    BulkCrystal,
    build_36_chain_fibril,
    build_bulk_crystal,
    build_chain_parallel_slab,
    build_slab_by_dimensions,
    build_unit_cell,
)
from cellulose_forge.glycam import GLYCAM_LEAPRC, generate_glycam_gromacs_reference
from cellulose_forge.glycam_mapping import (
    map_crystal_to_finite_glycam,
    map_crystal_to_periodic_glycam,
)
from cellulose_forge.gromacs import write_glycam_gromacs_system
from cellulose_forge.restraints import (
    ChainRestraints,
    group_mapped_chains,
    resolve_chain_restraints,
)


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
        choices=["36-chain"],
        help="named literature-defined fibril cross-section",
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
        help="total fibril padding along Cartesian y; overrides --vacuum for y",
    )
    build.add_argument(
        "--padding-z",
        type=float,
        help="total chain-direction padding for finite fibrils; periodic fibrils require zero",
    )
    build.add_argument(
        "--force-field",
        choices=["none", "glycam", "charmm36"],
        default="none",
        help="write coordinates only or a simulation-ready GROMACS system",
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
        "--restrain-chains",
        nargs="+",
        type=int,
        metavar="CHAIN_ID",
        help="sequential CelluloseForge chain IDs to restrain",
    )
    build.add_argument(
        "--restraint-force",
        nargs=3,
        type=float,
        metavar=("FX", "FY", "FZ"),
        help="position-restraint force constants in kJ mol^-1 nm^-2",
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
        )
        print(render_json(checks) if args.json else render_text(checks))
        return 0 if all(check.ok for check in checks if check.name != "gmx") else 1
    if args.command == "build":
        allomorph_key = {
            "I-alpha": "i_alpha",
            "I-beta": "i_beta",
            "II": "ii",
            "III-I": "iii_i",
        }[args.allomorph]
        definition = load_allomorph(allomorph_key)
        if args.shape != "fibril" and any(
            value is not None for value in (args.padding_x, args.padding_y, args.padding_z)
        ):
            build_parser().error(
                "--padding-x, --padding-y, and --padding-z are valid only with --shape fibril"
            )
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
                    args.restraint_force,
                )
            ):
                build_parser().error("chain restraints require --force-field glycam or charmm36")
            args.output.parent.mkdir(parents=True, exist_ok=True)
            written_path = write_coordinates(crystal, args.output)
            if written_path != args.output:
                print(
                    f"Requested structure has {crystal.atom_count} atoms; "
                    f"wrote mmCIF instead of PDB: {written_path}"
                )
    return 0


def _write_charmm36_build(
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
    gro_path, itp_path, topology_path, force_field_path = write_charmm36_system(
        mapped,
        topology,
        force_field_directory,
        args.output,
        restrained_chain_count=len(restraints.chain_ids) if restraints else 0,
        restraint_force=restraints.force_constants if restraints else None,
    )
    report = {
        "cellulose_forge_version": __version__,
        "allomorph": crystal.allomorph,
        "shape": args.shape,
        "force_field": "CHARMM36",
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
        "outputs": {
            "coordinates": coordinate_path.name,
            "gromacs_coordinates": gro_path.name,
            "gromacs_molecule_topology": itp_path.name,
            "gromacs_topology": topology_path.name,
            "gromacs_force_field": force_field_path.name,
        },
        "geometry_annotations": list(crystal.annotations),
    }
    if restraints is not None:
        report["position_restraints"] = {
            "define": "POSRES_CELLULOSE",
            "selected_chain_ids": list(restraints.chain_ids),
            "force_constants_kj_mol_nm2": list(restraints.force_constants),
            "atom_selection": "all atoms in each selected chain",
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


def _write_glycam_build(
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
        restraints = _resolve_restraints(args, crystal, parser)
        chain_mapping = None
        if restraints is not None:
            mapped, chain_mapping = group_mapped_chains(mapped, restraints.chain_ids)
        output_directory.mkdir(parents=True, exist_ok=True)
        coordinate_path = write_coordinates(crystal, output_directory / "structure.pdb")
        gro_path, topology_path = write_glycam_gromacs_system(
            mapped,
            files.gromacs_topology,
            output_directory,
            restrained_chain_count=len(restraints.chain_ids) if restraints else 0,
            restraint_force=restraints.force_constants if restraints else None,
        )

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
        "outputs": {
            "coordinates": coordinate_path.name,
            "gromacs_coordinates": gro_path.name,
            "gromacs_topology": topology_path.name,
        },
        "geometry_annotations": list(crystal.annotations),
        "reference_generation": reference_metadata,
    }
    if restraints is not None:
        report["position_restraints"] = {
            "define": "POSRES_CELLULOSE",
            "selected_chain_ids": list(restraints.chain_ids),
            "force_constants_kj_mol_nm2": list(restraints.force_constants),
            "atom_selection": "all atoms in each selected chain",
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
        msg = "GLYCAM output requires equal, nonempty cellulose chain lengths"
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
    if args.fibril_model is None:
        parser.error("--fibril-model 36-chain is required for fibrils")
    if args.degree_polymerization is None:
        parser.error("--degree-polymerization N is required for fibrils")
    if args.chain_boundary is None:
        parser.error("--chain-boundary finite|periodic is required for fibrils")
    if args.chain_boundary == "periodic" and args.padding_z not in {None, 0}:
        parser.error("--padding-z must be zero with --chain-boundary periodic")
    return build_36_chain_fibril(
        definition,
        degree_polymerization=args.degree_polymerization,
        vacuum=args.vacuum,
        padding_x=args.padding_x,
        padding_y=args.padding_y,
        padding_z=args.padding_z,
        chain_boundary=args.chain_boundary,
    )


def _build_slab(
    args: argparse.Namespace,
    definition: AllomorphDefinition,
    parser: argparse.ArgumentParser,
) -> BulkCrystal:
    """Resolve mutually exclusive slab sizing modes and build the slab."""
    if args.cells is not None:
        parser.error("--cells is not valid with --shape slab")
    if args.face is None:
        parser.error("--face H K L is required for slabs")
    if args.chain_boundary is None:
        parser.error("--chain-boundary finite|periodic is required for slabs")
    count_values = (args.layers, args.chains_per_layer, args.degree_polymerization)
    size_values = (args.slab_width, args.slab_thickness, args.slab_length)
    has_count_value = any(value is not None for value in count_values)
    has_size_value = any(value is not None for value in size_values)
    if has_count_value and has_size_value:
        parser.error("exact-count and physical-size slab options cannot be combined")
    if has_count_value:
        if not all(value is not None for value in count_values):
            parser.error(
                "--layers N, --chains-per-layer N, and --degree-polymerization N "
                "must be supplied together"
            )
        assert args.layers is not None
        assert args.chains_per_layer is not None
        assert args.degree_polymerization is not None
        return build_chain_parallel_slab(
            definition,
            MillerIndex(*args.face),
            layers=args.layers,
            chains_per_layer=args.chains_per_layer,
            degree_polymerization=args.degree_polymerization,
            vacuum=args.vacuum,
            chain_boundary=args.chain_boundary,
        )
    if has_size_value:
        if not all(value is not None for value in size_values):
            parser.error(
                "--slab-width A, --slab-thickness A, and --slab-length A must be supplied together"
            )
        assert args.slab_width is not None
        assert args.slab_thickness is not None
        assert args.slab_length is not None
        return build_slab_by_dimensions(
            definition,
            MillerIndex(*args.face),
            width=args.slab_width,
            thickness=args.slab_thickness,
            length=args.slab_length,
            vacuum=args.vacuum,
            chain_boundary=args.chain_boundary,
        )
    parser.error(
        "supply either exact counts (--layers, --chains-per-layer, "
        "--degree-polymerization) or physical dimensions (--slab-width, "
        "--slab-thickness, --slab-length)"
    )


if __name__ == "__main__":
    raise SystemExit(main())
