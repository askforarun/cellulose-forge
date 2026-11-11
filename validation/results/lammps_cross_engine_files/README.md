# LAMMPS cross-engine validation files

These folders retain the inputs and outputs used to compare the experimental
CelluloseForge LAMMPS backend with independent OpenMM or GROMACS evaluations.
All engines evaluated the same coordinates and force-field topology within a
case.

## Cases

- `glycam_periodic_ibeta`: periodic I-beta slab evaluated by GROMACS, OpenMM,
  and LAMMPS.
- `charmm36_finite_ibeta`: finite, longitudinally padded I-beta slab evaluated
  by GROMACS and LAMMPS.
- `charmm36_periodic_ibeta`: covalently periodic I-beta slab evaluated by
  GROMACS and LAMMPS.
- `glycam_periodic_ii`: periodic cellulose II slab evaluated by GROMACS,
  OpenMM, and LAMMPS.
- `charmm36_periodic_ii`: covalently periodic cellulose II slab evaluated by
  GROMACS and LAMMPS.

Each case separates common build files from engine-specific files. The
CHARMM36 GROMACS directories contain the generated cellulose ITP, topology,
run parameters, log, and extracted energy table. The externally distributed `charmm36.ff` directory is not copied
into this bundle; install the checksum-pinned release with
`cellulose-forge forcefield install charmm36` before reproducing those runs.

The LAMMPS directories contain the generated `run 0` input and the LAMMPS log;
the OpenMM directories contain the standalone runner. Binary and regenerable
outputs (coordinates, `system.data`, compiled `.tpr` and `.edr` files, and
serialized OpenMM XML) are not retained; regenerate them with
`validation/validate_lammps_engine.py`. Local file-system paths recorded by
the engines are redacted as `<local-path>`. `comparison.json` records the extracted
potential energies and relative differences.

For covalently periodic chains, the LAMMPS log can report inconsistent image
flags because the polymer graph winds through the box. The energy comparison
includes the boundary-crossing interactions, but unwrapped-coordinate and
direct MSD analyses should not use these image flags.
