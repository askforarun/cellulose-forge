# NAMD validation bundle

This folder keeps the reusable scripts for the experimental NAMD backend in
one place. Generated simulation folders remain self-contained; copy a build
folder plus `run_namd.sbatch` to a NAMD installation or Slurm cluster.

- `compare_bundle.py` independently parses `system.psf`, expands the generated
  GROMACS topology with ParmEd, and compares atom properties and bonded term
  sets. It also records the parameter checksum and periodic cell widths.
- `run_cross_engine.sbatch` runs NAMD, then performs a matching GROMACS t=0
  calculation when the build contains both engines. The generated build folder
  also contains this script, its matching MDP and `compare_t0.py`.

Compare one build locally:

```bash
python validation/namd/compare_bundle.py BUILD_DIRECTORY \
  --output comparison.json
```

Run a generated folder on Slurm:

```bash
cp validation/namd/run_cross_engine.sbatch BUILD_DIRECTORY/
cd BUILD_DIRECTORY
sbatch run_cross_engine.sbatch
```

The stored release validator will eventually build the full matrix, run NAMD,
and compare term-resolved energies against double-precision GROMACS and
OpenMM. Until NAMD runtime validation is complete, the CLI and documentation
label this backend experimental.
