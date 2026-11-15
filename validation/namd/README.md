# NAMD validation bundle

This folder keeps the reusable scripts for the experimental NAMD backend in
one place. Generated simulation folders remain self-contained and do not
depend on Slurm or another scheduler.

- `compare_bundle.py` independently parses `system.psf`, expands the generated
  GROMACS topology with ParmEd, and compares atom properties and bonded term
  sets. It also records the parameter checksum and periodic cell widths.
- `run_t0_comparison.sh` runs NAMD, then performs a matching GROMACS t=0
  calculation when the build contains both engines. The generated build folder
  also contains this script, its matching MDP and `compare_t0.py`.

Compare one build locally:

```bash
python validation/namd/compare_bundle.py BUILD_DIRECTORY \
  --output comparison.json
```

Run a generated folder on any workstation or cluster login/compute node where
the two executables are available:

```bash
cp validation/namd/run_t0_comparison.sh BUILD_DIRECTORY/
cd BUILD_DIRECTORY
NAMD_EXE=namd3 GMX_EXE=gmx ./run_t0_comparison.sh
```

The stored release validator will eventually build the full matrix, run NAMD,
and compare term-resolved energies against double-precision GROMACS and
OpenMM. Until NAMD runtime validation is complete, the CLI and documentation
label this backend experimental.
