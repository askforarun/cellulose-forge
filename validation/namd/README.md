# NAMD validation bundle

This folder keeps the reusable scripts for the experimental NAMD backend in
one place. Generated simulation folders remain self-contained and do not
depend on Slurm or another scheduler.

- `compare_bundle.py` independently parses `system.psf`, expands the generated
  GROMACS topology with ParmEd, and compares atom properties and bonded term
  sets. It also records the parameter checksum and periodic cell widths.
- Every generated NAMD/GROMACS build contains the sole maintained runner,
  `run_single_point_comparison.sh`, its two engine inputs, and
  `single_point_energy.py`. There is deliberately no second copy here that can
  drift from the generated script.

Compare one build locally:

```bash
python validation/namd/compare_bundle.py BUILD_DIRECTORY \
  --output comparison.json
```

Run a generated folder on any workstation or cluster login/compute node where
the two executables are available:

```bash
cd BUILD_DIRECTORY
NAMD_EXE=namd3 GMX_EXE=gmx ./run_single_point_comparison.sh
```

Generate comparison-capable folders with both `--md-engine namd` and
`--md-engine gromacs`. The recommended starting bulk cells are I-beta
cellulose `4 4 3`, triclinic I-alpha cellulose `4 5 6`, and alpha-chitin
`7 2 3`; each clears the strict 28 Angstrom NAMD pair-list gate.

Run all three cases with real executables using:

```bash
python validation/namd/run_acceptance_matrix.py \
  --parameter-file /path/to/par_all36_carb.prm \
  --audit-python /path/to/ambertools/bin/python \
  --output validation/namd_acceptance_matrix
```

The launcher refuses missing executables or parameters and stops at the first
failed build, topology audit, or energy comparison. It does not create mock
NAMD output.

The stored acceptance matrix now records successful NAMD 3.0.3 single-point
comparisons against GROMACS for cellulose Ibeta, triclinic cellulose Ialpha,
and alpha-chitin. These comparisons use matched, unswitched 12 angstrom
Lennard-Jones settings to isolate topology and engine consistency; production
NAMD inputs retain the native CHARMM36 10--12 angstrom force switch. The
backend remains experimental until short production-input dynamics have also
been validated.
