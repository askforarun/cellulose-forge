# NAMD validation bundle

This folder keeps the reusable scripts for the NAMD backend in one place.
Generated simulation folders remain self-contained and do not depend on Slurm
or another scheduler.

- `compare_bundle.py` independently parses `system.psf`, expands the generated
  GROMACS topology with ParmEd, and compares atom properties and bonded term
  sets. It also records the parameter checksum and periodic cell widths.
- `analyze_dynamics.py` checks completion, finite energies, temperature,
  minimum-image bond lengths, and heavy-atom RMSD after the generated
  minimization and 1 ps NVT smoke test.
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
failed build, topology audit, energy comparison, or dynamics smoke test. It
does not create mock NAMD output. Its live single-point check now runs the
shared schema-1 GROMACS and NAMD energy drivers on a neutral coordinate
reference, requires complete canonical decompositions and independent
potential-sum checks, and then compares every canonical term. The generated
shell runner remains available as a portable compatibility entry point until
the consolidated downloadable runner replaces it.

The stored acceptance matrix now records successful NAMD 3.0.3 single-point
comparisons against GROMACS for cellulose Ibeta, triclinic cellulose Ialpha,
and alpha-chitin. These comparisons use matched, unswitched 12 angstrom
Lennard-Jones settings to isolate topology and engine consistency; production
NAMD inputs retain the native CHARMM36 10--12 angstrom force switch. The same
matrix runs each generated production input through minimization and 1 ps NVT
and requires finite energies, intact bonds, and a heavy-atom RMSD below 2
angstrom.
