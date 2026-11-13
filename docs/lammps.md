# Experimental LAMMPS output

The experimental implementation on `main` can write finite or covalently periodic
GLYCAM06j-1 and CHARMM36 cellulose systems for LAMMPS:

```bash
cellulose-forge build \
  --allomorph I-beta --shape slab --face 1 1 0 \
  --layers 2 --chains-per-layer 4 --degree-polymerization 8 \
  --chain-boundary periodic --vacuum 30 \
  --force-field glycam --md-engine lammps \
  --output ibeta_lammps
```

Select `--force-field charmm36` after installing the pinned CHARMM36 files to
write the corresponding CHARMM system:

Experimental OPLS-AA and OPLS-AA-SEI builds (`--force-field opls-aa` or
`opls-aa-sei`) write their LAMMPS files directly, without OpenMM; see
[`opls.md`](opls.md) for their styles and the tabulated SEI pair terms.

```bash
cellulose-forge forcefield install charmm36
cellulose-forge build \
  --allomorph I-beta --shape slab --face 1 1 0 \
  --layers 2 --chains-per-layer 6 --degree-polymerization 8 \
  --chain-boundary periodic --vacuum 30 \
  --force-field charmm36 --md-engine lammps \
  --output ibeta_charmm36_lammps
```

The output contains `system.data`, `in.cellulose` (`in.chitin` for
`--polymer chitin`), crystallographic
coordinates, and `build_report.json`. The generated LAMMPS model uses `real`
units, `atom_style full`, harmonic bonds and angles, CHARMM-style periodic
torsions with zero dihedral 1--4 weights,
`lj/cut/coul/long`, PPPM electrostatics, and full GLYCAM 1--4 interactions.
CHARMM36 output instead uses the CHARMM-specific LAMMPS styles
`angle_style charmm`, `dihedral_style charmmfsw`, and
`pair_style lj/charmmfsw/coul/long`. It preserves Urey--Bradley terms,
pair-specific 1--4 Lennard--Jones parameters, full 1--4 electrostatics, and
the 10--12 angstrom CHARMM force switch.

OpenMM is used as the topology parser, so install the `openmm` extra or pass
`--openmm-python`. A LAMMPS executable is required only to run the generated
input; GROMACS is not required for file generation. `cellulose-forge doctor`
reports OpenMM-based file-generation readiness separately from LAMMPS runtime
readiness.

## Current boundary

- GLYCAM and CHARMM36 finite and periodic chains are implemented.
- Position restraints are not yet implemented.
- `in.cellulose` (or `in.chitin`) performs `run 0` for initial-energy validation. It is not a
  minimization or molecular-dynamics protocol; users must supply simulation
  controls appropriate to their scientific application.
- The backend rejects unsupported OpenMM force classes and non-equivalent 1--4
  exceptions rather than approximating them.

The periodic I-beta GLYCAM engineering case passes `run 0` with LAMMPS
22 Jul 2025. Its initial potential energy agrees with the equivalent OpenMM
system within approximately 4.5 × 10^-5 relative. Finite and periodic I-beta
CHARMM36 cases also pass `run 0`; their bond, angle/Urey--Bradley, and torsion
energies agree with GROMACS within 2 × 10^-5 relative and their total
potential energies within 5 × 10^-4 when both engines use the official
CHARMM36 10--12 angstrom force-switch settings.

That raw 1.8--2.5 × 10^-4 total-potential gap is almost entirely three LAMMPS
energy conventions, audited term by term in
`validation/validate_charmm36_engines.py --polymer cellulose`
(results in `validation/cellulose_charmm36_engines.json`) for I-beta (1 1 0)
finite and periodic, I-alpha (0 1 0) periodic, II (1 1 0) periodic, and
III-I (1 1 0) periodic slabs:

- LAMMPS uses the CHARMM Coulomb constant (332.0716 vs GROMACS 332.0637,
  a 2.4 × 10^-5 proportional difference in electrostatic energies and forces);
- `dihedral_style charmmfsw` adds a constant energy offset to each 1--4 pair;
- CHARMM vfswitch (`lj/charmmfsw`) and GROMACS force-switch differ by a
  constant potential shift.

The two Lennard--Jones differences are constant offsets and do not affect
forces. Together the three account for 85--95% of the raw gap; after removing
them, total potentials agree within 2.7 × 10^-5 and nonbonded totals within
1.3 kJ/mol, and bonds, angles, Urey--Bradley, proper dihedrals, 1--4
Lennard--Jones, and 1--4 electrostatics each agree within 0.2 kJ/mol. The
same audit for chitin is in `validation/chitin_charmm36_engines.json`. Each
results file records the CelluloseForge commit, whether the tracked tree was
clean, and the GROMACS, OpenMM, and LAMMPS versions used.

Reproduce the cross-engine checks with:

```bash
uv run python validation/validate_lammps_engine.py \
  --openmm-python /path/to/python-with-openmm \
  --lammps-executable /path/to/lmp
```

The retained engine-specific inputs, logs, energy files, and comparison summary
are under `validation/results/lammps_cross_engine_files/`.

LAMMPS may report `Inconsistent image flags` for a covalently periodic chain:
the molecular graph winds around the periodic box and therefore cannot be
assigned one globally consistent set of unwrapped image flags. Bonded forces
still use the correct minimum image; the cross-engine energy check includes
these boundary-crossing interactions. However, analyses based on unwrapped
coordinates, including `dump xu yu zu` and a direct `compute msd`, are not
reliable for these chain-spanning molecules. Reconstruct or unwrap the
trajectory with a polymer-aware procedure before using those analyses.

The generated input ends with `quit 0` to terminate explicitly after the
single-point check. This avoids an exit-time fault observed in some packaged
LAMMPS builds and does not alter the computed energy.
