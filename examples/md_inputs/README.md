# Simulation-ready example inputs

Small systems for trying each MD engine. Every folder holds the files that
`cellulose-forge build` wrote, including its starting run inputs:
`minim.mdp` and `nvt.mdp` for GROMACS, the configurable `run_openmm.py` for
OpenMM, and `min.in` and `nvt.in` for Amber. Every force field follows the
same rule for long-range dispersion correction: GLYCAM06 and OpenFF Sage use it
only for bulk crystals and unit cells; CHARMM36 never uses it; slabs and
fibrils never use it. `build_report.json` in each folder records the exact
build and its validation; machine-local paths are redacted.

| Folder | System | Atoms | Chains x DP | GROMACS | OpenMM | LAMMPS | Amber |
|---|---|---|---|---|---|---|---|
| `cellulose_ibeta_110_charmm36_periodic` | I-beta (1 1 0) slab, CHARMM36, periodic chains | 1008 | 8 x 6 | yes | yes | yes | - |
| `chitin_beta_100_charmm36_periodic` | beta-chitin (1 0 0) slab, CHARMM36, periodic chains | 1296 | 8 x 6 | yes | yes | yes | - |
| `chitin_alpha_100_charmm36_finite` | alpha-chitin (1 0 0) slab, CHARMM36, finite chains | 1320 | 8 x 6 | yes | yes | yes | - |
| `cellulose_ibeta_110_glycam_finite` | I-beta (1 1 0) slab, GLYCAM06j-1, finite chains | 1548 | 12 x 6 | yes | yes | yes | yes |
| `cellulose_ibeta_openff_sage_periodic_bulk` | I-beta bulk crystal, OpenFF Sage 2.3, periodic chains | 1512 | 18 x 4 | yes | yes | - | - |

Amber output is written only for finite chains, because Amber does not apply
minimum images to bonded terms across a periodic chain boundary.

These are short demonstration systems, not production models: the slabs are
two layers thick and DP 6, and the run inputs are smoke-test protocols. Choose
system size, equilibration, and production settings for your own study.
The OpenFF bulk folder additionally contains a 1 ns gently coupled
Parrinello--Rahman anisotropic NPT
equilibration template and a 100 ns Parrinello--Rahman anisotropic production
template. These are starting protocols for the planned stability assessment;
their presence is not a stability result. The isotropic compressibility of
water (4.5e-5 bar^-1) is used only to set the barostat response, as is
conventional; it is not a property of crystalline cellulose.

## Rebuilding

Each folder can be regenerated with `cellulose-forge build`; its
`build_report.json` records the allomorph, face, layers, chains, DP, chain
boundary, force field, and engines. For example,
`chitin_beta_100_charmm36_periodic`:

```bash
cellulose-forge build --polymer chitin --allomorph beta --shape slab \
  --face 1 0 0 --layers 2 --chains-per-layer 4 --degree-polymerization 6 \
  --chain-boundary periodic --vacuum 30 --force-field charmm36 \
  --md-engine gromacs --md-engine openmm --md-engine lammps \
  --output chitin_beta_100_charmm36_periodic
```

## CHARMM36 folders: install the force field first

`system.top` includes `charmm36.ff/forcefield.itp`. The CHARMM36 GROMACS port
is third-party and is not redistributed here. Install the pinned release once
and link it into the example folder:

```bash
cellulose-forge forcefield install charmm36
cd chitin_beta_100_charmm36_periodic
ln -s "$(cellulose-forge forcefield path charmm36)" charmm36.ff
```

(`cellulose-forge build` copies `charmm36.ff` into new outputs automatically.)
The OpenMM `system.xml` and LAMMPS `system.data` already contain every
parameter and do not need this step.

## Running

GROMACS (minimization, then 1 ps NVT):

```bash
gmx grompp -f minim.mdp -c structure.gro -p system.top -o minim.tpr
gmx mdrun -deffnm minim
gmx grompp -f nvt.mdp -c minim.gro -p system.top -o nvt.tpr
gmx mdrun -deffnm nvt
```

For the OpenFF periodic-bulk NPT example, continue from the generated NVT
checkpoint:

```bash
gmx grompp -f npt_equil.mdp -c nvt.gro -t nvt.cpt -p system.top -o npt_equil.tpr
gmx mdrun -deffnm npt_equil
gmx grompp -f npt_production.mdp -c npt_equil.gro -t npt_equil.cpt \
  -p system.top -o npt_production.tpr
gmx mdrun -deffnm npt_production
```

OpenMM (minimization, then Langevin dynamics):

```bash
python run_openmm.py --minimize-iterations 1000 --steps 5000
python run_openmm.py --help   # platform, temperature, timestep, reporters
```

LAMMPS (`in.cellulose` or `in.chitin` evaluates the initial energy with
`run 0`; add your own minimization or dynamics controls):

```bash
lmp -in in.chitin      # chitin examples
lmp -in in.cellulose   # cellulose examples
```

Amber (GLYCAM folder only; minimization, then 1 ps NVT):

```bash
sander -O -i min.in -p structure.prmtop -c structure.inpcrd -o min.out -r min.rst7
sander -O -i nvt.in -p structure.prmtop -c min.rst7 -o nvt.out -r nvt.rst7 -x nvt.nc
```

## Checked

Each folder was built with CelluloseForge on this branch, and every engine
input was run successfully: GROMACS preprocessing and minimization, OpenMM
minimization and 500 dynamics steps, LAMMPS `run 0`, and, for the GLYCAM
folder, Amber `sander` minimization. Term-by-term cross-engine agreement for
the CHARMM36 chitin systems is documented in
`validation/chitin_charmm36_engines.json`.
