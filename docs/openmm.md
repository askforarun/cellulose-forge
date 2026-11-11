# OpenMM output

CelluloseForge writes complete OpenMM systems for GLYCAM06j-1 and CHARMM36
builds, with finite or covalently periodic chains:

```bash
cellulose-forge build \
  --allomorph I-beta --shape slab --face 1 1 0 \
  --layers 3 --chains-per-layer 4 --degree-polymerization 8 \
  --chain-boundary periodic --vacuum 30 \
  --force-field glycam --md-engine openmm \
  --output ibeta_110_openmm
```

Use `--force-field charmm36` for CHARMM36. Repeat `--md-engine` to write GROMACS
files alongside the OpenMM system, for example
`--md-engine gromacs --md-engine openmm`.

## Generated files

- `system.xml`: the complete serialized OpenMM `System`, including every
  force-field term, X--H bond constraints, any position restraints, and the
  periodic settings required by covalently periodic chains;
- `state.xml`: starting positions and periodic box vectors;
- `topology.cif`: PDBx/mmCIF atom, residue, and chain identities of the
  force-field-mapped system, including terminal atoms added for finite chains;
- `run_openmm.py`: a standalone minimization and dynamics runner;
- `build_report.json`: OpenMM version, per-force initial energies, and checks.

The runner needs only OpenMM. It does not read GROMACS or Amber files and does
not require GROMACS to be installed.

## How the system is constructed

CelluloseForge first writes the GROMACS-format topology it already validates
for that build, then reads it with OpenMM's `GromacsTopFile` reader. When only
`--md-engine openmm` is requested, these GROMACS files are temporary and are not
kept. GROMACS itself is never run.

For covalently periodic chains, the O4--C1 bond joining the last residue to the
first crosses the chain-direction box boundary, as do the angles, torsions, 1--4
pairs, and exclusions around it. CelluloseForge therefore enables periodic
treatment for every bonded force and for `NonbondedForce` exceptions:

- `HarmonicBondForce`, `HarmonicAngleForce`, and `PeriodicTorsionForce` use
  `setUsesPeriodicBoundaryConditions(True)`;
- `NonbondedForce` uses `setExceptionsUsePeriodicBoundaryConditions(True)`.

The exception setting matters even though it causes no obvious failure when
omitted: without it, boundary-crossing 1--4 pairs are evaluated across the whole
box, their energy silently approaches zero, and the PME exclusion corrections
are wrong. Because the settings are stored in `system.xml`, users cannot forget
them.

## Build-time checks

Every build evaluates the unconstrained system on OpenMM's Reference platform
and records the energy of each force, plus the 1--4 exception energy alone. The
build fails if

- any energy is non-finite;
- any bonded distance, evaluated with the minimum image, deviates by more than
  0.05 nm from its equilibrium length;
- the initial bonded energy exceeds 50 kJ mol^-1 per atom, the same criterion
  as the GROMACS integration tests;
- periodic exceptions were not enabled for a covalently periodic build; or
- the box is narrower than twice the 1 nm nonbonded cutoff.

A ten-iteration minimization of the constrained system is also recorded.

## Nonbonded settings

The serialized system uses PME with a 1 nm cutoff and OpenMM's default PME
accuracy (`ewaldErrorTolerance` 5 × 10^-4). The isotropic long-range dispersion
correction is enabled for bulk crystals and unit cells, and disabled for slabs
and fibrils, whose vacuum regions violate its homogeneous-system assumption.
The build report records the cutoff, the dispersion setting, and both PME
tolerances. Build-time energies are evaluated with a tighter tolerance
(1 × 10^-5, GROMACS's default `ewald-rtol`) so that they can be compared with
other engines.

## Running the system

```bash
python run_openmm.py --platform CUDA --minimize-iterations 1000 --steps 50000
```

The runner uses a Langevin-middle integrator (default 300 K, 1 ps^-1 friction,
2 fs time step) and writes `trajectory.dcd`, `state_data.csv`,
`final_state.xml`, and `final.cif`. Position restraints selected with
`--restrain-chains` or `--restrain-layers` are active by default; switch them
off with `--no-restraints`. These defaults are a starting point rather than a
production protocol.

### Pressure coupling

`--pressure-bar P` adds a `MonteCarloAnisotropicBarostat`. By default it scales
every Cartesian box axis except the chain axis recorded in the build report;
`--pressure-axes` selects axes explicitly, for example `x` for a slab whose
vacuum lies along `z`.

OpenMM's Monte Carlo barostats move molecules as rigid units. A covalently
periodic chain spans the box along the chain axis, so scaling that axis would
stretch or compress the boundary-crossing O4--C1 bond. The runner therefore
refuses to scale the chain axis for periodic chains, and it requires a
rectangular box.

Position restraints hold atoms to fixed Cartesian reference coordinates, which
a barostat does not rescale. The runner therefore refuses `--pressure-bar` while
restraints are active; add `--no-restraints` to run with pressure coupling.

| Ensemble | Status |
| --- | --- |
| NVT, with or without restraints | Supported |
| NPT with the chain axis fixed (rectangular box, no restraints) | Supported; short-run validated |
| Fully anisotropic NPT of periodic crystals | Use GROMACS |

## Validation

`validation/validate_openmm_engine.py` builds eleven systems: finite and
periodic Iβ (110) slabs with GLYCAM06j-1 and CHARMM36, a periodic triclinic
GLYCAM Iβ bulk crystal, and periodic Iα (010), II (110), and III_I (110) slabs
with both force fields. For each it compares OpenMM and GROMACS zero-step
energies and runs 1 ps of NVT; the periodic Iβ slabs also run 1 ps of NPT with
the chain axis fixed.

With OpenMM 8.5 and GROMACS 2024.2:

- bonded terms (bonds, angles, Urey-Bradley, and torsions) and the 1--4 terms
  agreed within 8 × 10^-6 relative. These are force-field definitions, so the
  agreement is limited only by GROMACS output precision;
- the total nonbonded energy (short-range Lennard-Jones and Coulomb, PME
  reciprocal, 1--4, and, for the bulk crystal, dispersion correction) agreed
  within 1.3 × 10^-4 relative when GROMACS used matching settings: plain
  cutoffs without potential shifts, `ewald-rtol` 1 × 10^-5, and the same
  dispersion correction. The remaining difference comes from the engines'
  different PME grids;
- after dynamics, every bonded distance stayed within 0.026 nm of its
  equilibrium length, including the boundary-crossing glycosidic bonds.

The results file records the CelluloseForge commit, OpenMM and GROMACS
versions, and every per-category value. These are engine-equivalence and
short-run integration checks, not material-property validation.

The OpenMM integration tests run only when OpenMM is importable by the test
interpreter or by the interpreter named in `CELLULOSE_FORGE_OPENMM_PYTHON`.

## Installation

```bash
conda env create -f environment.yml
# or
pip install 'cellulose-forge[openmm]'
```

Pass `--openmm-python /path/to/python` if OpenMM is installed under a different
interpreter.
