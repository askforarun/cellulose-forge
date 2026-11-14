# Validation

This directory will contain reproducible comparisons with the legacy builders and GROMACS preprocessing, minimization, and simulation checks. Generated binary or large reference artifacts should be archived externally and identified here by version and checksum.

## High-throughput coordinate benchmark

`benchmark_coordinate_matrix.py` exercises the public CLI over a matrix of all
four supported allomorphs and three bulk sizes. Each of the 12 systems is built
twice as mmCIF; SHA-256 checksums test byte-for-byte determinism while elapsed
times and atom counts are recorded in
`results/high_throughput_coordinate_matrix.json`. Run it from the repository
root with:

```bash
uv run python validation/benchmark_coordinate_matrix.py
```

This is a coordinate-construction benchmark, not a force-field or molecular-
dynamics performance benchmark. Its purpose is to demonstrate reproducible
batch generation without requiring AmberTools or GROMACS.

## Chitin GLYCAM06j-1 validation

`validate_chitin_glycam.py` builds alpha- and beta-chitin periodic bulk
crystals and finite/periodic (100) slabs. All six cases run through GROMACS
preprocessing with zero warnings, minimization, and a one-ps NVT smoke test.
The validator also minimizes and runs a finite beta-chitin slab for one ps
with AmberTools `sander`. It checks the written periodic O4--C1 topology and
records the explicit amide 1--4 scaling pairs. On identical coordinates, it
also compares bonds, angles, proper plus improper dihedrals, 1--4 Lennard-Jones,
and 1--4 electrostatic energies between double-precision GROMACS and `sander`.
The largest difference after correcting for their documented Coulomb constants
is 0.000124 kcal/mol. Compact evidence is stored in
`results/chitin_glycam_engines.json`. These are topology and runtime tests, not an
assessment of GLYCAM force-field accuracy for crystalline chitin.

`validate_chitin_glycam_cross_engines.py` independently compares alpha- and
beta-chitin periodic bulk crystals and `(100)` slabs with finite and covalently
periodic chains in double-precision GROMACS, OpenMM, and LAMMPS. Bonds, angles,
torsions, ordinary
1--4 terms, and the amide-specific scaled 1--4 terms are compared separately.
The largest bonded or explicit-pair relative difference is 8.4 × 10^-6; the
largest combined nonbonded difference is 1.6 × 10^-4, within the recorded
PME--PPPM tolerance. All six OpenMM systems also complete a two-ps NVT smoke
test. Results are in `results/chitin_glycam_cross_engines.json`.

## GROMACS protocols

- `periodic_pme_minim.mdp`: energy minimization for cellulose chains bonded
  through the periodic chain-direction boundary.
- `finite_pme_minim.mdp`: energy minimization for chemically capped finite
  chains. The simulation box must leave sufficient separation between chain
  ends when three-dimensional periodic boundaries are used.
- `periodic_pme_nvt.mdp`: 10 ps, 300 K NVT smoke test for periodic-chain
  systems. It uses a 1 fs timestep, PME electrostatics, hydrogen-bond LINCS
  constraints, and the velocity-rescale thermostat. This short run tests
  numerical and topology stability; it is not an equilibration or production
  protocol.
- `periodic_pme_npt.mdp`: 10 ps, 300 K and 1 bar bulk-crystal pressure-
  equilibration smoke test. It continues from NVT velocities and uses
  isotropic stochastic cell rescaling, preserving the experimental monoclinic
  box shape while allowing its volume to respond. Independent relaxation of
  lattice lengths and angles requires a longer, validated anisotropic
  production protocol.
- `position_restraint_minim.mdp`: minimization smoke test with selective chain
  restraints enabled through `POSRES_CELLULOSE`.
- `position_restraint_pme_nvt_1ps.mdp`: deterministic 1 ps, 300 K PME/NVT
  integration test for selectively restrained cellulose chains.

## GLYCAM06j-1 and CHARMM36 force-field matrix

`validate_force_field_matrix.py` provides one reproducible validation matrix
for the two production force fields. It covers periodic bulk crystals of all
four cellulose allomorphs and the documented I-alpha (010), I-beta (110), II
(110), and III-I (110) slabs with finite and covalently periodic chains. The
matrix therefore contains 8 bulk systems and 16 slab systems. Every case must
pass `gmx grompp -maxwarn 0`, minimization, and a 1 ps PME/NVT smoke test;
periodic builds must also pass the mapped-geometry and independently parsed
written-topology O4--C1 bond checks.

```bash
uv run python validation/validate_force_field_matrix.py \
  --tleap-executable /path/to/ambertools/bin/tleap \
  --parmed-python /path/to/ambertools/bin/python \
  --charmm-force-field-directory /path/to/charmm36.ff
```

The compact, provenance-bearing result is
`results/force_field_matrix.json`. Cellulose II records its construction of
the missing hydroxyl hydrogens. `validate_cellulose_ii_hydroxyls.py` separately
requires the exact major donor/acceptor assignments reported in Table 3 of
Langan, Nishiyama and Chanzy (1999) and records periodic contact distances and
angles.

## CHARMM36 validation

The CHARMM36 adapter reads `BGLC` from the checksum-pinned official GROMACS
force-field distribution and generates finite or arbitrary even-DP cyclic
topologies directly. A 4x4x6 I-beta DP12 build (32 chains,
8,064 atoms) passed `gmx grompp`, converged below the configured minimization
force threshold, and completed 1 ps of NVT integration. A generated DP8 system
also passed preprocessing and minimization.

The reference chemistry can be audited without copying external parameters:

```bash
uv run python validation/validate_charmm36_reference.py \
  --rtf /path/to/top_all36_carb.rtf \
  --itp /path/to/CARB_12.itp \
  --output validation/results/charmm36_reference.json
```

This validates the periodic repeat against CHARMM36 `BGLC` plus cyclic `14bb`,
checks neutrality and atom modifications, and compares the complete regenerated
DP12 bonded graph with the source ITP. See `docs/charmm36.md` for
the exact build CLI, engineering results, and current scientific limitations.
The companion `audit_charmm36_coefficients.py` exhaustively checks the
nonbonded, bond, angle/Urey--Bradley, and proper-dihedral parameter types used
by that topology against the upstream carbohydrate PRM with explicit unit and
functional-form conversions. Its compact result is stored in
`results/charmm36_coefficients.json`.

The periodic I-beta shape matrix additionally contains a (110) slab and the
literature-defined 36-chain fibril at DP8. Both passed PME preprocessing and
minimized below the configured force threshold without `-maxwarn`. Compact
results are stored in `results/charmm36_ibeta_shapes.json`; generated GROMACS
run files remain under the ignored `validation/work/` directory.

The current periodic CHARMM36 bulk and slab validation is part of the unified
24-system matrix in `results/force_field_matrix.json`.

The periodic 36-chain CHARMM36 fibril matrix covers the I-alpha, I-beta, and II
presets at DP8 (6,048 atoms each). Every fibril passes preprocessing, PME
minimization, and 1-ps NVT without LINCS warnings or non-finite energies. The
cellulose II mapping follows each chain's measured periodic O4--C1 path; this
avoids assuming that fibril chain direction can be inferred from identifier
parity. Compact evidence is stored in
`results/charmm36_fibril_matrix.json`.

Finite DP6 I-alpha, I-beta, and II fibrils and an I-beta (110) slab pass tight
PME minimization followed by 1-ps NVT integration without LINCS or non-finite
energy warnings. Compact results are stored in
`results/charmm36_finite_pme_nvt.json`. The independent DP15 CHARMM-GUI
comparison can be reproduced with:

```bash
uv run python validation/validate_finite_charmm36_reference.py \
  --force-field-directory /path/to/charmm36.ff \
  --itp /path/to/CHARMM-GUI/CARA.itp \
  --output validation/results/charmm36_finite_reference.json
```

The comparison matches all atom names, types, charges, and 2,883 bonded
interaction records exactly after normalization by residue and atom identity.

Nanopolysaccharide Builder was reviewed as an independent CHARMM36 workflow
reference. It exposes CHARMM36 primarily through PSF generation with `psfgen`;
CelluloseForge's experimental path instead targets native GROMACS GRO/ITP/TOP
output without making `psfgen` a runtime dependency.

## Slab allomorph matrix

The 16 finite/periodic slab cases are now generated by the unified validator
described above and stored with the bulk cases in
`results/force_field_matrix.json`. These are construction, topology, and
runtime smoke tests, not equilibrated material-property predictions.

The I-beta `(110)` lattice construction additionally has a site-assignment
regression requiring origin and center conformations to alternate within each
layer, no interchain contact below 1.7 angstrom, and rejection of incomplete
odd-chain transverse repeats. A four-layer, four-chain-per-layer periodic DP12
GLYCAM slab minimized in nine steps and completed a one-ps PME/NVT test without
LINCS warnings or non-finite energies. The two central layers retained a
minimum-image RMSD of 0.780 angstrom from the minimized structure. Results are
stored in `results/ibeta_110_lattice_validation.json`; this is a short
stability check rather than an equilibrium material-property prediction.

## I-beta fibril smoke test

The literature-defined DP12 36-chain I-beta fibril was mapped to the periodic
GLYCAM chain reference and tested with GROMACS 2024.2. `gmx grompp` completed
with no warnings using `periodic_pme_minim.mdp`; the resulting system contained
9,072 atoms and 432 residues. A version 1.1.0 rerun converged in nine steps to
a maximum force of 768.4 kJ mol-1 nm-1 (configured threshold: 1000).

Selective-restraint validation used the same DP12 fibril with chain IDs 6, 11,
16, 21, 26, and 31 assigned force constants of 1000 kJ mol-1 nm-2 along each
Cartesian direction. The generated topology contained 30 free and six
restrained molecules, covering 792 restrained heavy atoms. This workflow was
rerun with CelluloseForge 1.1.0 after the heavy-atom restraint and orientation
changes. GROMACS 2024.2 preprocessing succeeded without `-maxwarn`;
minimization reached the force threshold in 15 steps; and the 1 ps PME/NVT
test completed without LINCS warnings or non-finite energies. The trajectory
recorded a mean position-restraint energy of 269.591 kJ mol-1, confirming that
the generated interaction was active.
Compact results are retained in
`results/position_restraint_validation.json`; trajectory and run files remain
untracked because they are reproducible binary artifacts.

## Large I-beta slab preprocessing test

An I-beta `(110)` slab containing four layers, 100 chains per layer, and 12
anhydroglucose units per chain was built with finite chain ends and 30 Å of
vacuum. The crystallographic model contained 100,800 atoms. Addition of the
GLYCAM terminal atoms produced a 102,000-atom, 5,200-residue GROMACS system
comprising 400 cellulose chains.

CelluloseForge wrote an 8.8 MB `structure.cif` with unrestricted identifiers, a
4.4 MB `structure.gro`, `system.top`, and `build_report.json`. GROMACS 2024.2
processed the GRO/TOP pair using `finite_pme_minim.mdp`, PME electrostatics, and
no `-maxwarn` allowance. `gmx grompp` completed with zero warnings, generated
all 21 nonbonded and 21 one-four parameter combinations, and produced a valid
run-input file. This validates preprocessing and file interoperability; it is
not a claim that a production simulation of the large slab has been completed.

## OpenMM engine validation

`validate_openmm_engine.py` builds eleven systems with both GROMACS and OpenMM
output from the same mapped structure: finite and covalently periodic I-beta
`(110)` slabs (three layers of four DP8 chains) with GLYCAM06j-1 and CHARMM36, a
periodic triclinic GLYCAM I-beta bulk crystal (3 × 3 × 2 cells, DP4), and
periodic I-alpha `(010)`, II `(110)`, and III-I `(110)` slabs (three layers of
six DP8 chains) with both force fields. For each system it compares OpenMM and
GROMACS zero-step energies and runs the generated `run_openmm.py` for 1 ps of
NVT; the periodic I-beta slabs also run 1 ps of NPT with the barostat scaling
only the in-plane `x` axis.

```bash
uv run python validation/validate_openmm_engine.py \
  --openmm-python /path/to/python-with-openmm \
  --tleap-executable /path/to/tleap --parmed-python /path/to/python-with-parmed
```

Two kinds of energy agreement are recorded separately:

- Bonded and 1-4 categories (bonds, angles, torsions, and combined
  Lennard-Jones and Coulomb 1-4 terms; CHARMM36 bonds, angles, and Urey-Bradley
  terms are compared together) agreed within 8 × 10^-6 relative.
- The total nonbonded energy agreed within 1.3 × 10^-4 relative. For this
  comparison GROMACS used plain cutoffs (`vdw-modifier` and `coulomb-modifier`
  `None`), `ewald-rtol` 1 × 10^-5, and the same dispersion correction as the
  OpenMM system: on for the bulk crystal and off for slabs. The residual
  difference reflects the engines' different PME grids.

GROMACS used `periodic-molecules = yes` for the periodic systems; OpenMM used
periodic bonded forces and periodic nonbonded exceptions. After NVT and
fixed-chain-axis NPT, every bonded distance remained within 0.026 nm of its
equilibrium length, so the boundary-crossing glycosidic bonds stayed intact.
The results file, `results/openmm_engine_validation.json`, records the
CelluloseForge commit (`c0574c3`, clean working tree), OpenMM 8.5, and GROMACS
2024.2. These are engine-equivalence and short-run integration checks, not
material-property validation. GitHub Actions installs OpenMM and GROMACS in a
dedicated Linux job and runs the four finite/periodic, GLYCAM/CHARMM36
cross-engine integration tests for every push and pull request.

## Experimental OPLS-AA and OPLS-AA-SEI matrix

`validate_opls_cellulose.py` covers both variants for finite and periodic
slabs of I-alpha `(010)`, I-beta, II and III-I `(110)`, periodic bulk crystals
of all four allomorphs, and finite and periodic 36-chain fibrils of I-alpha,
I-beta and II (36 systems in total):

```bash
python validation/validate_opls_cellulose.py \
  --openmm-python /path/to/python-with-openmm \
  --gmx-double /path/to/gmx_d
```

It requires mixed-precision `gmx`, double-precision `gmx_d`, `lmp_serial` and
an OpenMM interpreter, and checks:

- transcription of the Kony et al. (2002) Table 3 OPLS-AA torsions and every
  explicit bond and angle parameter against the installed GROMACS `oplsaa.ff`;
- torsion, 1-4 and SEI-pair energies recomputed in NumPy against
  double-precision GROMACS, and the OPLS-AA-SEI minus OPLS-AA nonbonded energy
  against the analytic effect of the 1-5/1-6 scaling;
- bond, angle, torsion, 1-4/SEI-pair, Lennard-Jones and total energies of
  double-precision GROMACS, OpenMM and LAMMPS on identical, unrounded
  coordinates (GROMACS reads a full-precision `.g96`; a GRO file would place
  pairs exactly on the cutoff) with unshifted cutoffs, no dispersion
  correction and tight Ewald settings: terms within 10^-6 and totals within
  10^-5 relative;
- `gmx grompp -maxwarn 0`, minimization and 1 ps NVT of the shipped GROMACS
  files, minimization and 1 ps with the shipped `run_openmm.py`, and the
  shipped LAMMPS input.

I-alpha bulk boxes have a tilted third vector. For these boxes the OpenMM 8.6.1
Reference platform disagreed with an exact nearest-image Lennard-Jones sum by
0.14--0.40 % in local measurements, while the CPU platform agreed within
10^-6 and untilted controls agreed on Reference within 10^-13
(`measure_openmm_tilted_box.py`, `results/openmm_tilted_box.json`). Those
cases are evaluated on the CPU platform with a 10^-5 tolerance for the OpenMM
comparison.

LAMMPS is compared in every category it reports: harmonic bonds and the SEI
pseudo-bond tables (separated with `compute bond`), angles, torsions, pair-style
Lennard-Jones (ordinary pairs plus 0.5-scaled 1-4 pairs), pair-style plus
long-range Coulomb, and the total. Its references are the GROMACS terms minus
the SEI pair contributions, which are recomputed independently in NumPy and
checked against GROMACS LJ-14/Coulomb-14. Ewald-dependent categories (Coulomb
and total) are normalized by the magnitude of the total potential, because
Ewald errors are absolute and a net Coulomb energy can be a small difference
of large terms. Results: `results/opls_cellulose_matrix.json`, with provenance.
These checks validate topology construction and engine translation, not
force-field accuracy.
