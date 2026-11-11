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

## Periodic GLYCAM allomorph matrix

`validate_glycam_allomorphs.py` builds periodic bulk crystals for I-alpha,
I-beta, II, and III-I, then runs `gmx grompp`, PME minimization, and a
configurable NVT stability test. The default is 1 ps. III-I uses 6x4x4 cells
because its short `a` lattice vector requires six repeats for compatibility
with the validation cutoff.

```bash
uv run python validation/validate_glycam_allomorphs.py \
  --tleap-executable /path/to/ambertools/bin/tleap \
  --parmed-python /path/to/ambertools/bin/python
```

Compact evidence is stored in `results/glycam_allomorph_matrix/`. Cellulose II
records its force-field construction of the missing hydroxyl hydrogens and
uses opposite topology-to-crystal residue mappings for its antiparallel origin
and center chains.

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

The periodic CHARMM36 bulk matrix covers I-alpha, I-beta, II, and III-I at
DP8. Every system passes preprocessing, PME minimization, and 1-ps NVT without
LINCS or non-finite-energy failures. Cellulose II records its constructed
hydroxyl hydrogens and antiparallel residue mapping. Compact evidence is stored
in `results/charmm36_allomorph_matrix.json`.

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

Chain-parallel slab construction is exercised for I-alpha `(010)`, I-beta
`(110)`, II `(110)`, and III-I `(110)`. For each allomorph, finite and periodic
DP6 systems containing three layers and six chains per layer were generated
with both GLYCAM06j-1 and CHARMM36. Finite chains used 12 angstrom of total
chain-end padding. All 16 systems passed PME preprocessing without `-maxwarn`.
Finite systems minimized below 100 kJ mol-1 nm-1 and periodic systems below
1000 kJ mol-1 nm-1; every system then completed 1 ps NVT without LINCS
warnings or non-finite energies. Periodic O4--C1 linkages were checked before
output and were shorter than 2 angstrom. Compact results are stored in
`results/slab_allomorph_matrix.json`. These are topology and runtime smoke
tests, not equilibrated material-property predictions. The original complete
matrix was initially run at commit `73f6543`; the four I-beta finite/periodic
GLYCAM and CHARMM36 entries were regenerated and rerun with the final 1.1.2
implementation committed as `63675a2`, after the alternating lattice-site
validation was added. The remaining allomorph results exercise code paths not
changed by that I-beta-specific correction.

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
