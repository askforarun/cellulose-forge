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

## Matched GLYCAM and experimental OPLS-CM5 comparison

`compare_force_fields.py` automates a matched engineering comparison using the
same periodic cellulose I-beta 4x4x4 bulk geometry (32 chains, DP8, 5,376
atoms) for GLYCAM and the experimental OPLS-CM5 backend. It builds both
systems, runs PME minimization, NVT, and isotropic NPT, checks the total charge
stored in the compiled TPR, records the force field's 1--4 scaling, and writes
complete- and late-window summaries for temperature, pressure, density,
volume, box lengths, and potential energy.

The default simulation stages are short 10 ps engineering checks. Use
`--npt-steps 100000` for the 100 ps comparison protocol. The OPLS inputs are
external templates and must be supplied explicitly; they are not vendored or
treated as a public validation standard.

```bash
uv run python validation/compare_force_fields.py \
  --opls-template /path/to/CNC.itp \
  --opls-reference-gro /path/to/solute.gro \
  --tleap-executable /path/to/ambertools/bin/tleap \
  --parmed-python /path/to/ambertools/bin/python \
  --npt-steps 100000
```

Large trajectories and binary GROMACS run files are written under
`validation/work/` and ignored by git. Compact evidence from the September 20,
2026 development run is retained in
`results/force_field_comparison/summary.json`. OPLS completed 100 ps. The
matched GLYCAM run was intentionally stopped at 38.48 ps after runtime
stability had been established, so its late window is 20--38.48 ps rather than
50--100 ps. These trajectories demonstrate topology and runtime integrity;
they do not establish force-field accuracy or converged material properties.

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

## Experimental periodic OPLS-CM5 allomorph matrix

`validate_opls_allomorphs.py` applies the externally supplied OPLS-CM5 repeat
template to periodic bulk I-alpha, I-beta, II, and III-I crystals, then runs
the same `gmx grompp`, PME minimization, and configurable NVT engineering
checks as the GLYCAM matrix. The external finite DP20 topology and its matching
coordinates are required explicitly:

```bash
uv run python validation/validate_opls_allomorphs.py \
  --opls-template /path/to/CNC.itp \
  --opls-reference-gro /path/to/solute.gro
```

Compact evidence is stored in `results/opls_allomorph_matrix/`. Successful
preprocessing, minimization, and 1-ps integration establish that the generated
topologies are internally executable. They do not independently validate the
transferability or physical accuracy of the I-beta-derived OPLS-CM5 parameters
for other cellulose allomorphs.

## Matched GLYCAM/OPLS-CM5 allomorph NPT comparison

`compare_allomorph_force_fields.py` builds every periodic bulk allomorph with
both force fields, minimizes it, and runs matched NVT and isotropic NPT stages.
The default duration is 10 ps for each ensemble:

```bash
uv run python validation/compare_allomorph_force_fields.py \
  --opls-template /path/to/CNC.itp \
  --opls-reference-gro /path/to/solute.gro \
  --tleap-executable /path/to/ambertools/bin/tleap \
  --parmed-python /path/to/ambertools/bin/python
```

The driver separately records Lennard--Jones and electrostatic 1--4 scaling.
The current generated systems use `fudgeLJ = fudgeQQ = 1.0` for GLYCAM and
`fudgeLJ = fudgeQQ = 0.5` for OPLS-AA. Compact results are stored in
`results/allomorph_force_field_comparison/`; trajectories and binary run files
remain under ignored `validation/work/` paths.

## I-beta fibril smoke test

The literature-defined DP12 36-chain I-beta fibril was mapped to the periodic
GLYCAM chain reference and tested with GROMACS 2024.2. `gmx grompp` completed
with no warnings using `periodic_pme_minim.mdp`; the resulting system contained
9,072 atoms and 432 residues. Steepest descent converged in 10 steps to a
maximum force of 949.3 kJ mol-1 nm-1 (configured threshold: 1000).

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
