# Validation

This directory will contain reproducible comparisons with the legacy builders and GROMACS preprocessing, minimization, and simulation checks. Generated binary or large reference artifacts should be archived externally and identified here by version and checksum.

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
