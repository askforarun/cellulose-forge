# Experimental OPLS-CM5 backend

This backend is under development on the `experimental/opls-cm5` branch. It is
not part of the validated CelluloseForge 1.0 feature set or the current
SoftwareX claims.

The supported constructions are pristine cellulose I-beta slabs and fibrils
with finite DP20 chains, plus even-DP periodic unit cells, bulk crystals, slabs,
and fibrils. They consume the pristine `CNC.itp` and matching `solute.gro`
from the reference data associated with Elahi et al. (2025). CelluloseForge does
not vendor those files because their repository does not currently state a
software/data license.

The pristine topology is available in the history of the authors' CNCstruc
repository, immediately before commit `6f20944` removed it:

```bash
git clone https://github.com/arashelahi/CNCstruc.git
git -C CNCstruc checkout 6f20944^ -- \
  data/simulation_traj_topol/CNC.itp \
  data/simulation_traj_topol/solute.gro
```

An experimental slab can then be built with:

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape slab \
  --face 1 1 0 \
  --layers 2 \
  --chains-per-layer 4 \
  --degree-polymerization 20 \
  --chain-boundary finite \
  --vacuum 30 \
  --force-field opls-cm5 \
  --opls-template CNCstruc/data/simulation_traj_topol/CNC.itp \
  --opls-reference-gro CNCstruc/data/simulation_traj_topol/solute.gro \
  --output ibeta_110_opls_cm5
```

The adapter verifies the exact 423-atom DP20 layout and zero total charge before
writing GROMACS files. Other allomorphs and malformed/non-neutral references
fail explicitly.

For a periodic bulk crystal, CelluloseForge derives the 21-atom internal AGU
repeat and its complete local bond, angle, Ryckaert--Bellemans dihedral, and
explicit 1--4-pair pattern from central residues of the finite reference. It
then instantiates those interactions cyclically, including the DP-to-1 boundary:

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape bulk \
  --cells 4 4 4 \
  --chain-boundary periodic \
  --force-field opls-cm5 \
  --opls-template CNCstruc/data/simulation_traj_topol/CNC.itp \
  --opls-reference-gro CNCstruc/data/simulation_traj_topol/solute.gro \
  --output ibeta_bulk_opls_cm5
```

The published rounded internal-repeat charges sum to +0.00163 e. For cyclic
chains, the adapter distributes the residual uniformly over all 21 repeat atoms
instead of assigning it to a chemically arbitrary atom. The source residual,
per-atom correction, corrected charge, and generated interaction counts are
recorded in `build_report.json`. This policy preserves equivalent-atom symmetry
but remains an experimental derivation requiring independent scientific review.

The proof-of-concept 2-layer by 4-chain slab passes `gmx grompp` with GROMACS
2024.2. A 100-step steepest-descent check completed without non-finite forces,
ending at 17,484.379 kJ/mol with a maximum force of 2,324.137 kJ mol-1 nm-1.
This is an engineering smoke test, not yet a scientific validation result.

A periodic 4x4x4-cell I-beta bulk crystal (32 chains, DP8, 5376 atoms) also
passes `gmx grompp` with OPLS 1--4 scaling (`fudge = 0.5`). PME minimization
converged in eight steps to a maximum force of 702.985 kJ mol-1 nm-1. One-ps
NVT and one-ps isotropic NPT smoke tests both completed normally. These short
runs establish topology and runtime integrity only; they do not validate bulk
thermodynamics or reproduce the published force-field study.
