# Experimental OPLS-CM5 backend

This backend is under development on the `experimental/opls-cm5` branch. It is
not part of the validated CelluloseForge 1.0 feature set or the current
SoftwareX claims.

The first supported construction is a pristine cellulose I-beta slab with
finite DP20 chains. It consumes the pristine `CNC.itp` and matching `solute.gro`
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
writing GROMACS files. Unsupported degrees of polymerization, periodic chains,
other allomorphs, and malformed/non-neutral references fail explicitly.

The proof-of-concept 2-layer by 4-chain slab passes `gmx grompp` with GROMACS
2024.2. A 100-step steepest-descent check completed without non-finite forces,
ending at 17,484.379 kJ/mol with a maximum force of 2,324.137 kJ mol-1 nm-1.
This is an engineering smoke test, not yet a scientific validation result.
