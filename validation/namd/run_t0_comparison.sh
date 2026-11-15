#!/bin/bash
set -euo pipefail

NAMD_EXE="${NAMD_EXE:-namd3}"
GMX_EXE="${GMX_EXE:-gmx}"
NAMD_THREADS="${NAMD_THREADS:-8}"

"${NAMD_EXE}" +p"${NAMD_THREADS}" run_namd.conf > namd.log
"${GMX_EXE}" grompp -f gromacs_namd_t0.mdp -c structure_namd.pdb \
  -p system.top -o gromacs_t0.tpr -maxwarn 0
"${GMX_EXE}" mdrun -s gromacs_t0.tpr -deffnm gromacs_t0 -nsteps 0
printf 'Bond\nAngle\nU-B\nProper Dih.\nImproper Dih.\nLJ-14\nCoulomb-14\nLJ (SR)\nCoulomb (SR)\nCoul. recip.\nPotential\n0\n' | \
  "${GMX_EXE}" energy -f gromacs_t0.edr -o gromacs_t0_energy.xvg
python compare_t0.py namd.log gromacs_t0_energy.xvg > t0_energy_comparison.json
