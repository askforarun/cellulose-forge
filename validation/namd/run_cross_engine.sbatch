#!/bin/bash
#SBATCH --job-name=cf_namd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

# Set these only if modules are required by the target cluster.
[[ -z "${NAMD_MODULE:-}" ]] || module load "${NAMD_MODULE}"
[[ -z "${GROMACS_MODULE:-}" ]] || module load "${GROMACS_MODULE}"

NAMD_EXE="${NAMD_EXE:-namd3}"
GMX_EXE="${GMX_EXE:-gmx}"

"${NAMD_EXE}" +p"${SLURM_CPUS_PER_TASK}" run_namd.conf > namd.log

if [[ -f system.top && -f structure.gro ]]; then
  "${GMX_EXE}" grompp -f gromacs_namd_t0.mdp -c structure_namd.pdb \
    -p system.top -o gromacs_t0.tpr -po gromacs_t0_processed.mdp -maxwarn 0
  "${GMX_EXE}" mdrun -s gromacs_t0.tpr -deffnm gromacs_t0 -nsteps 0
  printf '%s\n' Bond U-B Angle Proper-Dih. Improper-Dih. LJ-14 Coulomb-14 \
    LJ-SR Coulomb-SR Coul.-recip. Potential 0 | \
    "${GMX_EXE}" energy -f gromacs_t0.edr -o gromacs_t0_energy.xvg
  python compare_t0.py namd.log gromacs_t0_energy.xvg > t0_energy_comparison.json
fi

