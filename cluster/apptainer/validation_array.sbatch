#!/bin/bash
#SBATCH --job-name=cf-validate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-9%3
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -euo pipefail
flight start
flight env activate gridware
module load apps/apptainer/1.2.5
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs results
export APPTAINERENV_OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
apptainer exec --cleanenv --bind "${SLURM_SUBMIT_DIR}/results:/results" cellulose-forge-1.3.0.sif /opt/cellulose-forge/cluster/apptainer/run_validation_task.sh "${SLURM_ARRAY_TASK_ID}" /results

