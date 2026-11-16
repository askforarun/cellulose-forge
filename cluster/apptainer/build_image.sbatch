#!/bin/bash
#SBATCH --job-name=cf-image
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail
flight start
flight env activate gridware
module load apps/apptainer/1.2.5
cd "${SLURM_SUBMIT_DIR}"
apptainer build --fakeroot cellulose-forge-1.3.0.sif CelluloseForge.def
apptainer inspect cellulose-forge-1.3.0.sif

