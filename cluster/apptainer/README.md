# Private Apptainer validation bundle

This directory creates a self-contained CelluloseForge 1.3.0 validation image.
The prepared folder includes the exact Git commit and the user's local CHARMM36
distribution. It is private and must not be published or redistributed.

From the repository root on the Mac:

```bash
bash cluster/apptainer/prepare_bundle.sh
```

Copy the complete `cluster/apptainer/` directory to the cluster. Then run:

```bash
cd apptainer
mkdir -p logs results
image_job=$(sbatch --parsable build_image.sbatch)
sbatch --dependency=afterok:${image_job} validation_array.sbatch
```

The image job downloads and installs dependencies once. The validation array
runs at most three jobs concurrently and writes JSON and JUnit files under
`results/task_*`. No host Conda environment or molecular-simulation package is
used at run time; only the cluster's Apptainer module is required.

Copy `results/` and the Slurm logs back to the Mac after completion. Do not
publish the image because it embeds the private CHARMM36 distribution.
