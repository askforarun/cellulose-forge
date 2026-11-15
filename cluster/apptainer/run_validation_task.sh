#!/bin/bash
set -euo pipefail
task_id="${1:?task id required}"
output_root="${2:?output directory required}"
scratch="${SLURM_TMPDIR:-/tmp}/cellulose-forge-${SLURM_JOB_ID:-local}-${task_id}"
repo="${scratch}/repo"
mkdir -p "${scratch}" "${output_root}/task_${task_id}"
cp -a /opt/cellulose-forge "${repo}"
cd "${repo}"
core=/opt/envs/core/bin/python
openff=/opt/envs/openff/bin/python
gmx=/opt/envs/openff/bin/gmx
gmx_d=/opt/gromacs-double/bin/gmx_d
lmp=/opt/envs/openff/bin/lmp
ff=/opt/cellulose-forge-data/forcefields/charmm36-feb2026_cgenff-5.0.ff
work="${scratch}/work"
mkdir -p "${work}"
collect_results() {
  cp validation/results/*.json "${output_root}/task_${task_id}/"
  git restore validation/results
}

case "${task_id}" in
  0) "${core}" -m pytest --junitxml="${output_root}/task_0/core.xml" ;;
  1) "${openff}" -m pytest --junitxml="${output_root}/task_1/openff.xml" ;;
  2) "${core}" validation/validate_force_field_matrix.py --gmx "${gmx}" --tleap-executable /opt/envs/core/bin/tleap --parmed-python "${core}" --charmm-force-field-directory "${ff}" --threads "${SLURM_CPUS_PER_TASK:-8}" --work-directory "${work}/matrix"; collect_results ;;
  3) "${core}" validation/validate_openmm_engine.py --openmm-python "${openff}" --tleap-executable /opt/envs/core/bin/tleap --parmed-python "${core}" --gmx "${gmx}" --work-directory "${work}/openmm"; collect_results ;;
  4)
    "${core}" validation/validate_chitin_glycam.py --work-directory "${work}/chitin-glycam" --gmx-double "${gmx_d}"
    collect_results
    "${core}" validation/validate_chitin_glycam_cross_engines.py --work-directory "${work}/chitin-cross" --openmm-python "${openff}" --gmx-double "${gmx_d}" --lammps "${lmp}"
    collect_results
    ;;
  5) "${openff}" validation/validate_chitin_openff.py; collect_results ;;
  6)
    "${openff}" validation/validate_openff_fibrils.py
    collect_results
    "${openff}" validation/validate_openff_finite_slab.py --charge-method nagl
    collect_results
    "${openff}" validation/validate_openff_periodic_bulk.py
    collect_results
    "${openff}" validation/validate_openff_periodic_slab.py
    collect_results
    "${openff}" validation/validate_openff_shipped_outputs.py
    collect_results
    "${openff}" validation/measure_openmm_tilted_box.py
    collect_results
    "${openff}" validation/measure_openff_mixed_precision.py
    collect_results
    ;;
  7) "${openff}" validation/validate_opls_cellulose.py --openmm-python "${openff}" --gmx-double "${gmx_d}"; collect_results ;;
  8)
    "${core}" validation/validate_charmm36_engines.py --polymer cellulose --work-directory "${work}/charmm-cellulose" --openmm-python "${openff}" --lammps "${lmp}"
    collect_results
    "${core}" validation/validate_charmm36_engines.py --polymer chitin --work-directory "${work}/charmm-chitin" --openmm-python "${openff}" --lammps "${lmp}"
    collect_results
    ;;
  9)
    "${core}" validation/validate_cellulose_ii_hydroxyls.py
    collect_results
    "${core}" validation/benchmark_coordinate_matrix.py
    collect_results
    ;;
  *) echo "unknown validation task: ${task_id}" >&2; exit 2 ;;
esac
git status --short > "${output_root}/task_${task_id}/git-status.txt"
