# Matched force-field comparison evidence

These files record the September 20, 2026 engineering validation of the same
periodic cellulose I-beta 4x4x4 bulk crystal with GLYCAM06j-1 and the
experimental OPLS-CM5 adapter.

- `summary.json` contains the comparison and interpretation limits.
- Each force-field directory contains its CelluloseForge build report,
  minimization/NVT/NPT logs, compact observable tables, and the last coordinate
  frame as GRO.
- OPLS completed 100 ps NPT; its late analysis window is 50--100 ps.
- GLYCAM was intentionally stopped at 38.48 ps after stability was established;
  its late analysis window is 20--38.48 ps.

The XTC trajectories, EDR energy files, TPR run inputs, and checkpoints are not
versioned. They are reproducible with `validation/compare_force_fields.py` and
should be archived separately if needed for a publication dataset.

These results demonstrate file interoperability, correct force-field-specific
1--4 scaling, and stable GROMACS integration. They do not establish converged
thermodynamic properties or relative force-field accuracy. Absolute potential
energies must not be compared between GLYCAM and OPLS-CM5.
