# Periodic GLYCAM allomorph validation

This directory records matched GROMACS engineering checks for periodic bulk
cellulose I-alpha, I-beta, II, and III-I. Each system was minimized with PME
and then integrated for 1 ps in the NVT ensemble at 300 K. The systems use DP8
chains. III-I uses a 6x4x4 replication because its shorter crystallographic
`a` axis makes 4x4x4 incompatible with the validation protocol's cutoff.

Each allomorph directory contains its build report, minimization and NVT logs,
the final GRO structure, and the NVT observable table. Large trajectories and
binary GROMACS inputs are excluded.

Cellulose II requires special treatment because its diffraction dataset omits
HO2, HO3, and HO6. CelluloseForge constructs only those force-field hydrogen
coordinates using the hydroxyl directions obtained by aligning the GLYCAM
residue to crystallographic C1, C3, and C5 while preserving the source O--H
bond lengths. The origin and center chains are mapped in opposite residue directions
to preserve their antiparallel beta(1->4) connectivity. The generated O--H
distances are validated before GROMACS integration.

These short runs validate topology construction and numerical stability. They
do not establish equilibrium lattice parameters or comparative force-field
accuracy.
