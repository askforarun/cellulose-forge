# Experimental periodic OPLS-CM5 allomorph validation

This directory records matched GROMACS engineering checks for periodic bulk
cellulose I-alpha, I-beta, II, and III-I using the experimental OPLS-CM5
adapter. Each system was minimized with PME and integrated for 1 ps in the NVT
ensemble at 300 K. All systems use DP8 chains. III-I uses a 6x4x4 replication
because its shorter crystallographic `a` axis makes 4x4x4 incompatible with the
validation cutoff.

Each allomorph directory contains its sanitized build report, minimization and
NVT logs, final GRO structure, and NVT observable table. Large trajectories and
binary GROMACS inputs are excluded. The source `CNC.itp` and `solute.gro` are
also excluded because they are external inputs without a stated redistribution
license.

Cellulose II uses opposite topology-to-crystal residue mappings for its
antiparallel origin and center chains. Its diffraction data omit HO2, HO3, and
HO6, so the adapter constructs those force-field coordinates from the aligned
OPLS reference while preserving the source O--H bond lengths.

These runs establish preprocessing, topology, and short-integration stability.
They do not establish equilibrium properties or independently validate the
physical transferability of an I-beta-derived OPLS-CM5 topology to the other
allomorphs.
