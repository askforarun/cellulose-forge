# Examples

Validated example systems will be added with the first cellulose I-beta vertical prototype. Each example will record its build command, force-field version, GROMACS version, and output checksums.

The [`i_beta_unit_cell`](i_beta_unit_cell/) example is the first geometry-only vertical slice. It records one complete crystallographic unit cell before molecular connectivity or force-field assignment.

The [`i_alpha_unit_cell`](i_alpha_unit_cell/) example records the complete
triclinic P1 asymmetric unit, which is also the full I-alpha unit cell.

The [`iii_i_unit_cell`](iii_i_unit_cell/) example applies P2_1 symmetry to the
single fully hydrogenated III-I asymmetric residue.

The [`ii_unit_cell`](ii_unit_cell/) example records the antiparallel origin and
center chains. Its experimentally unresolved hydroxyl hydrogens remain absent.

The [`fibrils/ibeta_36chain_dp12.pdb`](fibrils/ibeta_36chain_dp12.pdb)
example is the literature-defined I-beta elementary fibril: 36 chains, six
cellobiose repeats per chain (DP 12), periodic chain-boundary metadata, and 30
angstrom of total padding along each transverse direction.

The [`md_inputs`](md_inputs/) folder holds small simulation-ready systems
for GROMACS, OpenMM, LAMMPS, and Amber; see its README for run commands.
