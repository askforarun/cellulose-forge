# CHARMM36 backend

CelluloseForge provides a CHARMM36-to-GROMACS path alongside the GLYCAM
implementation. It does not use `psfgen`: CelluloseForge reads `BGLC` from the
installed force field, applies finite or cyclic `14bb` beta(1->4) linkages, generates
the complete bonded graph, maps the crystallographic coordinates, and writes a
matched GROMACS system.

The present target supports pristine finite chains and periodic cellulose
I-alpha, I-beta, II, and III-I. Periodic chains require an even degree of
polymerization of at least four. Bonds, pairs,
angles, and proper dihedrals are generated from the patched molecular graph,
including the final-to-first glycosidic boundary.
CelluloseForge uses the pinned February 2026 CHARMM36/CGenFF 5.0 GROMACS port.
Install it once from the official MacKerell Laboratory distribution:

```bash
cellulose-forge forcefield install charmm36
```

The installer downloads the exact supported archive over HTTPS, verifies its
SHA-256 checksum, and places it in the user-local CelluloseForge data directory.
An upstream update therefore cannot silently alter an existing workflow. Use
`cellulose-forge forcefield path charmm36` to display the installed location.

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape bulk \
  --cells 4 4 6 \
  --chain-boundary periodic \
  --force-field charmm36 \
  --output ibeta_charmm36
```

The output directory contains crystallographic PDB coordinates, matched GRO
coordinates, the cellulose ITP, `system.top`, a local copy of the installed
force-field directory, and a machine-readable build report. The optional
`--charmm-template` argument retains the earlier reference-driven route for
regression studies; normal periodic builds do not require it.

The first 4x4x6 I-beta engineering test contained 32 periodic DP12 chains and
8,064 atoms. GROMACS 2024.2 accepted the system with CHARMM 1--4 scaling of
`fudgeLJ = 1` and `fudgeQQ = 1`. PME steepest-descent minimization converged in
1,115 steps to a maximum force of 923.17 kJ mol-1 nm-1, followed by a completed
1-ps NVT integration without LINCS or non-finite-energy failures.

The minimization input produced an excluded-atom-distance note before
relaxation; the subsequent NVT preprocessing reported a maximum excluded
distance of 0.406 nm. A separately generated periodic DP8 system also passed
`gmx grompp` and converged in 574 minimization steps to a maximum force of
909.10 kJ mol-1 nm-1. These runs are topology and runtime smoke tests, not
equilibration or validation of material properties.

`validation/validate_charmm36_reference.py` independently compares the supplied
ITP repeat with the `BGLC` residue after cyclic application of the CHARMM36
`14bb` beta(1->4) patch. For the audited reference it confirms removal of O1,
HO1, and HO4; the documented C1, C4, and O4 atom-type/charge changes; a neutral
21-atom repeat; and semantic equality of all 2,304 source and regenerated
bond, pair, angle, and dihedral records at DP12. The compact result, including
SHA-256 identifiers for both external inputs, is recorded in
`validation/results/charmm36_reference.json`. In addition,
`validation/audit_charmm36_coefficients.py` checks every parameter type actually
used by the periodic chain against the carbohydrate PRM: 10 nonbonded atom
types, 16 bond types, 38 angle/Urey--Bradley types, and 68 proper-dihedral type
tuples. All pass the CHARMM-to-GROMACS energy, distance, harmonic-factor, and
Lennard--Jones conversions; the result is recorded in
`validation/results/charmm36_coefficients.json`. External force-field files
remain unvendored. For finite chains, the reducing end retains `O1/HO1`, the
nonreducing end retains `HO4`, and intervening linkages receive the same `14bb`
atom changes. Terminal coordinates are aligned from a documented finite
cellulose reference; topology parameters come from the installed CHARMM36 release.

Periodic shape validation also covers an I-beta (110) slab and the published
36-chain I-beta fibril at DP8. Both passed `gmx grompp` with PME and converged
below 1000 kJ mol-1 nm-1 during steepest-descent minimization. The 1,344-atom
slab converged in 348 steps (Fmax 944.65 kJ mol-1 nm-1), and the 6,048-atom
fibril converged in 647 steps (Fmax 927.06 kJ mol-1 nm-1). Preprocessing issued
the already documented initial 0.952-nm excluded-distance note for each system;
neither run required `-maxwarn`. Compact evidence is stored in
`validation/results/charmm36_ibeta_shapes.json`.

Periodic bulk support now covers I-alpha, I-beta, II, and III-I. DP8 systems
for all four allomorphs pass `gmx grompp`, PME minimization, and 1-ps NVT
integration without LINCS or non-finite-energy failures. Cellulose II maps its
antiparallel chain class in the opposite residue direction and constructs the
HO2, HO3, and HO6 coordinates missing from its diffraction dataset by aligning
an I-beta O--H reference geometry on C1/C3/C5. These generated coordinates are
reported explicitly in `build_report.json`. Compact matrix results are stored
in `validation/results/charmm36_allomorph_matrix.json`.

The literature-defined 36-chain periodic fibrils are also exercised for every
allomorph for which that preset exists: I-alpha, I-beta, and II. Each DP8
fibril contains 6,048 atoms, and all three pass preprocessing, minimization,
and 1-ps NVT smoke tests. Cellulose II residue order is derived from the
periodic crystallographic O4--C1 connectivity of each chain rather than from
chain-number parity, which is not valid for the fibril selection. Compact
results are stored in `validation/results/charmm36_fibril_matrix.json`.

Finite DP6 versions of the same I-alpha, I-beta, and II 36-chain fibrils also
pass `gmx grompp` and converge below 1000 kJ mol-1 nm-1. Residue order is
derived independently for every open chain from its O4--C1 connectivity before
the reducing `O1/HO1` and nonreducing `HO4` coordinates are added. This avoids
placing a terminal group on an internal crystallographic residue in staggered
fibrils. Compact results are stored in
`validation/results/charmm36_finite_fibril_matrix.json`.

After tighter PME minimization to below 100 kJ mol-1 nm-1, the finite I-beta
(110) slab and all three finite fibrils complete 1 ps NVT integrations without
LINCS or non-finite-energy warnings. The short runs are integration smoke tests,
not equilibrated trajectories; their temperatures and energies are recorded
without interpreting them as material properties in
`validation/results/charmm36_finite_pme_nvt.json`.

`validation/validate_finite_charmm36_reference.py` provides an independent
topology check against the DP15 cellulose example processed by CHARMM-GUI v1.7.
After normalization by residue number and atom name, all 318 atom identities,
types, and charges match, as do all 2,883 bond, pair, angle, and dihedral
records. File hashes and counts are stored in
`validation/results/charmm36_finite_reference.json`; the external reference is
not vendored as a runtime dependency.
