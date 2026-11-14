# 0003: GLYCAM06j-1 parameterization of chitin (audit)

- Status: accepted for GROMACS, OpenMM, LAMMPS, and finite-chain Amber output
  after validation
- Branch: `feature/chitin-glycam`
- Base release: `v1.2.0` (`652ff76`)

## Scope

This decision covers alpha- and beta-chitin built from the crystallographic
definitions already shipped by CelluloseForge.  The intended first milestone
is finite and covalently periodic GLYCAM06j-1 output for bulk crystals and
chain-parallel slabs.  Chitin fibrils remain outside scope because no
authoritative fibril construction is implemented.

## GLYCAM residue definitions

AmberTools 25, loaded through `leaprc.GLYCAM_06j-1`, contains the
N-acetyl-beta-D-glucosamine residue family required for beta(1->4) chitin:

| Position in a finite chain | GLYCAM residue | Atoms | Residue charge (e) |
|---|---:|---:|---:|
| Reducing-end cap | `ROH` | 2 | -0.194 |
| Internal beta(1->4) GlcNAc | `4YB` | 27 | 0.000 |
| Non-reducing GlcNAc | `0YB` | 28 | +0.194 |

A finite DP6 reference generated independently with `tleap` therefore uses
`ROH 4YB 4YB 4YB 4YB 4YB 0YB`, contains 165 atoms, and is neutral. Its
topology contains 170 bonds, 313 angles, 539 dihedral records (including 12
impropers), and 456 unique 1--4 pairs.

A covalently periodic DP6 reference uses six neutral `4YB` residues followed
by an explicit `O4(DP)--C1(1)` closure.  It contains 162 atoms, 168 bonds,
312 angles, 540 dihedral records (including 12 impropers), and 456 unique
1--4 pairs. The total charge is zero. These counts are audit observations, not hard-coded production
expectations; the production comparison will inspect interaction identities
and parameters as well as counts.

Most GlcNAc proper torsions retain GLYCAM's unity electrostatic and
Lennard-Jones 1--4 scaling.  The amide `H2N--N2--C2N--O2N` terms instead carry
the Amber scaling stored by GLYCAM06j-1 (`SCEE=1.2`, `SCNB=2.0`).  The
zero-valued scale entries belong to the two acetyl impropers per linked
residue and do not generate ordinary 1--4 pairs.  Consequently, chitin output
must preserve per-dihedral scaling from the reference topology rather than
assuming the uniform unity scaling used by cellulose.

## Atom-name mapping

The crystallographic heavy-atom coordinates are retained.  GLYCAM atom names
map to the chitin crystal labels as follows:

| GLYCAM | Crystal |
|---|---|
| `H2N` | `HN2` |
| `C2N` | `C7` |
| `O2N` | `O7` |
| `CME` | `C8` |
| `H3M`, `H2M`, `H1M` | `H81`, `H82`, `H83` after geometry-based assignment |
| `H3O`, `H6O` | `HO3`, `HO6` |

The methyl-hydrogen correspondence must be assigned by the minimum rigid
geometry match after aligning the acetyl-group atoms. The three methyl
hydrogens have identical atom types and zero charges, so their individual
labels do not change the energy. Every other GLYCAM name
is identical to its crystallographic label.  Finite-only `O1`, `HO1`, and
`H4O` coordinates are transferred from the corresponding `tleap` reference by
the same proper rigid transformations used by the cellulose mapper.

## Validation outcome

The GROMACS, OpenMM, LAMMPS, and finite-chain Amber output was enabled after the following
checks passed:

1. A finite DP6 chain matches the independent `tleap`/ParmEd reference in
   atom names, types, charges, masses, bonds, angles, proper and improper
   dihedrals, 1--4 pairs, exclusions, and both terminal compositions.
2. A periodic DP6 chain matches the independently closed `tleap` reference,
   including every bond, angle, dihedral, pair, and exclusion crossing the
   final-to-first residue boundary.
3. Alpha- and beta-chitin bulk and `(100)` slab cases pass coordinate-to-
   topology bond checks, periodic-image clash checks, `gmx grompp -maxwarn 0`,
   minimization, and a one-ps NVT smoke test.
4. Finite slabs use at least 12 angstrom of longitudinal padding for the
   validation cases and contain no chain-end contact with their periodic
   images.
5. Native Amber output is enabled only for finite chains. Periodic chitin is
   emitted through the validated GROMACS topology rather than as a periodic
   `prmtop` that Amber engines could misinterpret.
6. A finite beta-chitin `(100)` slab on identical coordinates gives matching
   bonded and 1--4 energies in double-precision GROMACS and Amber `sander`.
   The largest corrected absolute difference is 0.000124 kcal/mol. The raw
   1--4 electrostatic difference follows the ratio between the GROMACS and
   Amber Coulomb constants (332.0637 and 332.0522 kcal A mol-1 e-2).
7. OpenMM stores the mixed amide pairs as pair-specific nonbonded exceptions.
   LAMMPS excludes each affected ordinary pair and restores its exact scaled
   Lennard-Jones plus Coulomb energy through `bond_style table`. This tabulated
   1--4 correction is not chemical connectivity. Alpha/beta finite/periodic
   slabs and periodic bulk crystals agree with double-precision GROMACS within
   8.4e-6 for bonded and
   explicit-pair categories and 1.6e-4 for the combined nonbonded energy.

The compact six-case GROMACS matrix, finite-chain Amber smoke test, and
term-by-term Amber/GROMACS energy comparison are in
`validation/results/chitin_glycam_engines.json`. The GROMACS conversion is
created directly from the validated Amber interaction graph; the only graph-
preserving modification is replacement of the six amide 1--4 rows per DP6
chain with explicit function-2 parameters carrying SCEE=1.2 and SCNB=2.0.

The short contacts in the ideal coordinates written by `tleap` are irrelevant
to the crystal mapping: `tleap` supplies the reference topology and terminal
geometry, while CelluloseForge replaces the internal coordinates with the
experimental chitin coordinates before simulation validation.
