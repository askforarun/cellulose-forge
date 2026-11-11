# 0002: CHARMM36 parameterization of chitin (audit)

- Status: proposed
- Branch: `feature/chitin-charmm36`
- Force field: CHARMM36 February 2026 GROMACS release
  (`charmm36-feb2026_cgenff-5.0.ff`), as already pinned for cellulose.

## Residue

The release defines N-acetyl-beta-D-glucosamine as `[ BGLCNA ]` in `carb.rtp`
("2-acetyl-2-deoxy-beta-D-glucosamine"): 30 atoms, net charge 0. It differs
from `BGLC` only at C2, where O2/HO2 are replaced by an N-acetyl group, and it
adds two improper dihedrals:

- `C CT N O` (carbonyl planarity)
- `N C C2 HN` (amide hydrogen planarity)

Both impropers have function-2 parameters in `ffbonded.itp`
(`CC2O1 CC331 NC2D1 OC2D1`, k = 1004.16 kJ mol^-1 rad^-2, and
`NC2D1 CC2O1 CC3161 HCP1`, k = 167.36 kJ mol^-1 rad^-2). The amide bonds and
angles are also present.

## beta(1->4) linkage

CelluloseForge applies the CHARMM `14bb` patch in code (`charmm36.py`):
delete O1, HO1, and HO4, and set C1 to CC3162/0.29 e, C4 to CC3161/0.09 e, and
O4 to OC301/-0.36 e. The patch only touches C1, C4, O4, O1, HO1, and HO4, and
these atoms have identical types and charges in `BGLC` and `BGLCNA`. The same
patch therefore applies to chitin, and the patched GlcNAc repeat has 27 atoms
and zero net charge, matching the 27-atom crystallographic residue.

## Atom-name mapping

| CelluloseForge (crystal) | CHARMM36 `BGLCNA` |
|---|---|
| N2, HN2 | N, HN |
| C7, O7 | C, O |
| C8, H81, H82, H83 | CT, HT1, HT2, HT3 |
| C1-C6, O3-O6, H1-H5, H61, H62, HO3, HO6 | same names |

Finite chains add the CHARMM terminal O1/HO1 (reducing end) and HO4
(non-reducing end); the existing rigid ring-alignment used for cellulose
termini can be reused because the ring atoms are shared.

## Required code changes

1. Generalize the residue reader from `BGLC` to a named residue, and carry
   `[ impropers ]` rows from the RTP into the topology as function-2 dihedrals.
2. Map the chitin atom names above; keep cellulose behaviour unchanged (the
   cellulose hash test and existing GROMACS tests must pass unmodified).
3. GROMACS: impropers are resolved by `grompp` from `ffbonded.itp`.
4. OpenMM: `GromacsTopFile` creates a `CustomTorsionForce` for function-2
   impropers. The OpenMM writer already enables periodic boundary conditions
   on every bonded force that supports them; this must be confirmed for the
   improper force across the periodic boundary.
5. LAMMPS: the converter currently rejects `CustomTorsionForce`. CHARMM
   harmonic impropers map to `improper_style harmonic` (K in kcal mol^-1
   rad^-2 with LAMMPS's E = K(chi - chi0)^2 convention); this must be added
   before chitin LAMMPS output is enabled.

## Validation plan

Bulk first, then slabs, finite and periodic chains separately: charge
neutrality; topology bonds checked with minimum-image distances; improper
count (two per residue); GROMACS preprocessing without warnings, zero-step
energies, minimization, and short NVT; then term-by-term agreement with
OpenMM and LAMMPS, including the improper term.
