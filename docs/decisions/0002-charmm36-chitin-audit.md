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
these atoms have identical types and charges in `BGLC` and `BGLCNA`, and the
patched GlcNAc repeat has 27 atoms and zero net charge, matching the 27-atom
crystallographic residue.

This is necessary but not sufficient evidence. The GROMACS distribution does
not contain the CHARMM patch definitions; the `14bb` values in CelluloseForge
were transcribed and then confirmed for cellulose against an independent
CHARMM-GUI topology (atom types, charges, and every bonded term of a DP15
chain). Applying them to GlcNAc is therefore an assumption until the same
check is made for chitin: generate a finite beta(1->4)-linked GlcNAc oligomer
with CHARMM-GUI Glycan Reader (or from the CHARMM `top_all36_carb.rtf`
`14bb` patch) and require exact agreement of atom types, charges, bonds,
angles, proper and improper dihedrals, and 1-4 pairs, including both
termini. Implementation may start, but chitin CHARMM36 output must not be
enabled until this comparison passes.

## Terminal chemistry

Finite chains use the CHARMM free-sugar termini, as for cellulose: the
reducing end keeps O1/HO1 (beta anomeric hydroxyl, C1 unpatched) and the
non-reducing end keeps HO4 (C4/O4 unpatched). No other terminal patches are
applied. The reducing-end anomer follows the crystal (beta).

## Periodic linkage

For covalently periodic chains every residue is patched on both sides: O4 of
the last residue bonds C1 of the first across the chain-axis boundary, and
angles, dihedrals, 1-4 pairs, and impropers are generated from that cyclic
bond graph, exactly as for cellulose. Impropers are local to the acetyl group
and never span the boundary.

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
