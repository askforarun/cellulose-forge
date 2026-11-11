# 0002: CHARMM36 parameterization of chitin (audit)

- Status: implemented; GROMACS output enabled after the reference comparison
  and GROMACS validation passed. OpenMM and LAMMPS remain gated.
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

### Canonical patch definition

The canonical CHARMM topology confirms these values. `top_all36_carb.rtf` in
`toppar_c36_feb26.tgz` (the February 2026 release matching the pinned GROMACS
port; SHA-256
`7d4d21ebdb216a48f600fadee862dc77d369a10d5753150635f35e09c3209e08`,
downloaded from mackerell.umaryland.edu) defines

```
PRES 14bb   ! (i)1->4(i-1) equatorial at C1 and equatorial at C4
dele atom 1HO4
dele atom 2HO1
dele atom 2O1
ATOM 1C4  CC3161    0.09
ATOM 1O4  OC301    -0.36
ATOM 2C1  CC3162    0.29
BOND 1O4  2C1
```

which matches the constants in `charmm36.py`. The patch refers only to atoms
named C1, C4, O4, O1, HO1, and HO4, so it applies to any hexopyranose with
those atoms, including `BGLCNA`. (The header's "0.02" is a label; the listed
deletions and charge changes sum to zero.) `RESI BGLCNA` in the same file has
the same 30 atoms, types, charges, and the two `IMPR` entries as the GROMACS
`carb.rtp`.

This confirms the linkage definition. The topology-level gate still applies:
the generated chitin topology must be compared with an independently built
reference (CHARMM-GUI, or psfgen/ParmEd applied to this RTF) before output is
enabled.

## Terminal chemistry

Finite chains use the CHARMM free-sugar termini, as for cellulose: the
reducing end keeps O1/HO1 (beta anomeric hydroxyl, C1 unpatched) and the
non-reducing end keeps HO4 (C4/O4 unpatched). No other terminal patches are
applied. The reducing-end anomer follows the crystal (beta).

Terminal residues are not individually neutral, which is normal for CHARMM
free-sugar termini:

| Residue | Atoms | Charge (e) | Changes from BGLCNA |
|---|---|---|---|
| Internal | 27 | 0.00 | O1, HO1, HO4 deleted; C1, C4, O4 patched |
| Reducing end | 29 | -0.18 | HO4 deleted; C4, O4 patched |
| Non-reducing end | 28 | +0.18 | O1, HO1 deleted; C1 patched |

A finite chain of DP N therefore has 27N + 3 atoms and zero net charge. The
reference comparison must check these per-residue charges, not require each
residue to be neutral.

## Periodic linkage

For covalently periodic chains every residue is patched on both sides: O4 of
the last residue bonds C1 of the first across the chain-axis boundary.
Bonds, angles, proper dihedrals, and 1-4 pairs are generated from that cyclic
bond graph, exactly as for cellulose. The two acetyl impropers are not
generated from the graph: they are residue-local terms imported from the
`BGLCNA` definition, two per residue, and never cross the periodic boundary.

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
   impropers. Because the impropers never span the boundary, the check is
   that both improper terms per residue survive parsing and serialization into
   `system.xml` and that their energy matches GROMACS term by term.
5. LAMMPS: the converter currently rejects `CustomTorsionForce`. GROMACS
   function-2 impropers use E = (1/2) k (xi - xi0)^2, whereas LAMMPS
   `improper_style harmonic` uses E = K (chi - chi0)^2, so
   K_LAMMPS = k_GROMACS / 2 after converting kJ mol^-1 to kcal mol^-1
   (divide by 4.184). Omitting the factor of one half would double the
   improper energy. This must be added and validated against GROMACS before
   chitin LAMMPS output is enabled.

## Validation plan

Bulk first, then slabs, finite and periodic chains separately: charge
neutrality; topology bonds checked with minimum-image distances; improper
count (two per residue); GROMACS preprocessing without warnings, zero-step
energies, minimization, and short NVT; then term-by-term agreement with
OpenMM and LAMMPS, including the improper term.
