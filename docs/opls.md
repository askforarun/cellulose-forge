# Experimental OPLS-AA and OPLS-AA-SEI support

CelluloseForge can parameterize cellulose with two experimental all-atom
carbohydrate force fields:

| `--force-field` | Parameters |
|---|---|
| `opls-aa` | OPLS-AA for carbohydrates (Damm, Frontera, Tirado-Rives and Jorgensen, *J. Comput. Chem.* 1997, 18, 1955), with torsions from the "OPLS-AA" column of Kony et al. Table 3 |
| `opls-aa-sei` | OPLS-AA-SEI (Kony, Damm, Stoll and van Gunsteren, *J. Comput. Chem.* 2002, 23, 1416; doi:10.1002/jcc.10139): the "OPLS-AA-SEI" torsions of Table 3 plus scaled intramolecular electrostatics |

No optional environment, AmberTools or downloaded force-field file is
needed: the parameters are written explicitly into every output.

```bash
cellulose-forge build \
  --allomorph I-beta --shape slab --face 1 1 0 \
  --layers 3 --chains-per-layer 4 --degree-polymerization 6 \
  --chain-boundary periodic --vacuum 30 \
  --force-field opls-aa-sei \
  --md-engine gromacs --md-engine openmm --md-engine lammps \
  --output opls_sei_ibeta_slab
```

OpenMM output needs a Python interpreter with OpenMM (`--openmm-python`, or
the current interpreter). GROMACS and LAMMPS output need only NumPy.

## Scope

Finite and covalently periodic slabs and fibrils, and covalently periodic bulk
crystals, of cellulose I-alpha, I-beta, II and III-I (fibrils: the
literature 36-chain fibril for I-alpha, I-beta and II, and the I-beta 18- and
24-chain presets and custom profiles). GROMACS, OpenMM and LAMMPS output.
Chitin, mannan, unit cells, finite bulk crystals, Amber output and position
restraints are rejected.

## Parameter sources

- **Atom types, charges, Lennard-Jones, bonds and angles** are the OPLS-AA
  hexopyranose parameters of Damm et al. (1997), as distributed in GROMACS
  `oplsaa.ff` (types `opls_140`, 157, 158, 169, 170, 183, 185--188,
  193--196). Internal residues are exactly neutral; the reducing-end
  hemiacetal residue carries -0.2 e and the non-reducing-end residue +0.2 e,
  so each finite chain is neutral.
- **Combination rule and 1-4 scaling** follow OPLS: geometric means for both
  sigma and epsilon, and 1-4 Coulomb and Lennard-Jones scaled by 0.5.
- **Torsions** are Kony et al. Table 3 Fourier coefficients V1--V3 (kcal/mol).
  GROMACS receives them as Ryckaert--Bellemans coefficients, OpenMM as an
  `RBTorsionForce`, and LAMMPS directly as `dihedral_style opls` K1--K3.
- **OPLS-AA-SEI electrostatic scaling**: 1-5 Coulomb interactions within each
  glucose residue are divided by 1.26, and the 1-6 interactions between the
  hydroxymethyl group and a vicinal C4 hydroxyl (HO6--O4 and O6--HO4) by
  1.22. Lennard-Jones remains at full strength for these pairs.

## Interpretations

Kony et al. parameterized monosaccharides. Applying the parameters to a
beta-(1->4) polymer needs three choices, all recorded in `build_report.json`:

1. **OS-C-C-OS** (O5-C5-C4-O4 when O4 is a glycosidic ether) does not occur in
   a free hexopyranose and is absent from Table 3. Both variants use the
   OPLS-AA polyether term (V1 = -0.55 kcal/mol, GROMACS `oplsaa.ff`).
2. **C-C-C-OS_ace** is assigned to C3-C2-C1-O with O the exocyclic acetal
   (glycosidic) oxygen on C1; C3-C2-C1-O5 uses C-C-C-OS.
3. **SEI scaling is applied within each residue only**, as published. The 1-6
   term applies only where O4 is a hydroxyl, that is, at the non-reducing end
   of a finite chain; in internal residues O4 is a glycosidic ether. A
   covalently periodic chain therefore has 47 scaled 1-5 pairs per residue and
   no 1-6 pairs.

The Table 3 OPLS-AA column agrees with the hexopyranose dihedral types in
GROMACS `oplsaa.ff` to the three printed decimals, which the validation
checks.

## Engine representations

**GROMACS.** `cellulose_opls.itp` lists every bond, angle and dihedral with
explicit parameters and includes `oplsaa.ff/forcefield.itp` only for the atom
types and defaults. OPLS-AA-SEI pairs are `[ pairs ]` function 2 entries
(fudgeQQ, both charges, geometric sigma and full epsilon) plus matching
`[ exclusions ]`.

**OpenMM.** OpenMM's GROMACS reader rejects function-2 pairs, so the `System`
is built directly. Coulomb uses a PME `NonbondedForce`; Lennard-Jones uses a
`CustomNonbondedForce` with the geometric combination rule and a plain 1.0 nm
cutoff. Every 1-4 and SEI pair is a `NonbondedForce` exception with its own
charge product, sigma and epsilon and is excluded from the custom force.

The build-time energy check normally runs on OpenMM's double-precision
Reference platform. For boxes whose third vector is tilted (cellulose I-alpha
bulk crystals) it runs on the CPU platform instead. In local measurements with
OpenMM 8.6.1 (`validation/measure_openmm_tilted_box.py`,
`validation/results/openmm_tilted_box.json`), the Reference platform's
Lennard-Jones energy for two I-alpha bulk boxes differed from an exact
nearest-image sum by 0.14 and 0.40 %, unchanged when coordinates were wrapped
into the cell, while the CPU platform agreed within 10^-6 and untilted
I-beta, II and III-I boxes agreed on Reference within 10^-13. The serialized
system is unaffected; the platform used is recorded as `validation_platform`.

**LAMMPS.** The data file is written directly (`pair_modify mix geometric`,
`special_bonds lj/coul 0.0 0.0 0.5`, `dihedral_style opls`). LAMMPS has no
per-pair 1-5 scaling, and `bond_style special`, which exists for this purpose,
refuses long-range Coulomb. Each OPLS-AA-SEI pair is therefore written as a
tabulated pseudo-bond (`sei_pairs.table`, 0.001 A spacing) whose energy is the
pair's full Lennard-Jones plus its scaled Coulomb term. The pseudo-bond
excludes the pair from the pair style, which removes the pair's full Ewald
interaction exactly, and `special_bonds ... angle yes dihedral yes` stops the
pseudo-bonds from creating false 1-3 or 1-4 exclusions. Consequences:

- `ebond` includes the SEI pair energies;
- LAMMPS prints two expected warnings: table forces flagged at inflection
  points, and "Bond/angle/dihedral extent > half of periodic box length", a
  heuristic that treats a ~5 A pseudo-bond as part of a dihedral;
- a pair that ever moved outside 1--8 A would stop LAMMPS with an error (the
  largest intra-residue 1-5/1-6 distance is about 6 A).

## Nonbonded policy

Real-space cutoff 1.0 nm (PME/PPPM) in every engine. As for GLYCAM06 and
Sage, the long-range dispersion correction (OPLS-AA was fitted with it) is
used for bulk crystals and disabled for slabs and fibrils, whose transverse
vacuum violates its homogeneous-fluid assumption: `DispCorr` in GROMACS, the
custom force's long-range correction in OpenMM, and `pair_modify tail` in
LAMMPS. The generated GROMACS inputs use the default potential-shift modifier;
OpenMM and LAMMPS use a plain cutoff. Like every force-field build, the
periodic box must exceed 21 A (twice the 1.0 nm cutoff plus a 0.5 A pair-list
buffer on each side) in every direction.

## Validation boundary

`validation/validate_opls_cellulose.py` checks both variants for finite and
periodic slabs and 36-chain fibrils and periodic bulk crystals of every
supported allomorph:

- **transcription**: the Table 3 OPLS-AA torsions and every explicit bond and
  angle parameter against the installed GROMACS `oplsaa.ff`;
- **independent recomputation**: torsion, 1-4 and SEI-pair energies recomputed
  in NumPy from the written topology against double-precision GROMACS, and the
  OPLS-AA-SEI minus OPLS-AA nonbonded energy against the analytic effect of the
  scaling;
- **cross-engine agreement** on identical, unrounded coordinates with
  unshifted cutoffs and tight Ewald settings, against double-precision
  GROMACS: OpenMM bond, angle, torsion, 1-4/SEI exception, Lennard-Jones and
  total energies; LAMMPS harmonic-bond, SEI-table, angle, torsion,
  pair-style Lennard-Jones, Coulomb and total energies (the SEI tables are
  compared with the independently recomputed SEI pair energy);
- **runtime**: `gmx grompp -maxwarn 0`, minimization and 1 ps NVT in GROMACS;
  minimization and 1 ps of dynamics with the generated `run_openmm.py`; the
  generated LAMMPS input.

These checks establish that the topologies reproduce the published parameters
and that the three engines evaluate the same energy function. They do **not**
establish the accuracy of OPLS-AA or OPLS-AA-SEI for crystalline cellulose.
