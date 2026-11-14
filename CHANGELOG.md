# Changelog

All notable changes to CelluloseForge are documented in this file.

## Unreleased

### Added

- GLYCAM06j-1 topology for finite and covalently periodic alpha- and
  beta-chitin slabs, with GROMACS, OpenMM and LAMMPS output and native Amber
  output for finite chains. Mixed GlcNAc 1-4 scaling is represented explicitly
  in GROMACS, as native OpenMM exceptions, and as pair-specific tabulated 1-4
  corrections in LAMMPS. Alpha/beta finite/periodic cross-engine energies are
  validated term by term against double-precision GROMACS.

### Changed

- Corrected the modelled cellulose-II HO2 and HO6 orientations to follow the
  chain-specific neutron-diffraction network in Langan, Nishiyama and Chanzy
  (JACS 1999, Table 3): O2o-H...O6o, O6o-H...O6c, O2c-H...O2o and
  O6c-H...O2c. GLYCAM06j-1, CHARMM36, OpenFF Sage and OPLS-AA/SEI now share
  the same explicit acceptor assignments. Heavy-atom crystallographic
  coordinates and bonded topologies are unchanged.
- Native Amber smoke-test inputs now place `vdwmeth=0` in Amber's `&ewald`
  namelist. Version 1.2.0 incorrectly placed it in `&cntrl`, which AmberTools
  `sander` rejects. Native Amber output remains finite-chain-only and applies
  to slabs and fibrils, where the isotropic dispersion correction is disabled.

## 1.2.0 - 2026-11-14

### Added

- Experimental OPLS-AA (Damm et al. 1997) and OPLS-AA-SEI (Kony et al. 2002)
  for finite and periodic cellulose slabs and fibrils and periodic bulk
  crystals, with GROMACS, OpenMM and LAMMPS output (`--force-field opls-aa`
  or `opls-aa-sei`). Parameters are written explicitly; OPLS-AA-SEI scales
  1-5 electrostatics within each residue by 1/1.26 and hydroxymethyl/C4
  hydroxyl 1-6 pairs by 1/1.22. OpenMM systems are built directly; LAMMPS SEI
  pairs are tabulated pseudo-bonds. Validated across the three engines in
  `validation/validate_opls_cellulose.py`.
- Optional experimental OpenFF Sage 2.3 output for finite and periodic
  cellulose slabs and periodic bulk crystals through GROMACS and OpenMM.
- Direct OpenFF glucan graph construction without AmberTools, NAGL charges,
  and build-time stereochemistry and periodic charge-equivalence checks.
- Geometry-specific Sage dispersion policy: enabled for bulk crystals and
  disabled for vacuum slabs and fibrils.
- Experimental OpenFF Sage fibrils: the literature 36-chain fibril for
  I-alpha, I-beta and II and the I-beta 18-chain, 24-chain and custom row
  profiles, with finite or covalently periodic chains, validated in
  `validation/validate_openff_fibrils.py`.
- Every force-field build now writes starting run inputs for each selected
  engine: `minim.mdp` and a 1 ps `nvt.mdp` for GROMACS, and `min.in` and a
  1 ps `nvt.in` for Amber. OpenFF OpenMM output now uses the same configurable
  `run_openmm.py` and `state.xml` as GLYCAM and CHARMM36. Generated inputs are
  smoke-test protocols, not production settings.

### Changed

- OpenMM build-time energy checks for GLYCAM06j-1 and CHARMM36 now use the CPU
  platform for boxes whose third vector is tilted (cellulose I-alpha bulk
  crystals), as the OPLS builder already did; OpenMM 8.6.1 Reference
  nonbonded energies disagreed with exact sums for such boxes. Only the
  energies recorded in `build_report.json` change; serialized systems are
  unchanged. The platform used is recorded as `validation_platform`.
- One long-range dispersion policy for every force field and engine: GLYCAM06
  and OpenFF Sage use it only for bulk crystals and unit cells; CHARMM36 never
  uses it; slabs and fibrils never use it. **CHARMM36 OpenMM bulk-crystal and
  unit-cell output previously enabled it** and now matches the LAMMPS output,
  the cross-engine validators, and standard CHARMM36 practice. Amber inputs set
  `vdwmeth` accordingly.
- Every force-field build now rejects periodic boxes no wider than twice the
  nonbonded cutoff plus a 0.5 A pair-list buffer on each side (21 A for a
  1.0 nm cutoff, 19 A for Sage) in any perpendicular direction, and warns when a slab or
  fibril has less transverse vacuum than the longest interaction range.
  The 0.5 A margin is a conservative CelluloseForge policy: the GROMACS Verlet
  buffer depends on the MD settings, and with the generated inputs a 1.0 nm
  cutoff was rejected by `grompp` on a 20.62 A box (above the previous 20 A
  threshold) and accepted on 21.04 A. This tightens the rule for every force
  field.

- Cellulose II hydroxyl hydrogens (HO2, HO3, HO6), which are absent from the
  diffraction data, are now placed on an exact tetrahedral cone (0.96 A O-H,
  108.5 degree C-O-H) for GLYCAM and OpenFF output, class by class (HO3, HO2,
  HO6) so equivalent hydroxyls are treated identically and the result does not
  depend on atom order. HO3 takes the intrachain O3-H...O5 bond; the others take
  the clash-free torsion with the best geometric hydrogen bond (linearity and
  distance), otherwise the reference torsion. Placement uses a periodic cell
  list and scales linearly (approximately 0.76 ms per hydroxyl; 4.36 s for 34,560 atoms). This
  replaces copying the O-H direction from an aligned reference residue, which
  gave C6-O6-H angles near 128 degrees. GLYCAM cellulose II coordinates
  therefore differ from earlier development builds; topologies are unchanged.
  The recorded GLYCAM II entries in `validation/results/slab_allomorph_matrix.json`
  predate this change and are regenerated with the release validation set.

- Coordinate-only alpha- and beta-chitin unit cells and periodic bulk crystals
  with `--polymer chitin`. Definitions are converted from the refined structures
  of Sikorski, Hori and Wada (2009) and Nishiyama, Noishiki and Wada (2011) by
  `tools/convert_chitin_cif.py`; `--o6-conformation` selects the major or minor
  alpha-chitin hydroxymethyl position. Chirality is checked at all five ring
  carbons. Fibrils are not yet available for chitin.
- Chain-parallel alpha- and beta-chitin slabs for any `(h k 0)` face with finite
  or periodic chains, built by the general reciprocal-lattice slab path. Slab
  constructors now carry the definition's residue name, so chitin slabs are
  written as `NAG`.
- CHARMM36 chitin topologies for GROMACS (`--polymer chitin --force-field
  charmm36`): the `BGLCNA` residue with the `14bb` linkage, residue-local acetyl
  impropers, free-sugar termini for finite chains, and covalently periodic
  chains. Output is `chitin_charmm36.itp`; the build report records the
  polymer, the CHARMM36 residue, and the validation scope. The topology matches
  a psfgen reference built from `top_all36_carb.rtf`, and alpha and beta bulk
  and slab systems pass GROMACS preprocessing, minimization, and a 1 ps NVT
  smoke test.
- OpenMM and LAMMPS output for CHARMM36 chitin. The LAMMPS converter now
  writes GROMACS function-2 impropers as `improper_style harmonic` with
  K = k/2, and carries CHARMM 1-4 pairs whose dihedral has zero force constant
  on a zero-force-constant dihedral. Energies match GROMACS term by term for
  alpha and beta, finite and periodic slabs
  (`validation/validate_charmm36_engines.py`).
- Cellulose CHARMM36 LAMMPS audit: the same term-by-term comparison
  (`--polymer cellulose`) for I-beta, I-alpha, II and III-I slabs shows that
  the CHARMM Coulomb constant and the two charmmfsw Lennard-Jones energy
  conventions explain 85-95% of the earlier 1.8-2.5e-4 LAMMPS/GROMACS
  total-potential gap; corrected totals agree within 2.7e-5. Validation
  results record the commit, tree status and engine versions.
- A regression test that records hashes of fourteen coordinate-only cellulose
  builds, so cellulose output cannot change unnoticed.

### Changed

- Emit GLYCAM torsions with the standard LAMMPS `charmm` dihedral style and
  zero dihedral 1--4 weights, preserving the energy while avoiding the
  optional package required by `dihedral_style fourier`.

### Added

- Added experimental `--md-engine lammps` output for finite and periodic
  GLYCAM06j-1 and CHARMM36 systems. The exact converter writes `system.data`
  and `in.cellulose`, preserves CHARMM Urey--Bradley and pair-specific 1--4
  terms, and rejects unsupported force terms or non-equivalent semantics.

- Added repeatable `--md-engine gromacs|amber` selection and recorded the
  resolved engines in each build report. `--topology-format` remains as a
  deprecated compatibility alias.
- Added `--md-engine openmm` for GLYCAM06j-1 and CHARMM36 systems with finite
  or covalently periodic chains. Builds write a serialized OpenMM `system.xml`,
  `state.xml`, `topology.cif`, and a standalone `run_openmm.py`. Periodic chains
  enable periodic boundary conditions for every bonded force and for
  nonbonded 1-4 exceptions. Position restraints are included as a switchable
  force, and the runner keeps the chain axis fixed under pressure coupling.
- Added I-beta 18-chain (`2,3,4,4,3,2`) and 24-chain
  (`3,4,5,5,4,3`) fibril presets.
- Added `--fibril-profile` for deterministic custom I-beta row profiles within
  the validated 36-chain lattice envelope.
- Added complete multi-chain, finite-chain GLYCAM Amber `prmtop`/`inpcrd`
  output, with a ParmEd reload check recorded in the build report. The existing
  GROMACS output remains the default. Periodic Amber output is refused because
  Amber does not minimum-image a covalent bond crossing the chain-direction box
  boundary.

### Validation

- LAMMPS `run 0` passes for finite and periodic I-beta and periodic cellulose
  II slabs. GLYCAM bond, angle, torsion, Lennard--Jones, and electrostatic
  categories are compared independently with GROMACS. CHARMM36 bond,
  angle/Urey--Bradley, and torsion energies agree within 2e-5 relative and
  total potential energies within 5e-4 using matched official force-switch
  settings. Periodic GLYCAM totals agree with OpenMM within 5e-4.
- The 18-chain finite and 24-chain periodic I-beta presets pass GROMACS 2024.2
  preprocessing with both GLYCAM06j-1 and CHARMM36.
- A finite 18-chain DP6 dual-format GLYCAM build passes GROMACS preprocessing
  and completes an AmberTools23 `sander` smoke minimization with finite bonded
  energies.
- OpenMM 8.5 and GROMACS 2024.2 zero-step bonded and 1-4 energies agree within
  8e-6 relative, and total nonbonded energies within 1.3e-4 with matching
  cutoff, PME, and dispersion settings, for finite and periodic I-beta slabs,
  a periodic triclinic I-beta bulk crystal, and periodic I-alpha, II, and III-I
  slabs with GLYCAM06j-1 and CHARMM36. Short OpenMM NVT runs, and NPT runs with
  the chain axis fixed, keep every bonded distance within 0.026 nm of
  equilibrium.

### Documentation

- README reorganized: a consolidated "Validation" section gathers all
  validation evidence and states what is not validated; a new "Common errors"
  section lists every build error with its fix; a table records the cutoff and
  dispersion settings of the generated inputs for each force field and engine;
  the command reference covers every option.
- Tutorial: experimental OpenFF Sage systems, generated run inputs and
  restraint defines, the cellulose II hydrogen method, and new common errors.

## 1.1.2 - 2026-11-07

### Added

- Added a Conda `environment.yml` that installs CelluloseForge together with
  AmberTools and ParmEd for GLYCAM construction.

### Fixed

- Enforce a right-handed slab orientation matrix with determinant `+1`, while
  retaining chains along Cartesian `-y`, the requested surface normal along
  Cartesian `-z`, and low-to-high `x` chain numbering.
- Detect any nonzero `tleap` `Errors = N` report rather than only the literal
  values 1 or 2.
- Finite-chain GLYCAM and CHARMM36 slab or fibril builds now reject less than
  3 angstrom of longitudinal padding and warn below the documented 12
  angstrom baseline, preventing capped chain ends from overlapping their
  periodic images.
- Force-field builds now validate and report the closest same-chain contact
  across the longitudinal periodic boundary.
- Periodic CHARMM36 and GLYCAM mappings now derive residue order from the
  crystallographic O4--C1 connectivity for every supported allomorph.
- Periodic force-field builds reject any mapped O4--C1 glycosidic bond longer
  than 2 angstrom and record the maximum mapped linkage length in the build
  report.
- I-beta `(110)` slabs assign alternating origin and center conformations to
  their corresponding lattice sites. Explicit chain counts must be even, and
  physical-width requests round upward to a complete transverse lattice repeat.

### Changed

- Added `--padding-longitudinal` as the coordinate-independent chain-axis
  option. The historical `--padding-z` spelling remains as a deprecated alias.

### Validation

- Added residue-level D-cellulose chirality reporting and regression tests for
  unit-cell, bulk, slab, and fibril constructions.
- Added independent checks of mapped and written periodic glycosidic bonds,
  including a 50 kJ mol-1 per-atom initial bonded-energy ceiling in the
  dependency-gated GROMACS tests.
- Regenerated the I-beta `(110)` reference slab and reran finite and periodic
  GLYCAM06j-1 and CHARMM36 minimization and one-ps NVT smoke tests.
- Added an I-beta `(110)` site-assignment regression: origin and center chains
  alternate within each layer, no interchain contact is shorter than 1.7
  angstrom, and the two central layers retain a 0.780 angstrom minimum-image
  RMSD in the periodic GLYCAM one-ps NVT test.

- Added direct C1 stereocentre-volume tests for positive and negative reference
  faces of all four supported allomorphs.
- Added a finite-GLYCAM terminal-contact regression test.
- Added periodic-image contact tests and a CLI regression test that refuses an
  unpadded finite force-field slab.
- Added allomorph-wide periodic glycosidic-bond regression tests.
- Periodic builds now independently parse the emitted GROMACS molecule
  topology, expand its `[ molecules ]` inventory, and validate every
  instantiated inter-residue O4--C1 bond that it contains.
- Build reports preserve the mapped-geometry and written-topology bond checks
  as separate results.
- GROMACS single-point integration tests require an initial bond energy below
  50 kJ mol-1 per atom.
- Added an I-beta `(110)` interchain-contact regression and a one-ps periodic
  GLYCAM/NVT core-geometry stability test.
- Build reports record a residue-level C1 signed-volume chirality validation.

## 1.1.1 - 2026-11-06

### Added

- Canonical GPL-3.0-or-later license metadata and the CelluloseForge Zenodo
  concept DOI.
- A project overview image in the README.

### Fixed

- GROMACS scientific-validation CI now uses Ubuntu 24.04 and a bounded,
  checksum-verified download of the pinned CHARMM36 distribution.

This patch release does not change the scientific models or generated
topologies from version 1.1.0.

## 1.1.0 - 2026-10-23

### Added

- Coordinate-only builds write a sibling `<structure>.build_report.json`
  recording the interpreted geometry, chain count, derived degree of
  polymerization, and coordinate filename.
- Slab sizing is selected independently for each axis, allowing hybrid
  requests such as an exact layer count with minimum width and chain length.
  Requested and realized dimensions are recorded in the build report.

### Changed

- Slabs use a simulation-oriented Cartesian frame compatible with GROMACS
  semi-isotropic pressure coupling: chains lie along `-y`, and the outward
  surface normal and vacuum direction lie along `-z`. Slabs built with 1.0.0
  used a `-y` surface normal and chains along `z`; coordinates for identical
  requests therefore differ by a rigid rotation.
- Fibrils follow the same convention, with chains along Cartesian `y` and the
  cross-section in the `xz` plane. Slab and fibril layers are numbered from low
  to high Cartesian `z`. `--padding-y` now pads the second transverse direction
  (Cartesian `z`), and `--padding-z` pads along the chain direction; the option
  names are retained for compatibility.
- Position restraints apply only to heavy atoms of the selected chains;
  hydrogens remain unrestrained. The build report records this in its
  `atom_selection` field.

### Fixed

- Type annotation of the coordinate-only build-report path.

## 1.0.0 - 2026-10-15

Initial stable release.

### Added

- Crystallographic definitions for cellulose I-alpha, I-beta, II, and III-I.
- Unit-cell and periodic bulk-crystal construction.
- Chain-parallel slab construction for cellulose I-alpha, I-beta, II, and
  III-I, with validated reference workflows for I-alpha `(010)`, I-beta
  `(110)`, II `(110)`, and III-I `(110)`.
- Literature-defined 36-chain fibrils for cellulose I-alpha, I-beta, and II.
- Finite and covalently periodic chain-boundary treatments.
- PDB and mmCIF coordinate output, including automatic mmCIF fallback for
  systems exceeding the PDB atom-serial limit.
- GLYCAM06j-1 mapping and matched GROMACS GRO/TOP output through AmberTools and
  ParmEd.
- Native finite- and periodic-chain CHARMM36 GROMACS topology generation from
  a checksum-pinned official force-field release, without requiring VMD or
  `psfgen` at runtime.
- Deterministic chain and transverse-layer position-restraint selection for
  GLYCAM and CHARMM36 systems.
- Machine-readable build reports, deterministic identifiers, automated tests,
  and Linux/macOS continuous integration.

### Validation

- 137 automated tests cover crystallography, geometry, coordinate formats,
  GLYCAM and CHARMM36 mapping, finite and periodic interactions, restraints,
  and GROMACS output generation; eight external-dependency cases are skipped
  when their required tools are unavailable.
- Finite and periodic GLYCAM06j-1 and CHARMM36 reference slabs for all four
  supported allomorphs pass GROMACS preprocessing, minimization, and short NVT
  validation workflows.
- Periodic bulk systems for all four allomorphs and finite or periodic 36-chain
  fibrils for the implemented presets have reproducible GROMACS validation
  workflows.
- A finite 102,000-atom cellulose I-beta `(110)` slab passes GROMACS
  preprocessing with PME and zero warnings.
- The built wheel installs and generates representative structures in clean
  Python 3.11, 3.12, and 3.13 environments.

### Scope limitations

- OPLS, chemical functionalization, native CHARMM/PSF output, and native Amber
  topology output are not included.
- Slabs are restricted to faces parallel to the crystallographic chain axis.
  Faces outside the validated reference matrix should be assessed for their
  intended application.
- Fibril presets are limited to cellulose I-alpha, I-beta, and II.
