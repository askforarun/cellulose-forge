# Changelog

All notable changes to CelluloseForge are documented in this file.

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

### Changed

- Added `--padding-longitudinal` as the coordinate-independent chain-axis
  option. The historical `--padding-z` spelling remains as a deprecated alias.

### Validation

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
