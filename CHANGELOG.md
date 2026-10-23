# Changelog

All notable changes to CelluloseForge are documented in this file.

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
