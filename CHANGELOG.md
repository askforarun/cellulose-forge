# Changelog

All notable changes to CelluloseForge are documented in this file.

## 1.0.0 - Unreleased

Initial stable release candidate.

### Added

- Crystallographic definitions for cellulose I-alpha, I-beta, II, and III-I.
- Unit-cell and periodic bulk-crystal construction.
- Cellulose I-beta chain-parallel slab construction with a validated `(110)`
  reference workflow and explicitly marked exploratory alternative faces.
- Literature-defined 36-chain fibrils for cellulose I-alpha, I-beta, and II.
- Finite and covalently periodic chain-boundary treatments.
- PDB and mmCIF coordinate output, including automatic mmCIF fallback for
  systems exceeding the PDB atom-serial limit.
- GLYCAM06j-1 mapping and matched GROMACS GRO/TOP output through AmberTools and
  ParmEd.
- Machine-readable build reports, deterministic identifiers, automated tests,
  and Linux/macOS continuous integration.

### Validation

- Ninety automated tests cover crystallography, geometry, coordinate formats,
  GLYCAM mapping, periodic interactions, and GROMACS output generation.
- A periodic DP12 cellulose I-beta 36-chain fibril passes GROMACS preprocessing
  and energy minimization.
- A finite 102,000-atom cellulose I-beta `(110)` slab passes GROMACS
  preprocessing with PME and zero warnings.
- The built wheel installs and generates representative structures in clean
  Python 3.11, 3.12, and 3.13 environments.

### Scope limitations

- CHARMM36, OPLS, chemical functionalization, native CHARMM/PSF output, and
  native Amber topology output are not included.
- Slabs for cellulose I-alpha, II, and III-I are not supported in this release.
- Cellulose I-beta slab faces other than `(110)` remain exploratory.
