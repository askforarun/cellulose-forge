# Supported-scope audit

This file maps the Version 1 user-facing claims to executable or archived
evidence. Features not listed here must be described as exploratory or future
work.

| Capability | Supported scope | Evidence |
|---|---|---|
| Unit-cell geometry | I-alpha, I-beta, II, III-I | Unit and CLI tests for cell matrices, symmetry, labels, counts, and coordinate output |
| Periodic bulk geometry | I-alpha, I-beta, II, III-I | Replication tests and GLYCAM/CHARMM36 allomorph validation matrices |
| Slab geometry | I-beta chain-parallel `(h,k,0)` | `(110)` is the validated reference; other accepted faces remain exploratory |
| Fibril geometry | Literature-defined 36-chain I-alpha, I-beta, and II presets | Geometry tests and force-field fibril matrices |
| GLYCAM06j-1 | Finite and periodic pristine cellulose with GROMACS output | Mapping/topology tests, allomorph matrix, large-system preprocessing, and dependency-gated `grompp` integration test |
| CHARMM36 | Finite and periodic pristine cellulose with GROMACS output | Exact periodic and finite topology comparisons, allomorph/fibril matrices, PME/NVT smoke tests, and dependency-gated `grompp` integration tests |
| Position restraints | Selected deterministic chain IDs for GLYCAM and CHARMM36 builds | Unit tests and `position_restraint_validation.json` |
| Large coordinates | Automatic mmCIF fallback above the PDB atom-serial limit; GRO for parameterized systems | Unit tests and 102,000-atom GLYCAM validation |
| Clean installation | Python 3.11 and 3.13 on Linux and macOS | CI wheel-install matrix; local isolated Python 3.13 wheel test |

Explicit exclusions for Version 1 are functionalized cellulose, OPLS,
native CHARMM/PSF output, native Amber topology output, arbitrary branched
glycans, and slab construction for allomorphs other than I-beta.

The SoftwareX manuscript and README must remain no broader than this table.
