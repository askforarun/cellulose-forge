# Contributing

Contributions should preserve the separation between crystallographic geometry, molecular connectivity, force-field assignment, and GROMACS serialization.

## Development workflow

1. Open an issue describing the scientific or technical change.
2. Record significant design choices in `docs/decisions/`.
3. Add or update tests alongside implementation changes.
4. Run the test, lint, and type-check commands documented in `README.md`.
5. Document scientific data provenance and numerical tolerances.

Generated reference data must include the command, tool version, source force-field version, and checksum needed to reproduce it.

