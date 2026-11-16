# 0007: Cross-engine single-point energy drivers

- Status: accepted for incremental implementation
- Branch: `refactor/cross-engine-energy-drivers`
- Scope: evaluation and comparison only; topology generation is unchanged

## Context

CelluloseForge already validates GROMACS, OpenMM, LAMMPS, NAMD, and Amber
energies, but adapters and comparison rules are duplicated across validation
scripts. OpenFF Interchange provides a useful driver pattern: one engine
adapter returns normalized terms, and one report layer compares them.
Interchange is not adopted as the universal parameterized-system hub: its
non-OpenFF imports are experimental, it has no NAMD writer, and replacing the
term-matched CHARMM36, GLYCAM, and OPLS paths would be a scientific regression.

## Decision

1. Add an internal `cellulose_forge.energy_drivers` package. It evaluates files
   produced by existing writers and never changes a Hamiltonian.
2. Aggregate terms are `bond`, `angle`, `dihedral`, `improper`, `vdw`,
   `electrostatics`, and `potential`, in kJ/mol. `angle` includes
   Urey--Bradley; `dihedral` means proper torsions; `vdw` and `electrostatics`
   include their 1--4 terms. `potential` is independently reported.
3. A detailed CHARMM view may additionally expose `urey_bradley`, `lj_14`,
   `coulomb_14`, and separate ordinary nonbonded terms. The aggregate schema
   never folds impropers into proper torsions.
4. Corrections remain visible with raw, correction, and comparable values,
   reason, and reference. Raw engine output is never overwritten.
5. Comparisons use `max(absolute, relative * scale)`, retaining existing
   CelluloseForge tolerances rather than Interchange's molecular-scale default.
6. Required engines fail validation if unavailable or unsuccessful. Only
   engines declared optional before execution may be skipped, with a reason.
7. Reports converge on a versioned identity containing coordinate, atom-order,
   topology/system, box-vector, and physical-settings fingerprints, plus
   engine and force-field versions. Comparison refuses incompatible identities.
   Schema-1 reports remain readable only during migration.
8. OpenFF remains a parameterization route, not an engine driver. Its GROMACS
   and OpenMM representations use those engine adapters.
9. The package remains internal during migration and is not exported from the
   package root until a separate API-stability decision.

## Migration

1. Extract terms, parsing, and comparison from
   `validation/cross_engine_single_point` without changing its golden JSON.
2. Wrap GROMACS and OpenMM, then LAMMPS, NAMD, and finite-chain Amber.
3. Preserve LAMMPS convention corrections and detailed CHARMM diagnostics.
4. Migrate one validator at a time and reproduce its previous result before
   retiring duplicate code.
5. Adopt the strict identity schema after every adapter supplies all hashes.

## Non-goals

- importing CHARMM36, GLYCAM, or OPLS into Interchange;
- replacing existing topology writers;
- cloning Interchange's complete parameterized object model;
- silently changing settings, term aggregation, or tolerances;
- coarse-grained mapping, backmapping, or MCP behavior.

## Validation requirements

- Existing comparison JSON remains numerically and structurally unchanged
  during the first extraction.
- Missing or unknown terms, wrong units, nonnumeric/nonfinite values,
  coordinate mismatch, and missing potential energy fail explicitly.
- Required-engine failures cannot become skips.
- Convention corrections are tested before and after correction.
- Engine adapters use identical coordinates, atom order, cell vectors, and
  matched cutoff, switching, electrostatics, dispersion, and constraints.
