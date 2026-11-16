# Cross-engine single-point energy bundle

A release bundle has one shared coordinate source and one subfolder per MD
engine:

```text
cross_engine_single_point/
  manifest.json
  run_single_point_validation.py
  shared/
  gromacs/
  openmm/
  lammps/
  namd/
  amber/                 # finite GLYCAM systems only
```

Each engine subfolder contains its native inputs and a `run_single_point.py` adapter.
The adapter performs an energy evaluation at the supplied coordinates without
minimisation or dynamics and writes `energy.json`. Energies use kJ/mol and the
canonical names listed in `manifest.example.json`. The manifest's
`required_terms` prevents a false pass based on only a small overlapping
subset. An adapter is also where a
documented engine convention is normalized; it must record every correction
in the JSON rather than silently changing a value.

The canonical breakdown is deliberately aggregate and has no synonyms:
`bond`; `angle` (including Urey--Bradley); `dihedral` (proper torsions);
`improper`; `vdw` (ordinary plus 1--4 Lennard-Jones); `electrostatics`
(real-space, reciprocal-space, and 1--4); `restraint`;
`dispersion_correction`; and `potential`.

Every engine report must be a complete decomposition of `potential`, not only
the subset named in `required_terms`. The potential must equal the sum of all
other reported canonical terms within the configured numerical tolerance.
Restraint and dispersion-correction contributions must therefore be reported
separately whenever they are active and included in `potential`; an unknown
nonzero contribution is an error rather than being silently absorbed.

Run the bundle with:

```bash
python run_single_point_validation.py
```

Use `--no-run` to compare existing `energy.json` files. The output
`single_point_energy_comparison.json` contains the raw normalized energies, the common
terms compared against the reference engine, absolute and relative
differences, and pass/fail status. Missing executables cause an explicit
failure unless an engine is marked optional in the manifest.

Select only installed engines by repeating `--engine`; include the reference
engine named in the manifest. For example, an Amber/GROMACS bundle can be run
without installing OpenMM, LAMMPS, or NAMD:

```bash
python run_single_point_validation.py --engine gromacs --engine amber
```

Not every force field has every engine. Amber is finite-chain GLYCAM only;
native NAMD currently applies to CHARMM36. A downloadable validation release
therefore contains several system bundles rather than claiming that one
topology format is valid for every engine and force field.
