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

GROMACS, OpenMM, NAMD and LAMMPS use the built-in adapters selected in `manifest.json`; they
do not need a copied runner script. Other engines retain a subfolder runner
during the staged migration. Each adapter performs an energy evaluation at the
supplied coordinates without minimisation or dynamics and writes `energy.json`. Energies use kJ/mol and the
canonical names listed in `manifest.example.json`. The manifest's
`required_terms` prevents a false pass based on only a small overlapping
subset. An adapter is also where a
documented engine convention is normalized; it must record every correction
in the JSON rather than silently changing a value.

The OpenMM entry must point to an unconstrained single-point `System` and give
an explicit force classification. Force and interaction indices are zero-based.
Every harmonic-bond and torsion force must be represented, even when its list
of special interactions is empty. Production `system.xml` files with hydrogen
constraints are intentionally refused.

The NAMD adapter verifies that its binary coordinates reproduce the shared GRO
coordinates and that the input explicitly sets `mergeCrossterms no`,
`LJcorrection off`, `rigidBonds none` and `run 0`. Nonzero `BOUNDARY` or `MISC`
terms require an explicit classification in the manifest.

The LAMMPS adapter checks that `system.data` contains the shared GRO
coordinates before executing `run 0`. Periodic impropers and GLYCAM explicit
pair energies must be supplied from the topology when present; they cannot be
recovered uniquely from a LAMMPS thermo row.

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
The manifest may define per-term `[absolute_kj_mol, relative]` tolerances;
this is useful for net electrostatics, whose real-space and reciprocal-space
components can be much larger than their sum.

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
