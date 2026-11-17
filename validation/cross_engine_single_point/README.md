# Cross-engine single-point energy bundle

A release bundle has one shared coordinate source and one subfolder per MD
engine:

```text
cross_engine_single_point/
  manifest.json
  run_single_point_validation.py
  shared/coordinate_reference.json
  gromacs/
  openmm/
  lammps/
  namd/
  amber/                 # finite GLYCAM systems only
```

The shared JSON contains only ordered positions and full box vectors in
angstrom. It is owned by the validation bundle, not by any MD engine. Each
adapter reads its own native coordinate file and proves that the exact file
named by the engine input matches this neutral reference. Consequently an
Amber evaluation never needs a GRO file or a GROMACS installation.

All five engines use the built-in adapters selected in `manifest.json`; they
do not need copied runner scripts. Each adapter performs an energy evaluation at the
supplied coordinates without minimisation or dynamics and writes `energy.json`. Energies use kJ/mol and the
canonical names listed in `manifest.example.json`. The manifest's
`required_terms` prevents a false pass based on only a small overlapping
subset. An adapter is also where a
documented engine convention is normalized; it must record every correction
in the JSON rather than silently changing a value.

Each engine entry must select exactly one execution mode. Release bundles use
an `adapter` plus typed `inputs`. The legacy `command` mode remains available
only for externally maintained adapters during migration; mixing both modes or
providing neither is refused.

The OpenMM entry must point to an unconstrained single-point `System` and give
an explicit force classification. Force and interaction indices are zero-based.
Every harmonic-bond and torsion force must be represented, even when its list
of special interactions is empty. Production `system.xml` files with hydrogen
constraints are intentionally refused.

The NAMD adapter resolves the configuration's `binCoordinates`, verifies that
binary file and all three `cellBasisVector` values against the neutral reference,
and requires the input to set `mergeCrossterms no`,
`LJcorrection off`, `rigidBonds none` and `run 0`. Nonzero `BOUNDARY` or `MISC`
terms require an explicit classification in the manifest.

The LAMMPS adapter resolves the input's `read_data` command and checks that
native data file against the neutral positions and box before executing `run 0`.
Periodic impropers and GLYCAM explicit
pair energies must be supplied from the topology when present; they cannot be
recovered uniquely from a LAMMPS thermo row.

The Amber adapter is limited to finite-chain GLYCAM systems. It verifies the
native `inpcrd` positions and six-value box against the neutral reference and requires an `imin=1`,
`maxcyc=0` input. The periodic-improper contribution must be supplied
explicitly in kcal/mol, including an explicit zero for cellulose.

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
The example uses the live-tested global gate of 0.1 kJ/mol or 1e-4 relative,
with a 2e-4 relative electrostatics override for the documented PME/PPPM
summation difference. These are validation tolerances, not claims about
force-field accuracy.

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

To evaluate one installed engine without comparing it to another, use
`--evaluate-only`. For example, this requires only CelluloseForge and
AmberTools:

```bash
python run_single_point_validation.py --evaluate-only --engine amber
```

Any engine can be evaluated independently with
`--evaluate-only --engine ENGINE`; no other MD engine is required. If
`--engine` is omitted, the runner attempts every engine listed in the
manifest. A missing executable is then a failure for a required entry and a
recorded skip only for an entry explicitly marked optional.

For NAMD and LAMMPS, evaluation is deliberately limited to the declarative
single-point command subset emitted by CelluloseForge. The adapters validate
both command names and their energy-relevant arguments before execution and
refuse includes, scripting, coordinate mutation, dynamics, or unverified
settings. This is a validation contract for generated inputs, not a general
parser for arbitrary user-written NAMD or LAMMPS scripts.

Not every force field has every engine. Amber is finite-chain GLYCAM only;
native NAMD currently applies to CHARMM36. A downloadable validation release
therefore contains several system bundles rather than claiming that one
topology format is valid for every engine and force field.
