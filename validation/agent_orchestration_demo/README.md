# External-agent build demonstration

This directory records one end-to-end demonstration that a natural-language
build request can be translated into a correct `cellulose-forge build`
invocation and independently verified, rather than trusted on the agent's
narration alone.

CelluloseForge does not embed an LLM. This is a proof of concept showing how an
external agent can use the same documented CLI and machine-readable reports as
a human or conventional workflow manager; it is not a benchmark of general
natural-language accuracy.

## Natural-language request

> Build me the validated cellulose I-beta (110) slab reference structure with
> 4 layers, 20 chains per layer, DP 12, finite chains, and 30 angstrom of
> vacuum, ready for GROMACS with the GLYCAM force field.

## Agent-translated command

```bash
cellulose-forge doctor --tleap-executable /path/to/tleap --parmed-python /path/to/python
```

reports dependency availability first (see
[`cellulose-forge doctor`](../../docs/tutorial.md)),
then:

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape slab \
  --face 1 1 0 \
  --layers 4 \
  --chains-per-layer 20 \
  --degree-polymerization 12 \
  --chain-boundary finite \
  --vacuum 30 \
  --force-field glycam \
  --tleap-executable /path/to/tleap \
  --parmed-python /path/to/python \
  --output ibeta_110_slab_agent_build
```

This is the same validated `(110)` slab example documented in the tutorial and
README (Section 3 of the manuscript). Using an established reference request
makes the agent-produced command and its resulting metadata directly
checkable; this single case does not establish performance for arbitrary
natural-language requests.

## Independent verification

The command was executed against a real AmberTools/GLYCAM environment
(AmberTools with `leaprc.GLYCAM_06j-1`, ParmEd 4.3.0). The resulting
`build_report.json` (kept in this directory; the large `structure.pdb`,
`structure.gro`, and `system.top` outputs are regenerable and excluded per the
project's convention for bulk generated artifacts) shows. Local executable
paths in the retained report were replaced with `<AMBERTOOLS_ENV>` and
`<PARMED_ENV>` placeholders before publication; all scientific values are
unchanged.

- `force_field_atom_count: 20400` for 80 chains of a 255-atom finite GLYCAM
  chain topology (`255 * 80 = 20400`), matching `chain_count: 80`;
- `reference_generation.validation.total_charge` equal to zero to floating
  point precision, confirming charge conservation in the generated topology;
- `SCEE_SCALE_FACTOR` and `SCNB_SCALE_FACTOR` both `1.0`, matching GLYCAM's
  documented unity 1-4 electrostatic and Lennard-Jones scaling;
- the force-field release and conversion software recorded for provenance.

These are the same independent checks (charge conservation, correct scaling
factors, atom-count consistency) used elsewhere in the validation suite; they
do not depend on the agent's own report of success.

## Reproducing this demonstration

```bash
cellulose-forge doctor --json
```

should report `"glycam_gromacs_files_ready": true` before attempting the
build above. If `tleap` or ParmEd are not on `PATH`, pass
`--tleap-executable` and `--parmed-python` to both `doctor` and `build`.

`doctor` checks dependency availability; it does not replace construction-time
validation or `gmx grompp` for a generated system.
