# Cellulose II unit cell

`structure.pdb` is a geometry-only diagnostic structure generated from the
published monoclinic P2_1 fractional coordinates retained by
Cellulose-Builder. It contains antiparallel origin and center chains, with two
residues per chain after application of the P2_1 screw operation (72 atoms).

```bash
uv run cellulose-forge build \
  --allomorph II \
  --shape unit-cell \
  --output examples/ii_unit_cell/structure.pdb
```

The diffraction dataset does not provide `HO2`, `HO3`, or `HO6` coordinates.
They are intentionally absent rather than guessed. Consequently this
diagnostic structure cannot yet be mapped to a complete GLYCAM topology.
