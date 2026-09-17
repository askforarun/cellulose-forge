# Cellulose I-alpha unit cell

`structure.pdb` is a geometry-only diagnostic structure generated from the
published triclinic P1 fractional coordinates retained by Cellulose-Builder.
The unit cell contains one cellulose chain and two anhydroglucose residues
(42 atoms). The cellulose chain follows crystallographic `a`.

Regenerate it with:

```bash
uv run cellulose-forge build \
  --allomorph I-alpha \
  --shape unit-cell \
  --output examples/i_alpha_unit_cell/structure.pdb
```

The coordinates have not been energy-minimized. Force-field topology is not
embedded in this diagnostic PDB.
