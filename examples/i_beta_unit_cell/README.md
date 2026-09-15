# Cellulose I-beta unit cell

Generate the diagnostic unit-cell PDB with:

```bash
uv run cellulose-forge build \
  --allomorph I-beta \
  --shape unit-cell \
  --output examples/i_beta_unit_cell/i_beta_unit_cell.pdb
```

The structure contains 84 atoms: four 21-atom anhydroglucose residues arranged as one cellobiose repeat in each of two crystallographic chains. PDB chain `A` denotes the origin chain and chain `B` denotes the center chain.

Coordinates are intentionally unwrapped so that each crystallographic residue remains geometrically intact across cell boundaries. This is a geometry-only diagnostic PDB and does not yet include covalent connectivity or force-field parameters.

The fractional coordinates and unit-cell parameters are from Nishiyama, Langan, and Chanzy, *J. Am. Chem. Soc.* **2002**, 124, 9074-9082, DOI: [10.1021/ja0257319](https://doi.org/10.1021/ja0257319), transcribed through the retained Cellulose-Builder reference data.
