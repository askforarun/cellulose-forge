# CelluloseForge

CelluloseForge is a Python project for constructing crystalline cellulose structures and producing simulation-ready GROMACS systems parameterized with CHARMM36 or GLYCAM.

The project is in its foundation phase. Scientific functionality will be added incrementally according to [plan.md](plan.md), beginning with a validated cellulose I-beta crystallographic prototype.

## Version 1 scope

- Cellulose I-alpha, I-beta, II, and III-I.
- Periodic bulk crystals, Miller-index surface slabs, and published nanofibrils.
- Arbitrary crystallographic faces selected by signed Miller indices, plus published cross-section presets.
- Finite and covalently periodic cellulose chains.
- CHARMM36 and GLYCAM force-field adapters.
- GROMACS GRO, TOP, and ITP output.

Functionalization, Antechamber, OPLS, native PSF, and native Amber topology output are deferred until after Version 1.

## Development setup

Using `uv`:

```bash
uv sync --extra dev
uv run pytest
uv run ruff check .
uv run mypy src
```

Run the initial CLI:

```bash
uv run cellulose-forge --version
```

Generate an I-beta periodic bulk supercell by repeating the crystallographic
unit cell along its `a`, `b`, and `c` axes:

```bash
uv run cellulose-forge build \
  --allomorph I-beta \
  --shape bulk \
  --cells 2 2 2 \
  --output ibeta_2x2x2.pdb
```

This diagnostic PDB contains complete crystallographic unit cells. It does not
yet contain bonds or force-field parameters.

Generate a chain-parallel cellulose I-beta `(110)` visualization slab with four
layers, 20 chains per layer, DP 12, and 30 angstrom of vacuum:

```bash
uv run cellulose-forge build \
  --allomorph I-beta \
  --shape slab \
  --face 1 1 0 \
  --layers 4 \
  --chains-per-layer 20 \
  --degree-polymerization 12 \
  --vacuum 30 \
  --output ibeta_110_slab.pdb
```

Degree of polymerization is the number of anhydroglucose units (AGUs) per
chain. The slab implementation currently requires even DP and retains complete
crystallographic residues. Cellulose I-beta `(110)` is validated against the SI
wetting structure. Other chain-parallel I-beta `(h k 0)` faces can be generated
for exploration, but their surface termination and layer arrangement still
require independent validation. Exploratory status is recorded in the PDB
header. The SI wetting reference corresponds to four layers, 100 chains per
layer, and DP 12; that 100,800-atom cellulose substrate exceeds the standard
PDB atom-serial limit and will be written as GROMACS GRO once that writer is
implemented.

## Legacy references

The repository retains the original Cellulose-Builder and doGlycans distributions as validation references. They are not imported by the new Python package. See [legacy/README.md](legacy/README.md) for provenance and licensing notes.

## License

This project is distributed under the GNU General Public License, version 3 or later. See [LICENSE](LICENSE).
