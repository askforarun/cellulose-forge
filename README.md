# CelluloseForge

CelluloseForge is a Python project for constructing crystalline cellulose structures and producing simulation-ready GROMACS systems parameterized with GLYCAM.

The project is in its foundation phase. Scientific functionality will be added incrementally according to [plan.md](plan.md), beginning with a validated cellulose I-beta crystallographic prototype.

## Version 1 scope

- Cellulose I-alpha, I-beta, II, and III-I.
- Parallelepiped crystallites, surface-oriented slabs, fibrils, nanocrystallites, and supported monolayers.
- Arbitrary crystallographic faces selected by signed Miller indices, plus published cross-section presets.
- Finite and covalently periodic cellulose chains.
- GLYCAM parameterization represented correctly for GROMACS.
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

## Legacy references

The repository retains the original Cellulose-Builder and doGlycans distributions as validation references. They are not imported by the new Python package. See [legacy/README.md](legacy/README.md) for provenance and licensing notes.

## License

This project is distributed under the GNU General Public License, version 3 or later. See [LICENSE](LICENSE).
