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

This diagnostic coordinate file contains complete crystallographic unit cells.
It does not yet contain bonds or force-field parameters. Use a `.cif` or
`.mmcif` output extension to request PDBx/mmCIF explicitly.

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
  --chain-boundary finite \
  --vacuum 30 \
  --output ibeta_110_slab.pdb
```

Alternatively, request minimum physical dimensions in angstrom:

```bash
uv run cellulose-forge build \
  --allomorph I-beta \
  --shape slab \
  --face 1 1 0 \
  --slab-width 120 \
  --slab-thickness 25 \
  --slab-length 60 \
  --chain-boundary periodic \
  --vacuum 30 \
  --output ibeta_110_sized_slab.pdb
```

The exact-count and physical-size options are mutually exclusive. Physical
dimensions are rounded upward to complete chains, molecular layers, and
cellobiose repeats. Requested and realized dimensions are recorded in PDB
remarks or mmCIF comments.

Degree of polymerization is the number of anhydroglucose units (AGUs) per
chain. The slab implementation currently requires even DP and retains complete
crystallographic residues. Cellulose I-beta `(110)` is validated against the SI
wetting structure. Other chain-parallel I-beta `(h k 0)` faces can be generated
for exploration, but their surface termination and layer arrangement still
require independent validation. Exploratory status is recorded in the
coordinate-file metadata. Following the original Cellulose-Builder convention, `(100)` and
`(010)` expose origin-chain layers, whereas `(200)` and `(020)` expose
center-chain layers located halfway through the unit cell. The selected chain
type is recorded in the coordinate-file metadata. Their lateral periodic box lengths are
exact lattice multiples: `(100)/(200)` repeat along `b`, while `(010)/(020)`
repeat along `a`. The SI wetting reference corresponds to four layers, 100
chains per layer, and DP 12; that 100,800-atom cellulose substrate exceeds the
standard PDB atom-serial limit. When a requested `.pdb` would exceed 99,999
atoms, CelluloseForge writes a sibling `.cif` file instead and reports its path.
It never writes a silently wrapped or truncated PDB. PDBx/mmCIF also preserves
full atom, residue, and chain identifiers for large-system visualization.

Slab output uses a visualization-friendly Cartesian convention: the transverse
surface direction is `x`, the outward surface normal is `-y`, and the cellulose
chain direction is `z`. Consequently, the slab layers are visible directly in
the default `xy` viewing plane without manually rotating the coordinates.

Slab chains use deterministic one-based numbering. Layer 1 is the requested
low-`y` face; layers then increase inward toward high `y`. Within each layer,
chains are numbered from low `x` to high `x`. PDB one-character chain labels are
reused after 62 chains, but the four-character segment field records the stable
chain number as `C001`, `C002`, and so on. Restraint and topology generation must
use the numeric chain mapping rather than the recycled one-character label.

Generate the published cellulose I-beta 36-chain elementary fibril with six
cellobiose repeats (DP 12) per chain:

```bash
uv run cellulose-forge build \
  --allomorph I-beta \
  --shape fibril \
  --fibril-model 36-chain \
  --degree-polymerization 12 \
  --chain-boundary periodic \
  --vacuum 30 \
  --output ibeta_36chain_dp12.pdb
```

The `36-chain` model reproduces the explicit, allomorph-specific chain
selections published with Cellulose-Builder; it is not generated as an ideal
regular hexagon. Geometry presets are available for cellulose I-alpha,
I-beta, and II. Cellulose III-I is rejected with a clear error because the
reference work did not define a 36-chain III-I fibril. `--vacuum` is the total
padding added independently across each transverse box dimension. Chains are
numbered deterministically from 1 through 36 in the published selection order,
and the exact source fragment list is embedded in the output metadata.

The periodic DP12 I-beta fibril has also been mapped to the GLYCAM reference
chain topology as 36 repeated molecules. GROMACS 2024.2 accepts it with PME
electrostatics and zero preprocessing warnings, and steepest-descent
minimization reaches the configured 1000 kJ mol-1 nm-1 force threshold.

Generate the same fibril as a complete GLYCAM/GROMACS system:

```bash
uv run cellulose-forge build \
  --allomorph I-beta \
  --shape fibril \
  --fibril-model 36-chain \
  --degree-polymerization 12 \
  --chain-boundary periodic \
  --vacuum 30 \
  --force-field glycam \
  --output ibeta_fibril
```

With `--force-field glycam`, `--output` names a directory rather than a
coordinate file. CelluloseForge uses AmberTools `tleap` with
`leaprc.GLYCAM_06j-1` to generate and validate one chain, uses ParmEd to
convert its topology, maps the crystallographic coordinates into GLYCAM atom
order, and writes:

```text
ibeta_fibril/
|-- structure.pdb
|-- structure.gro
|-- system.top
`-- build_report.json
```

AmberTools and ParmEd are external requirements only for GLYCAM builds.
Coordinate-only builds (`--force-field none`, the default) do not require
them. If ParmEd is installed in another environment, pass its interpreter with
`--parmed-python /path/to/python`. `build_report.json` records the force-field
release, DP, boundary condition, chain count, output files, validation results,
and executable provenance.

## Legacy references

Cellulose-Builder and doGlycans are cited as external scientific references;
their source distributions, publications, and binary examples are not
vendored in this repository. See [legacy/README.md](legacy/README.md) for
provenance and citations.

## License

This project is distributed under the GNU General Public License, version 3 or later. See [LICENSE](LICENSE).
