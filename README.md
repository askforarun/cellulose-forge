# CelluloseForge

CelluloseForge is a crystallography-first Python tool for building atomistic
cellulose structures. It generates periodic bulk crystals, cellulose I-beta
surface slabs, and literature-defined 36-chain fibrils. Supported systems can
also be parameterized with GLYCAM06j-1 and written directly for GROMACS.

> **Development status:** CelluloseForge is pre-alpha research software. The
> supported and validated boundaries below are intentionally narrower than the
> long-term project plan.

## Key features

- Crystallographic datasets for cellulose I-alpha, I-beta, II, and III-I.
- Periodic bulk replication along the crystallographic `a`, `b`, and `c` axes.
- Chain-parallel I-beta slabs selected with signed Miller indices `(h,k,0)`.
- Literature-defined 36-chain fibrils for I-alpha, I-beta, and II.
- Finite or covalently periodic cellulose chains.
- Deterministic chain and residue numbering.
- PDB and PDBx/mmCIF coordinate output.
- Automatic mmCIF fallback when a requested PDB exceeds 99,999 atoms.
- GLYCAM06j-1 parameterization through AmberTools and ParmEd.
- GROMACS `.gro` and `.top` output plus a machine-readable build report.

## Current support

| Structure | Available allomorphs | Validation status |
|---|---|---|
| Unit cell | I-alpha, I-beta, II, III-I | Crystallographic cell, symmetry, and atom counts tested |
| Periodic bulk | I-alpha, I-beta, II, III-I | Replication and box vectors tested |
| Slab | I-beta chain-parallel `(h,k,0)` faces | `(110)` validated; other faces are exploratory |
| 36-chain fibril | I-alpha, I-beta, II | Geometry presets tested; periodic I-beta DP12 validated with GLYCAM/GROMACS |

CHARMM36, functionalization, OPLS, native CHARMM/PSF output, and arbitrary
branched glycans are not implemented in the current release.

## Table of contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Verify the installation](#verify-the-installation)
- [Quick start](#quick-start)
- [Full tutorial](#full-tutorial)
- [GLYCAM and GROMACS workflow](#glycam-and-gromacs-workflow)
- [Command reference](#command-reference)
- [Output files](#output-files)
- [Large systems](#large-systems)
- [Scientific conventions](#scientific-conventions)
- [Development and testing](#development-and-testing)
- [Citation and provenance](#citation-and-provenance)

## Requirements

### Coordinate-only builds

- Python 3.11, 3.12, or 3.13
- NumPy 1.26 or newer, below 2.3

### GLYCAM/GROMACS builds

In addition to the coordinate-only requirements:

- AmberTools with `tleap`
- GLYCAM06j-1, distributed with AmberTools as `leaprc.GLYCAM_06j-1`
- ParmEd
- GROMACS for preprocessing, minimization, and simulation

AmberTools and ParmEd are not required when `--force-field none` is used.
GROMACS is not required to generate files, but it is required to validate and
simulate the generated system.

## Installation

CelluloseForge is not yet released on PyPI. Install it from the GitHub
repository.

### Coordinate-only installation with uv

```bash
git clone https://github.com/askforarun/cellulose-forge.git
cd cellulose-forge
uv sync
uv run cellulose-forge --help
```

Install [uv](https://docs.astral.sh/uv/getting-started/installation/) first if
the `uv` command is unavailable.

### Coordinate-only installation with pip

```bash
git clone https://github.com/askforarun/cellulose-forge.git
cd cellulose-forge
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
cellulose-forge --help
```

On Windows PowerShell, activate the environment with
`.venv\Scripts\Activate.ps1`.

### Recommended GLYCAM installation

Create one conda environment containing CelluloseForge, AmberTools, and ParmEd:

```bash
conda create -n cellulose-forge -c conda-forge \
  python=3.12 ambertools parmed pip
conda activate cellulose-forge
git clone https://github.com/askforarun/cellulose-forge.git
cd cellulose-forge
python -m pip install -e .
```

Verify the external tools:

```bash
which tleap
python -c "import parmed; print(parmed.__version__)"
gmx --version
```

GROMACS is commonly installed separately using the package manager or module
system appropriate to the workstation or cluster. Consult the
[GROMACS installation guide](https://manual.gromacs.org/current/install-guide/index.html).

If ParmEd is installed in a different Python environment, pass its interpreter
with `--parmed-python /path/to/python`. Likewise, a non-default `tleap` can be
selected with `--tleap-executable /path/to/tleap`.

## Verify the installation

```bash
cellulose-forge --version
cellulose-forge build --help
```

For a source checkout with development dependencies:

```bash
uv sync --extra dev
uv run pytest
uv run ruff check .
uv run mypy src
```

## Quick start

### Unit cell

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape unit-cell \
  --output ibeta_unit_cell.pdb
```

### Periodic bulk crystal

The values supplied to `--cells` are repetitions along the crystallographic
`a`, `b`, and `c` axes. Bulk geometry contains complete crystallographic unit
cells.

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape bulk \
  --cells 4 4 5 \
  --output ibeta_4x4x5_bulk.pdb
```

### I-beta `(110)` slab using exact counts

Degree of polymerization is the number of anhydroglucose units in each chain.
The current slab builder requires an even value.

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
  --output ibeta_110_slab.pdb
```

### I-beta `(110)` slab using physical dimensions

Dimensions are specified in angstrom. CelluloseForge rounds upward to complete
chains, molecular layers, and cellobiose repeats, and records the realized
dimensions in the output metadata.

```bash
cellulose-forge build \
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

Exact-count options and physical-dimension options are mutually exclusive.

### Literature-defined 36-chain fibril

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape fibril \
  --fibril-model 36-chain \
  --degree-polymerization 12 \
  --chain-boundary periodic \
  --vacuum 30 \
  --output ibeta_36chain_dp12.pdb
```

The preset is an explicit published chain arrangement, not an ideal regular
hexagon. Cellulose III-I has no supported 36-chain reference preset.

## Full tutorial

The [step-by-step tutorial](docs/tutorial.md) covers installation, unit cells,
bulk crystals, exact-count and dimension-based slabs, alternative faces,
36-chain fibrils, GLYCAM generation, GROMACS preprocessing and minimization,
large systems, visualization, chain numbering, and troubleshooting.

## GLYCAM and GROMACS workflow

Select `--force-field glycam` and provide an output directory instead of a
coordinate filename:

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape fibril \
  --fibril-model 36-chain \
  --degree-polymerization 12 \
  --chain-boundary periodic \
  --vacuum 30 \
  --force-field glycam \
  --output ibeta_fibril
```

CelluloseForge asks `tleap` to construct one GLYCAM chain, converts that
reference with ParmEd, maps the crystallographic coordinates into GLYCAM atom
order, and repeats the topology for the generated chains.

Preprocess a generated system with GROMACS:

```bash
gmx grompp \
  -f validation/periodic_pme_minim.mdp \
  -c ibeta_fibril/structure.gro \
  -p ibeta_fibril/system.top \
  -o ibeta_fibril/minim.tpr
```

Run energy minimization:

```bash
gmx mdrun -deffnm ibeta_fibril/minim
```

Use `validation/finite_pme_minim.mdp` for finite-chain systems and ensure that
the periodic box provides sufficient separation between chain ends.

## Command reference

```text
cellulose-forge build --allomorph ALLOMORPH --shape SHAPE [OPTIONS] --output PATH
```

| Option | Values or units | Meaning |
|---|---|---|
| `--allomorph` | `I-alpha`, `I-beta`, `II`, `III-I` | Cellulose crystal allomorph |
| `--shape` | `unit-cell`, `bulk`, `slab`, `fibril` | Structure type |
| `--cells NA NB NC` | positive integers | Bulk repetitions along `a`, `b`, and `c` |
| `--face H K L` | signed integers | Slab Miller index; current slabs require I-beta and `L=0` |
| `--layers` | positive integer | Complete molecular layers in an exact-count slab |
| `--chains-per-layer` | positive integer | Chains in each exact-count slab layer |
| `--degree-polymerization` | positive even integer | Anhydroglucose units per chain |
| `--slab-width` | angstrom | Minimum transverse width |
| `--slab-thickness` | angstrom | Minimum thickness normal to the surface |
| `--slab-length` | angstrom | Minimum chain-direction length |
| `--fibril-model` | `36-chain` | Named literature-defined fibril preset |
| `--chain-boundary` | `finite`, `periodic` | Terminal groups or a covalent bond across the chain boundary |
| `--vacuum` | angstrom | Total transverse vacuum padding |
| `--force-field` | `none`, `glycam` | Coordinate-only or GLYCAM/GROMACS output |
| `--tleap-executable` | executable path | `tleap` used for GLYCAM builds |
| `--parmed-python` | Python path | Interpreter containing ParmEd |
| `--output` | file or directory | Coordinate file for `none`; directory for `glycam` |

Run `cellulose-forge build --help` for the authoritative parser help.

## Output files

### Coordinate-only mode

The filename extension selects PDB or PDBx/mmCIF:

```text
structure.pdb
structure.cif
structure.mmcif
```

### GLYCAM mode

```text
ibeta_fibril/
├── structure.pdb       # or structure.cif for more than 99,999 crystal atoms
├── structure.gro       # GROMACS coordinates
├── system.top          # standalone GROMACS topology
└── build_report.json   # inputs, provenance, counts, and generated files
```

## Large systems

PDB atom serial numbers are limited to five columns. If a requested PDB would
contain more than 99,999 atoms, CelluloseForge writes a sibling `.cif` file and
prints its path. It never silently truncates or wraps a PDB.

For a parameterized system exceeding this limit:

- use `structure.cif` for visualization and unambiguous identifiers;
- use `structure.gro` with `system.top` for GROMACS;
- use the compiled `.tpr` for simulation.

GRO atom and residue identifiers are also five-column fields and therefore wrap
for very large systems. GROMACS reads the coordinate records in topology order,
but mmCIF should be used when unique human-readable identifiers matter.

A finite I-beta `(110)` reference slab with 100,800 crystallographic atoms was
mapped to 102,000 GLYCAM atoms and passed GROMACS 2024.2 `gmx grompp` with PME,
zero warnings, and no `-maxwarn` allowance.

## Scientific conventions

### Slab orientation

Slabs use a visualization-friendly Cartesian frame:

- transverse surface direction: `x`;
- outward surface normal: `-y`;
- cellulose chain direction: `z`.

The requested face is the low-`y` surface. Layers increase toward high `y`, and
chains within a layer are ordered from low `x` to high `x`.

### Chain identifiers

Chains have stable one-based numeric identifiers. PDB one-character chain labels
are reused after 62 chains, while the segment field records `C001`, `C002`, and
so on. Use numeric chain mappings or mmCIF identifiers for restraints and
analysis of large systems.

### I-beta faces

The `(110)` slab is validated against the reference wetting structure. Other
chain-parallel I-beta `(h,k,0)` surfaces are generated deterministically but are
currently marked exploratory. Following the Cellulose-Builder convention,
`(100)` and `(010)` expose origin-chain layers, whereas `(200)` and `(020)`
expose center-chain layers halfway through the unit cell.

## Development and testing

```bash
uv sync --extra dev
uv run pytest --cov=cellulose_forge --cov-report=term-missing
uv run ruff check .
uv run mypy src
```

Continuous integration tests Python 3.11 and 3.13 on Linux and macOS.

## Citation and provenance

Citation metadata is provided in [CITATION.cff](CITATION.cff). A release DOI
will be added after the first stable archival release. The JOSS manuscript is
maintained in [paper/paper.md](paper/paper.md).

The crystallographic datasets are derived from the primary diffraction
literature for cellulose I-alpha, I-beta, II, and III-I. Cellulose-Builder,
doGlycans, and Nanopolysaccharide Builder are related scientific tools and are
cited in the manuscript and [legacy provenance record](legacy/README.md). Their
source distributions are not required at runtime and are not vendored here.

## License

CelluloseForge is distributed under the GNU General Public License, version 3
or later. See [LICENSE](LICENSE).
