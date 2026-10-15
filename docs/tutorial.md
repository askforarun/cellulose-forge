# CelluloseForge tutorial

This tutorial builds progressively from a single cellulose unit cell to
simulation-ready GLYCAM and CHARMM36 GROMACS systems. It is written for users who are new to
CelluloseForge but already have basic familiarity with molecular simulation
files.

## Table of contents

1. [Choose an installation mode](#1-choose-an-installation-mode)
2. [Inspect the command-line interface](#2-inspect-the-command-line-interface)
3. [Generate and inspect a unit cell](#3-generate-and-inspect-a-unit-cell)
4. [Replicate a periodic bulk crystal](#4-replicate-a-periodic-bulk-crystal)
5. [Build a finite I-beta slab](#5-build-a-finite-i-beta-110-slab)
6. [Request a slab by physical dimensions](#6-request-a-slab-by-physical-dimensions)
7. [Build other allomorph slabs](#7-build-other-allomorph-slabs)
8. [Build a 36-chain elementary fibril](#8-build-a-36-chain-elementary-fibril)
9. [Understand chain numbering](#9-understand-chain-numbering)
10. [Generate a GLYCAM/GROMACS fibril](#10-generate-a-glycamgromacs-fibril)
11. [Generate CHARMM36/GROMACS systems](#11-generate-charmm36gromacs-systems)
12. [Preprocess and minimize with GROMACS](#12-preprocess-and-minimize-with-gromacs)
13. [Work with systems above 99,999 atoms](#13-work-with-systems-above-99999-atoms)
14. [Visualize the results](#14-visualize-the-results)
15. [Common errors](#15-common-errors)
16. [Use an external workflow manager or LLM agent](#16-use-an-external-workflow-manager-or-llm-agent)
17. [What to record in a research workflow](#17-what-to-record-in-a-research-workflow)

## 1. Choose an installation mode

Coordinate-only builds require Python and NumPy. GLYCAM builds additionally
require AmberTools and ParmEd. CHARMM36 builds require the checksum-pinned
official GROMACS CHARMM36 distribution, installed with
`cellulose-forge forcefield install charmm36`, but do not require a reference
ITP, VMD, or `psfgen`. GROMACS is needed only when preprocessing or simulating
the generated files.

### Coordinate-only environment

```bash
git clone https://github.com/askforarun/cellulose-forge.git
cd cellulose-forge
uv sync
uv run cellulose-forge --version
```

When using this environment, prefix tutorial commands with `uv run`. For
example, `cellulose-forge build ...` becomes `uv run cellulose-forge build ...`.

### GLYCAM environment

```bash
conda create -n cellulose-forge -c conda-forge \
  python=3.12 ambertools parmed pip
conda activate cellulose-forge
git clone https://github.com/askforarun/cellulose-forge.git
cd cellulose-forge
python -m pip install -e .
```

### CHARMM36 environment and force-field installation

CHARMM36 construction needs the normal Python installation of CelluloseForge,
but does not need AmberTools, ParmEd, VMD, or `psfgen`. After installing the
package, install the pinned official GROMACS force-field release once:

```bash
cellulose-forge forcefield install charmm36
cellulose-forge forcefield path charmm36
```

The first command downloads the February 2026 CHARMM36/CGenFF 5.0 archive from
the MacKerell Laboratory, verifies the archive checksum, and validates the
installed files. The second command prints the directory that subsequent
builds will use.

On a cluster without internet access, copy the official archive to the cluster
and install it locally:

```bash
cellulose-forge forcefield install charmm36 \
  --archive /path/to/charmm36-feb2026_cgenff-5.0.ff.tgz
```

The default installation root is `~/.local/share/cellulose-forge`. To use a
project or shared-filesystem location, set the override consistently when
installing and building:

```bash
export CELLULOSE_FORGE_DATA_DIR=/path/to/cellulose-forge-data
cellulose-forge forcefield install charmm36
```

GROMACS is not needed to construct the files, but it is required for `grompp`,
minimization, and simulation validation.

Check that the required programs are visible with a single diagnostic
command:

```bash
cellulose-forge doctor
```

This reports Python, NumPy, `tleap`, ParmEd, and GROMACS availability, states
which build modes are ready, and prints a suggested fix for anything missing.
Pass `--json` for machine-readable output, which is convenient for scripts or
agents that need to check readiness programmatically before choosing a build
mode:

```bash
cellulose-forge doctor --json
```

If `tleap` or ParmEd live in a different environment than the one running
`cellulose-forge`, point `doctor` at them the same way `build` accepts them:

```bash
cellulose-forge doctor \
  --tleap-executable /path/to/tleap \
  --parmed-python /path/to/python
```

If GROMACS is unavailable, you can still generate GLYCAM files and transfer
them to a machine or cluster where GROMACS is installed; `doctor` reports
GROMACS preprocessing readiness separately from GLYCAM file-generation
readiness for exactly this reason.

## 2. Inspect the command-line interface

```bash
cellulose-forge build --help
```

Every build requires three choices:

1. an allomorph, such as `I-beta`;
2. a shape, such as `bulk`, `slab`, or `fibril`;
3. an output path.

Coordinate-only mode is the default. It is equivalent to explicitly supplying
`--force-field none`.

Create a directory for the tutorial results:

```bash
mkdir cellulose_forge_tutorial
```

## 3. Generate and inspect a unit cell

Start with cellulose I-beta:

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape unit-cell \
  --output cellulose_forge_tutorial/ibeta_unit_cell.pdb
```

The PDB contains a `CRYST1` record describing the crystallographic box. Inspect
the header and count the atom records:

```bash
head -n 8 cellulose_forge_tutorial/ibeta_unit_cell.pdb
grep -c '^HETATM' cellulose_forge_tutorial/ibeta_unit_cell.pdb
```

The same command writes
`cellulose_forge_tutorial/ibeta_unit_cell.build_report.json`. Coordinate-only
reports record the chain count, degree of polymerization, geometry annotations,
and the coordinate file actually written. For dimension-driven slabs, inspect
this report to find the derived DP and the requested versus realized dimensions.

The unit-cell file is a crystallographic coordinate model. In coordinate-only
mode it does not contain a molecular-mechanics topology.

The other available allomorph names are `I-alpha`, `II`, and `III-I`:

```bash
cellulose-forge build \
  --allomorph II \
  --shape unit-cell \
  --output cellulose_forge_tutorial/cellulose_ii_unit_cell.pdb
```

## 4. Replicate a periodic bulk crystal

Build a `4 × 4 × 5` I-beta supercell:

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape bulk \
  --cells 4 4 5 \
  --output cellulose_forge_tutorial/ibeta_4x4x5_bulk.pdb
```

`--cells 4 4 5` means four crystallographic unit cells along `a`, four along
`b`, and five along `c`. It does not select a surface face. A bulk build retains
the crystallographic basis and repeats complete unit cells.

Check the generated metadata:

```bash
head -n 8 cellulose_forge_tutorial/ibeta_4x4x5_bulk.pdb
```

Use bulk geometry when the intended system is periodic in all three crystal
directions. Surface simulations should instead use a slab.

## 5. Build a finite I-beta `(110)` slab

The exact-count interface directly specifies molecular layers, chains in each
layer, and degree of polymerization:

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
  --output cellulose_forge_tutorial/ibeta_110_finite.pdb
```

The important arguments are:

- `--face 1 1 0`: expose the `(110)` crystallographic face;
- `--layers 4`: retain four complete molecular layers;
- `--chains-per-layer 20`: place 20 chains in each layer;
- `--degree-polymerization 12`: retain 12 anhydroglucose units per chain;
- `--chain-boundary finite`: do not join the last residue to the first;
- `--vacuum 30`: add a total of 30 Å normal to the slab.

The coordinate-only structure contains crystallographic atoms. Force-field
mode later constructs the additional terminal atoms required for finite GLYCAM
chains.

Cellulose II is a separate documented case: its diffraction dataset omits
HO2, HO3, and HO6. GLYCAM output generates those hydrogen coordinates using
the aligned force-field hydroxyl directions while preserving GLYCAM O--H bond
lengths. It does not alter the crystallographic heavy atoms. Cellulose II's
origin and center chains are mapped in opposite residue directions so their
antiparallel beta(1->4) connectivity remains correct. The generated-coordinate
method is included in `build_report.json`.

CelluloseForge orients generated slabs consistently: the transverse surface
direction is `x`, the outward surface normal is `-y`, and the cellulose chain
direction is `z`.

## 6. Request a slab by physical dimensions

If physical size is more useful than molecular counts, request minimum width,
thickness, and length in angstrom:

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
  --output cellulose_forge_tutorial/ibeta_110_sized.pdb
```

The requested dimensions are minima. CelluloseForge rounds upward to retain
whole chains, whole molecular layers, and whole cellobiose repeats. Requested
and realized dimensions are written to the coordinate metadata.

Do not combine `--layers` or `--chains-per-layer` with physical-dimension
options in the same command.

## 7. Build other allomorph slabs

For example, generate `(100)` and `(200)` slabs:

```bash
cellulose-forge build \
  --allomorph I-beta --shape slab --face 1 0 0 \
  --layers 4 --chains-per-layer 20 --degree-polymerization 12 \
  --chain-boundary finite --vacuum 30 \
  --output cellulose_forge_tutorial/ibeta_100.pdb

cellulose-forge build \
  --allomorph I-beta --shape slab --face 2 0 0 \
  --layers 4 --chains-per-layer 20 --degree-polymerization 12 \
  --chain-boundary finite --vacuum 30 \
  --output cellulose_forge_tutorial/ibeta_200.pdb
```

`(100)` and `(200)` have parallel normals but represent different plane
spacings and chain-layer selections. CelluloseForge therefore preserves the
exact Miller indices instead of reducing `(200)` to `(100)`.

The allomorph validation matrix uses I-alpha `(010)`, I-beta `(110)`, II
`(110)`, and III-I `(110)`. For example:

```bash
cellulose-forge build \
  --allomorph I-alpha --shape slab --face 0 1 0 \
  --layers 3 --chains-per-layer 6 --degree-polymerization 6 \
  --chain-boundary finite --vacuum 24 --padding-z 12 \
  --force-field charmm36 --output ialpha_010_slab

cellulose-forge build \
  --allomorph II --shape slab --face 1 1 0 \
  --layers 3 --chains-per-layer 6 --degree-polymerization 6 \
  --chain-boundary periodic --vacuum 24 \
  --force-field glycam --output cellulose_ii_110_slab
```

For I-alpha the crystallographic chain axis is `a`, so a chain-parallel face
requires `h=0`. For I-beta, II, and III-I the chain axis is `c`, so it requires
`l=0`. Finite slabs can use `--padding-z` to separate explicit chain ends from
their periodic images. Periodic slabs reject nonzero chain-end padding.

## 8. Build a 36-chain elementary fibril

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape fibril \
  --fibril-model 36-chain \
  --degree-polymerization 12 \
  --chain-boundary periodic \
  --vacuum 30 \
  --output cellulose_forge_tutorial/ibeta_36chain_dp12.pdb
```

The model uses the published discrete 36-chain selection. It is not generated
by cutting an ideal regular hexagon. Presets are available for I-alpha, I-beta,
and II. A III-I 36-chain model is intentionally rejected because the reference
builder did not define that preset.

`--vacuum 30` adds 30 Å of total padding along both transverse box directions,
with 15 Å on each opposing side. Use `--padding-x` and `--padding-y` when the
two transverse dimensions need different padding. These options override
`--vacuum` only for their respective direction.

A fibril may also contain finite, chemically terminated chains. In that case,
`--padding-z` adds total space beyond the two chain ends:

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape fibril \
  --fibril-model 36-chain \
  --degree-polymerization 12 \
  --chain-boundary finite \
  --padding-x 30 \
  --padding-y 40 \
  --padding-z 20 \
  --output cellulose_forge_tutorial/ibeta_36chain_finite_padded.pdb
```

For a periodic fibril, CelluloseForge derives the z box length from the degree
of polymerization and crystallographic chain repeat. It rejects nonzero
`--padding-z`, which would otherwise separate atoms joined by the periodic
covalent bond.

## 9. Understand chain numbering

Slab chains are numbered deterministically:

1. layer 1 begins at the requested low-`y` face;
2. layer numbers increase toward high `y`;
3. within a layer, chains increase from low `x` to high `x`.

Fibril layers follow the stored literature preset from low to high Cartesian
`y`; within each layer, chain IDs increase from low to high Cartesian `x`.
These stable numeric identifiers are the appropriate basis for later restraint
or selection files.

### Selective chain position restraints

Force-field builds can create separate free and restrained cellulose molecule
types. Select either deterministic chain IDs or complete transverse layers;
the selectors are mutually exclusive:

```bash
uv run cellulose-forge build \
  --allomorph I-beta \
  --shape fibril \
  --fibril-model 36-chain \
  --degree-polymerization 40 \
  --chain-boundary periodic \
  --force-field glycam \
  --restrain-layers 1 3 \
  --restraint-force 1000 1000 1000 \
  --output restrained_fibril
```

Slab and fibril layers are numbered from low to high Cartesian `y`. Fibrils
are displayed with their cross-section rotated counterclockwise by 90 degrees;
their layers follow the explicit bottom-to-top rows of the literature-defined
preset, and chain IDs within each layer increase from low to high Cartesian
`x`. The generated topology places
free molecules before restrained molecules, and `build_report.json` records
the layer-to-chain expansion and mapping from original chain ID to GROMACS
molecule index and global atom range. The same selectors work with
`--force-field charmm36`.
Enable the restraints in an MDP file with
`define = -DPOSRES_CELLULOSE` and supply the generated coordinates to
`gmx grompp` with `-r structure.gro`. Without this define, the same topology
runs with the restraint section disabled.
One-character PDB chain labels repeat after 62 chains, so use segment identifiers
or mmCIF identifiers for larger systems.

## 10. Generate a GLYCAM/GROMACS fibril

Activate the environment containing AmberTools and ParmEd, then run:

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape fibril \
  --fibril-model 36-chain \
  --degree-polymerization 12 \
  --chain-boundary periodic \
  --vacuum 30 \
  --force-field glycam \
  --output cellulose_forge_tutorial/ibeta_fibril_glycam
```

Unlike coordinate-only mode, the output is a directory:

```text
ibeta_fibril_glycam/
├── structure.pdb
├── structure.gro
├── system.top
└── build_report.json
```

Inspect the report:

```bash
python -m json.tool \
  cellulose_forge_tutorial/ibeta_fibril_glycam/build_report.json
```

The report records the force-field release, degree of polymerization, boundary
condition, chain and atom counts, external executable provenance, charge check,
one-four scaling, and generated filenames.

If `tleap` or ParmEd is installed outside the active environment:

```bash
cellulose-forge build [the remaining build options] \
  --force-field glycam \
  --tleap-executable /path/to/tleap \
  --parmed-python /path/to/python \
  --output cellulose_forge_tutorial/custom_environment_build
```

## 11. Generate CHARMM36/GROMACS systems

CHARMM36 supports pristine finite and periodic chains. Install the pinned
February 2026 GROMACS force-field release once:

```bash
cellulose-forge forcefield install charmm36
```

Confirm which installed directory CelluloseForge will use:

```bash
cellulose-forge forcefield path charmm36
```

### Finite I-beta (110) slab

This example constructs explicit reducing and nonreducing chain termini:

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
  --force-field charmm36 \
  --output cellulose_forge_tutorial/ibeta_110_finite_charmm36
```

### Periodic bulk crystal

For a bulk crystal, the last residue of every chain is covalently linked to
the first residue through the periodic box. The corresponding boundary bonds,
angles, pairs, and dihedrals are generated explicitly:

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape bulk \
  --cells 4 4 6 \
  --chain-boundary periodic \
  --force-field charmm36 \
  --output cellulose_forge_tutorial/ibeta_bulk_charmm36
```

### Periodic 36-chain fibril

Build the fibril without a reference-chain topology:

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape fibril \
  --fibril-model 36-chain \
  --degree-polymerization 8 \
  --chain-boundary periodic \
  --vacuum 30 \
  --force-field charmm36 \
  --output cellulose_forge_tutorial/ibeta_fibril_charmm36
```

CelluloseForge reads `BGLC`, applies the cyclic `14bb` linkage, generates the
chain topology at DP8, maps the crystallographic
coordinates, and writes `structure.gro`, `cellulose_charmm36.itp`, `system.top`,
the local force-field directory, and `build_report.json`. This is direct
GROMACS topology generation: VMD and `psfgen` are not invoked.

Use `--chain-boundary finite` for explicit reducing and nonreducing termini, or
`--chain-boundary periodic` for a final-to-first glycosidic linkage. Periodic bulk models are tested
for I-alpha, I-beta, II, and III-I; periodic 36-chain fibrils are tested for
I-alpha, I-beta, and II.

Each CHARMM36 output directory contains:

```text
structure.pdb                 crystallographic/force-field coordinates
structure.gro                 GROMACS coordinates
cellulose_charmm36.itp        generated cellulose topology
system.top                    complete GROMACS system topology
charmm36-feb2026_cgenff-5.0.ff/ pinned force-field files
build_report.json             interpreted options, provenance, and identifiers
```

Inspect the report before simulation:

```bash
python -m json.tool \
  cellulose_forge_tutorial/ibeta_110_finite_charmm36/build_report.json
```

Preprocess and minimize the finite slab:

```bash
gmx grompp \
  -f validation/finite_pme_minim.mdp \
  -c cellulose_forge_tutorial/ibeta_110_finite_charmm36/structure.gro \
  -p cellulose_forge_tutorial/ibeta_110_finite_charmm36/system.top \
  -o cellulose_forge_tutorial/ibeta_110_finite_charmm36/minim.tpr

gmx mdrun \
  -deffnm cellulose_forge_tutorial/ibeta_110_finite_charmm36/minim
```

Use `validation/periodic_pme_minim.mdp` for periodic bulk crystals and
periodic fibrils. Do not add chain-direction vacuum to a periodic chain because
that would separate atoms joined by the final-to-first glycosidic bond.

## 12. Preprocess and minimize with GROMACS

From the repository root, preprocess the periodic fibril:

```bash
gmx grompp \
  -f validation/periodic_pme_minim.mdp \
  -c cellulose_forge_tutorial/ibeta_fibril_glycam/structure.gro \
  -p cellulose_forge_tutorial/ibeta_fibril_glycam/system.top \
  -o cellulose_forge_tutorial/ibeta_fibril_glycam/minim.tpr
```

Do not use `-maxwarn` to hide unexplained topology warnings. Investigate any
warning before simulation.

Run minimization:

```bash
gmx mdrun \
  -deffnm cellulose_forge_tutorial/ibeta_fibril_glycam/minim
```

Inspect the final potential energy and maximum force:

```bash
grep -E 'Potential Energy|Maximum force' \
  cellulose_forge_tutorial/ibeta_fibril_glycam/minim.log
```

The included MDP files are validation and smoke-test protocols. They are not a
substitute for selecting and reporting an application-specific equilibration
and production protocol.

The same commands apply to a CHARMM36 output directory by replacing the GRO,
TOP, and output paths.

## 13. Work with systems above 99,999 atoms

The published-size example below contains 100,800 crystallographic atoms:

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape slab \
  --face 1 1 0 \
  --layers 4 \
  --chains-per-layer 100 \
  --degree-polymerization 12 \
  --chain-boundary finite \
  --vacuum 30 \
  --output cellulose_forge_tutorial/large_slab.pdb
```

Although `.pdb` was requested, CelluloseForge writes `large_slab.cif` because
PDB atom serials cannot uniquely represent the system. This is reported on the
terminal.

For a GLYCAM build, add `--force-field glycam` and make `--output` a directory.
The resulting finite system contains 102,000 atoms after terminal atoms are
added. Use:

- `structure.cif` for visualization and unique identifiers;
- `structure.gro` and `system.top` for GROMACS;
- the compiled `.tpr` for simulation.

GRO display identifiers wrap after 99,999, but GROMACS reads coordinate records
in topology order. CelluloseForge has validated this 102,000-atom system with
GROMACS 2024.2 `gmx grompp`, PME electrostatics, zero warnings, and no
`-maxwarn` allowance.

## 14. Visualize the results

PDB and GRO files can be opened in common molecular viewers such as VMD. For
large mmCIF structures, use a viewer version with PDBx/mmCIF support. When a
viewer shows bonds inferred from distance, remember that coordinate-only files
do not yet carry the force-field topology; the authoritative simulation bonds
are in `system.top` after a GLYCAM or CHARMM36 build.

Useful views are:

- looking along `z` to inspect the slab or fibril cross-section;
- looking along `y` to inspect the slab surface plane;
- displaying the periodic box to confirm vacuum placement.

## 15. Common errors

### `--output must be a directory when --force-field glycam is selected`

Use a directory name without `.pdb`, `.cif`, or `.gro`:

```bash
--force-field glycam --output my_glycam_system
```

### `tleap` is not found

Run `cellulose-forge doctor` to confirm which executable is missing and where
CelluloseForge looked for it, then activate the AmberTools environment or
provide its full path:

```bash
--tleap-executable /path/to/amber/bin/tleap
```

### ParmEd cannot be imported

Run `cellulose-forge doctor --parmed-python /path/to/python` to check a
specific interpreter, then install ParmEd in the active environment or
identify a Python interpreter that contains it:

```bash
--parmed-python /path/to/amber/environment/bin/python
```

### CHARMM36 is not installed

Install the pinned release from the official source:

```bash
cellulose-forge forcefield install charmm36
```

Use `cellulose-forge forcefield path charmm36` to inspect its location. A local
official archive can be supplied with `--archive` for offline installation.

### A slab face is rejected

Slab faces must be parallel to the cellulose chain axis: `h=0` for I-alpha and
`l=0` for I-beta, II, and III-I. The zero index `(0,0,0)` is never a valid
crystallographic plane.

### A fibril model is rejected for III-I

The current literature reference does not provide a supported III-I 36-chain
preset. Use unit-cell or bulk geometry for III-I instead.

### The requested PDB was not created

Check the terminal message for a sibling `.cif` path. CelluloseForge switches
to mmCIF rather than writing invalid wrapped PDB identifiers.

## 16. Use an external workflow manager or LLM agent

CelluloseForge does not contain an LLM and does not require one. External
workflow managers and LLM-based software agents can orchestrate it through the
same deterministic CLI used by a human. Machine-readable diagnostics allow an
agent to inspect dependency availability before choosing a build mode:

```bash
cellulose-forge doctor --json
```

Consider this natural-language request:

> Build a validated cellulose I-beta (110) slab with 4 layers, 20 chains per
> layer, DP 12, finite chains, 30 angstrom vacuum, and GLYCAM output.

Its corresponding CLI action is:

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
  --output cellulose_forge_tutorial/ibeta_110_slab
```

The command is ordinary CelluloseForge syntax: it can be inspected, rerun, and
validated without the agent. After execution, check `build_report.json` rather
than relying on the agent's description of success. For GLYCAM output, also run
the GROMACS preprocessing procedure in Section 12. The retained
proof-of-concept evidence is documented in
[`validation/agent_orchestration_demo/`](../validation/agent_orchestration_demo/).

This example demonstrates interface compatibility, not general
natural-language accuracy. Applications that accept unrestricted user text
should validate allowed values and require human review before launching costly
simulation jobs.

## 17. What to record in a research workflow

For reproducibility, retain:

- the complete CelluloseForge command;
- the CelluloseForge version or Git commit;
- `build_report.json` for parameterized builds or the sibling
  `<structure>.build_report.json` for coordinate-only builds;
- the selected force field and its exact release;
- for GLYCAM, the AmberTools and ParmEd versions;
- for CHARMM36, the pinned release name and force-field archive checksum;
- the GROMACS version;
- the MDP files and exact `gmx grompp` command;
- checksums or an archive DOI for large generated reference files.

See the
[validation record](https://github.com/askforarun/cellulose-forge/blob/main/validation/README.md)
for currently completed checks and the
[project plan](https://github.com/askforarun/cellulose-forge/blob/main/plan.md)
for features that remain future work.
