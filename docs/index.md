# CelluloseForge

CelluloseForge constructs experimentally grounded crystalline cellulose and
chitin unit cells, bulk crystals, chain-parallel slabs, and literature-defined
fibrils, plus coordinate-only Mannan I crystals. Supported systems can be
parameterized with GLYCAM06j-1, finite- or periodic-chain CHARMM36, or
experimental OpenFF Sage 2.3, and written for GROMACS, OpenMM, LAMMPS, and
(finite GLYCAM chains) Amber, each with starting run inputs. The CHARMM36
workflow does not require VMD or `psfgen` at runtime.

Start with the [step-by-step tutorial](tutorial.md). The main repository README
contains the installation summary, command reference, a consolidated
"Validation" section, and a "Common errors" section. Engine-specific details
are in [`openmm.md`](openmm.md), [`lammps.md`](lammps.md), [`amber.md`](amber.md),
[`charmm36.md`](charmm36.md), [`openff.md`](openff.md), and [`opls.md`](opls.md); the validation
scripts and results are documented in `validation/README.md`; planned features
remain separate in `plan.md`.
