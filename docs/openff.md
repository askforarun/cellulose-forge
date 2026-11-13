# Experimental OpenFF Sage support

CelluloseForge can construct experimental OpenFF Sage 2.3 systems without
AmberTools. Install the optional environment from the repository root:

```bash
conda env create -f environment-openff.yml
conda activate cellulose-forge-openff
```

The backend constructs the beta-(1->4)-D-glucan molecular graph directly from
CelluloseForge atom names and connectivity. Charges are assigned at build time
with NAGL model `openff-gnn-am1bcc-1.0.0.pt`. The build fails unless the chain
is neutral, equivalent atoms in periodic residues have identical charges, and
every C1--C5 stereocentre has the independently recorded D-cellulose parity.

## Scope

Supported geometries are finite and covalently periodic slabs, finite and
covalently periodic fibrils, and covalently periodic bulk crystals. Slabs and
bulk crystals cover cellulose I-alpha, I-beta, II and III-I. Fibrils use the
same builders as GLYCAM and CHARMM36: the literature 36-chain fibril for
I-alpha, I-beta and II, and the I-beta 18- and 24-chain presets and custom
`--fibril-profile` row profiles. GROMACS and OpenMM outputs are supported.
Chitin, mannan, Amber, LAMMPS, restraints, and finite bulk crystals are
intentionally rejected.

```bash
cellulose-forge build \
  --allomorph I-beta --shape slab --face 1 1 0 \
  --layers 2 --chains-per-layer 4 --degree-polymerization 6 \
  --chain-boundary periodic --vacuum 30 \
  --force-field openff --md-engine gromacs --md-engine openmm \
  --output openff_ibeta_slab

cellulose-forge build \
  --allomorph I-beta --shape fibril --fibril-model 36-chain \
  --degree-polymerization 6 --chain-boundary periodic --vacuum 20 \
  --force-field openff --md-engine gromacs --md-engine openmm \
  --output openff_ibeta_fibril
```

The outputs match the other force fields: `minim.mdp` and `nvt.mdp` for
GROMACS (with Sage's 0.8--0.9 nm potential switch), and for OpenMM the same
configurable `run_openmm.py`, `system.xml`, `state.xml` and `topology.cif` as
GLYCAM and CHARMM36 builds.

OpenFF Sage was fitted with a homogeneous-fluid dispersion correction. It is
therefore enabled for bulk crystals and disabled for slabs and fibrils, whose
transverse vacuum violates its assumption, consistently in the serialized
OpenMM system and generated GROMACS inputs. This is the same rule used for
GLYCAM06; CHARMM36 never uses the correction. As for every force field, the
periodic box must exceed twice the nonbonded cutoff plus a 0.5 A pair-list
buffer on each side (19 A for Sage's 0.9 nm cutoff) in every perpendicular
direction, and a slab or fibril with less transverse vacuum than
the cutoff triggers a warning.

## Validation boundary

Model construction is validated: graph connectivity, stereochemistry,
periodic bonded terms, NAGL charge equivalence, and term-by-term GROMACS versus
OpenMM energies. Generated systems also undergo minimization and short NVT
smoke tests. The generated `run_openmm.py` is a smoke-test runner, not a
production simulation protocol. These checks do **not** establish the physical accuracy or
transferability of OpenFF Sage for crystalline cellulose; that scientific
assessment is outside the software validation scope.
