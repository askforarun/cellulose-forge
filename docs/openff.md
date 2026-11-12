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

Supported geometries are finite and covalently periodic slabs and covalently
periodic bulk crystals of cellulose I-alpha, I-beta, II and III-I. GROMACS and
OpenMM outputs are supported. Fibrils, chitin, mannan, Amber, LAMMPS,
restraints, and finite bulk crystals are intentionally rejected.

```bash
cellulose-forge build \
  --allomorph I-beta --shape slab --face 1 1 0 \
  --layers 2 --chains-per-layer 4 --degree-polymerization 6 \
  --chain-boundary periodic --vacuum 30 \
  --force-field openff --md-engine gromacs --md-engine openmm \
  --output openff_ibeta_slab
```

OpenFF Sage was fitted with a homogeneous-fluid dispersion correction. It is
therefore enabled for bulk crystals and disabled for vacuum-containing slabs,
consistently in the serialized OpenMM system and generated GROMACS inputs.

## Validation boundary

Model construction is validated: graph connectivity, stereochemistry,
periodic bonded terms, NAGL charge equivalence, and term-by-term GROMACS versus
OpenMM energies. Generated systems also undergo minimization and short NVT
smoke tests. These checks do **not** establish the physical accuracy or
transferability of OpenFF Sage for crystalline cellulose; that scientific
assessment is outside the software validation scope.
