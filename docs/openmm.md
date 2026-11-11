# OpenMM output

CelluloseForge supports OpenMM for finite-chain GLYCAM systems. The adapter
uses the same complete-system `structure.prmtop` and `structure.inpcrd` files
as native Amber output, so atom ordering, coordinates, charges, bonded terms,
and periodic box remain identical across the two paths.

```bash
cellulose-forge build \
  --allomorph I-beta --shape fibril --fibril-model 18-chain \
  --degree-polymerization 12 --chain-boundary finite \
  --vacuum 30 --padding-longitudinal 20 \
  --force-field glycam --md-engine openmm \
  --output ibeta_18chain_openmm
```

## Generated files

- `structure.pdb`: mapped cellulose coordinates;
- `structure.prmtop`: complete GLYCAM topology;
- `structure.inpcrd`: complete coordinates and periodic box;
- `run_openmm.py`: standalone minimization and dynamics runner;
- `build_report.json`: Amber and OpenMM validation metadata.

The build-time smoke test creates an OpenMM PME system with hydrogen-bond
constraints, checks the initial potential energy is finite, requests ten
minimization iterations, and checks the resulting energy is finite. The report
records the OpenMM version, selected platform, particle and force counts,
initial and minimized energies, and Python executable.

Run the generated workflow from its output directory:

```bash
python run_openmm.py --minimize-iterations 1000 --steps 10000
```

The runner uses a 300 K Langevin-middle integrator, a 2 fs time step, PME with
a 1 nm cutoff, and constrained hydrogen bonds. These defaults are a portable
starting point rather than a production protocol; users should select ensemble,
coupling, reporting, and simulation length for their scientific application.

## Installation and boundary restriction

Install OpenMM with the project environment or the optional dependency:

```bash
conda env create -f environment.yml
# or
pip install 'cellulose-forge[openmm]'
```

Pass `--openmm-python /path/to/python` if OpenMM is installed under a different
interpreter. Covalently periodic chains are not currently emitted through this
route because the underlying Amber topology represents a bonded term crossing
the chain-direction periodic boundary without the required minimum-image
treatment. Use the validated GROMACS engine path for those systems.
