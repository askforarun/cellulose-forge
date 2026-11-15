# NAMD output (experimental)

CelluloseForge can write native NAMD inputs for existing CHARMM36 cellulose
and anhydrous alpha/beta-chitin systems. This backend is experimental until
the stored NAMD/GROMACS/OpenMM energy matrix and dynamics checks are complete.

The first milestone supports the finite and covalently periodic chain
geometries already supported by CHARMM36. It writes:

- `system.psf`, an extended X-PLOR PSF with every chain expanded explicitly;
- `structure_namd.pdb`, in exactly the PSF atom order;
- `par_all36_carb.prm`, copied from a canonical CHARMM parameter file supplied
  by the user; and
- `run_namd.conf`, a conservative minimization and 1 ps NVT smoke test.

Example:

```bash
cellulose-forge build \
  --polymer cellulose --allomorph I-beta \
  --shape slab --face 1 1 0 --layers 4 --chains-per-layer 8 \
  --degree-polymerization 8 --chain-boundary periodic --vacuum 30 \
  --force-field charmm36 --md-engine namd \
  --charmm-parameter-file /path/to/par_all36_carb.prm \
  --output ibeta_namd
```

The parameter file is not redistributed. It must be from the same CHARMM36
release as the installed GROMACS port. CelluloseForge copies it into the build
folder and does not attempt to reconstruct a CHARMM parameter file from the
GROMACS port.

The generated input uses PME, `exclude scaled1-4`, `oneFourScaling 1.0`, the
native CHARMM 10--12 Angstrom force switch, and no isotropic Lennard-Jones tail
correction. `wrapAll off` preserves chains bonded through a periodic boundary.
The initial configuration is deliberately not GPU-resident because NAMD 3's
GPU-resident mode does not support this force-switching setting.

## Safety checks and current limits

Before writing, CelluloseForge requires the explicit CHARMM36 1-4 pair list to
match the unique bond-graph distance-three pairs exactly. Missing, extra, or
duplicate pairs stop the build. Full triclinic basis vectors are written in
Angstrom, so tilted bulk cells are retained.

NAMD restraints, NPT, GLYCAM, OpenFF, OPLS-AA, chitin dihydrate, and mannan are
outside this first milestone. The build report labels NAMD runtime validation
as pending until the NAMD executable has parsed and run the full acceptance
matrix.

