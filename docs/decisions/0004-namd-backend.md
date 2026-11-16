# NAMD backend: implementation and independent-review handoff

Status: in development on `feature/namd`  
Base: `release/1.3.0-rc` at `bb0c6d7`

## Scope

The first milestone adds native NAMD output for existing CHARMM36 cellulose
and anhydrous alpha/beta-chitin builds. It covers finite and covalently
periodic chains in the shapes already supported by the CHARMM36 backend.
GLYCAM, OpenFF, OPLS-AA, mannan, chitin dihydrate, NPT, and NAMD restraints are
outside this milestone.

Output consists of:

- an EXT X-PLOR `system.psf`;
- `structure_namd.pdb` in exactly the PSF atom order;
- a user-supplied canonical CHARMM carbohydrate parameter file copied as
  `par_all36_carb.prm`;
- `run_namd.conf`, initially a CPU/non-GPU-resident minimization and 1 ps NVT
  smoke test.

The canonical parameter file is required explicitly with
`--charmm-parameter-file`; CelluloseForge does not translate the installed
GROMACS port back into a CHARMM PRM and does not redistribute CHARMM files.

## Current implementation

- `src/cellulose_forge/namd.py` writes PSF/PDB/configuration/parameter files.
- `src/cellulose_forge/charmm36.py` exposes instantiated interaction indices
  so NAMD output does not parse a rendered GROMACS ITP.
- The PSF expands every chain into global atom indices and separates proper
  and improper dihedrals.
- Before writing, the explicit CHARMM36 pair list must equal the unique
  graph-distance-three endpoints after excluding 1-2 and 1-3 pairs. Output is
  refused for missing, extra, or duplicate pairs.
- Full triclinic basis vectors are written directly in angstrom.
- Native CHARMM settings use PME, `exclude scaled1-4`,
  `oneFourScaling 1.0`, a 10--12 angstrom force switch via
  `vdwForceSwitching on`, and no LJ tail correction.
- `wrapAll off` is deliberate for chains covalently bonded through the cell.
- NAMD builds require every perpendicular box width to exceed 28 Angstrom.
- The canonical February 2026 parameter file is enforced by SHA-256 and the
  digest is recorded in the build report.
- PSF interaction rows follow psfgen's grouping, and zero `NUMLP` and
  `NCRTERM` sections are present for VMD compatibility.
- Reusable scheduler-neutral comparison scripts live together in `validation/namd/`.

## Checks completed by Codex

- Focused writer, capability-profile and CHARMM36 regression suite: 23 tests
  passed; the full core suite also exits cleanly.
- Ruff and mypy pass for the new backend.
- A real NAMD-only I-beta (110) periodic CHARMM36 slab build produced the
  PSF/PDB/configuration/PRM bundle without retaining temporary GROMACS files.
- OpenMM 8.5 independently parsed both a 1,008-atom cellulose bundle and a
  1,296-atom beta-chitin bundle. In each case the PSF atom count, independently
  parsed PDB atom count, and created OpenMM System particle count agree, and
  the canonical February 2026 `par_all36_carb.prm` resolves every parameter.
- The CLI build report records NAMD files even when GROMACS was not requested.
- A real-build check caught and fixed the configuration's coordinate filename;
  it now names the PSF-ordered `structure_namd.pdb`.

NAMD 3.0.3 was subsequently installed on the development Mac. The stored
acceptance matrix now records real NAMD parsing, term-resolved single-point
comparison against GROMACS, minimization, and 1 ps NVT for periodic I-beta,
triclinic I-alpha, and alpha-chitin.

## Changes after independent review

The five review findings against `d674258` are addressed:

1. NAMD builds now require every perpendicular box width to exceed 28
   Angstrom, twice its 14 Angstrom pair-list distance.
2. PSF rows use psfgen-compatible grouping (four bonds, three angles, or two
   proper/improper dihedrals per line), with zero `NUMLP` and `NCRTERM`
   sections. Headless VMD loaded a real 2,592-atom beta-chitin PSF/PDB with all
   192 impropers.
3. The February 2026 `par_all36_carb.prm` SHA-256 is enforced and recorded.
4. Missing/wrong parameter files and NAMD restraint requests are rejected
   before an output directory is created. Supplying the parameter option
   without requesting NAMD is also rejected.
5. A real periodic CHARMM36 integration test compares the expanded NAMD PSF
   with the generated GROMACS topology. The independent script in
   `validation/namd/compare_bundle.py` additionally reproduced exact atom and
   interaction sets for real cellulose and chitin bundles.

Each generated folder also contains a portable shell runner, dedicated
unconstrained NAMD and GROMACS matched unswitched 12 Angstrom single-point
inputs, and a dependency-free comparison script. Both engines read coordinates
derived from the same GRO file. Missing terms, inconsistent energy sums, and
tolerance failures stop the comparison. Production inputs retain the native
CHARMM 10--12 Angstrom force switch and pass minimization plus 1 ps NVT in all
three stored cases, with finite energies and intact bonds.

## Independent review requested from Claude

Please review rather than rewrite first. In particular:

1. Parse the generated EXT PSF independently and compare every atom row and
   global bond/angle/proper/improper count with the expanded GROMACS topology.
2. Confirm that the graph-derived 1-4 equality test is correct for finite and
   periodic cellulose and GlcNAc chains, including the last-O4/first-C1
   closure.
3. Check PSF formatting against NAMD/psfgen output and identify any mandatory
   zero-count sections or EXT field-width issues.
4. Verify all nine components of triclinic cell vectors and the cell origin
   for I-alpha bulk and a rotated slab.
5. Check that `par_all36_carb.prm` contains every PSF atom type and every
   required bond, angle/Urey--Bradley, proper, and improper parameter.
6. Run NAMD parsing, minimization, and 1 ps NVT for the acceptance matrix below.
7. Compare NAMD energy components against double-precision GROMACS/OpenMM on
   identical high-precision coordinates. Isolate 1-4 energy if necessary.
8. Confirm that NAMD 3 is run outside GPU-resident mode because native CHARMM
   `vdwForceSwitching` is not supported by that mode.

## Minimum acceptance matrix

- Cellulose CHARMM36: finite and periodic I-beta (110) slabs; triclinic
  I-alpha bulk; cellulose-II periodic bulk or slab; one periodic 36-chain
  fibril.
- Chitin CHARMM36: alpha periodic bulk or slab and beta finite slab, including
  both acetyl impropers per residue.
- Every case: PSF/PDB order, neutral total charge, topology counts, NAMD parse,
  minimization, 1 ps NVT, finite energies, no constraint failure, and bond
  checks after dynamics.
- Periodic cases: every boundary O4--C1 bond and its associated angles and
  dihedrals must be present and remain geometrically valid by minimum image.

## Deferred items

- NAMD harmonic restraints require a separate constraint PDB and an explicit
  conversion from GROMACS `1/2 K r^2` to NAMD `k r^2`.
- GLYCAM has heterogeneous 1-4 scaling and cannot use this backend's single
  CHARMM exclusion/scaling policy without a separate representation.
- NPT is deferred until periodic-chain cell scaling is independently tested.
