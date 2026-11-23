# ADR 0008: CHARMM Drude polarizable cellulose

Status: experimental on `feature/polarizable-drude`. The native topology
generator is implemented and checked against psfgen and OpenMM. In review,
NAMD 3.0.3 minimized a finite native bundle and completed a short dynamics run
with the corrected writer. The same finite I-beta case is retained as a gated
regression test that passes with NAMD 3.0.3 and the audited official streams.
Matched NAMD 3.0.3/OpenMM 8.6.1 single points pass for finite I-beta (110),
I-beta (100), and cellulose II (110) slabs. The topology has not been compared
with a CHARMM-generated reference, and production-length validation remains
open.

## Scientific source

CelluloseForge uses the official MacKerell Laboratory CHARMM Drude 2023a
streams. The carbohydrate stream defines `BGLC`, the beta-(1->4) `14bb`
patch, atom polarizabilities and Thole factors, off-centre lone pairs, and
anisotropic response.

The audited files are obtained from the MacKerell Laboratory's CHARMM Drude
force-field page (`drude_toppar_2023.tgz`). They are user-supplied and are not
redistributed by CelluloseForge. Their required filenames and SHA-256 hashes
are pinned in `drude_topology.py`.

The backend must not approximate this model by changing only partial charges.
A supported engine must reproduce all of the following:

1. parent--Drude oscillator pairs and their redistributed charge and mass;
2. anisotropic Drude response;
3. excluded-pair Thole screening;
4. massless off-centre lone-pair sites;
5. the linkage-specific atom, lone-pair, and anisotropy changes in `14bb`.

The force field was parameterized on monosaccharides and glycosidic linkages
in solution. It has not been validated for crystalline cellulose, so every
build report labels the output experimental.

## Native topology construction

NAMD documents that psfgen cannot in general construct a complete Drude PSF
and recommends CHARMM for this task. Requiring a licensed CHARMM executable
would conflict with CelluloseForge's portable open-source workflow. The
backend therefore parses the limited, audited subset of the official streams
needed for `BGLC` and `14bb` and writes the topology natively
(`cellulose_forge.drude_topology`, `cellulose_forge.drude`). No CHARMM, VMD or
psfgen executable is used.

### Audit pin

The audit applies to two files only, identified by SHA-256:

| File | SHA-256 |
| --- | --- |
| `toppar_drude_main_protein_2023a.str` | `083fa30a7bb00fdf0737a0bcfceb687c4939da0fb53047ab23054ca1016836ab` |
| `toppar_drude_carbohydrate_2023a.str` | `bfd45af59bce51742aa3e6b5075089e526e63befc86942471168804763d69cfb` |

Any other file is refused. The streams are not redistributed in the
repository; the unit tests use a small invented residue and linkage, and the
tests that need the official files run only when
`CELLULOSE_FORGE_DRUDE_TOPPAR` names their directory.

### Parser

The parser is fail-closed. Inside the two selected definitions it implements
`ATOM` (with `ALPHA` and `THOLE`), `BOND`, `DELETE ATOM`, `LONEPAIR relative`,
`LONEPAIR bisector`, and `ANISOTROPY`. It skips only records that carry no
topology here: `GROUP`, `IC`, `DONOR`, `ACCEPTOR`, and
`PATCH FIRST NONE LAST NONE`. Everything else is an error, including `IMPR`,
`DOUBLE`, `TRIPLE`, `CMAP`, other deletions, other lone-pair modes,
continuation lines, duplicate atoms, and duplicate definitions.

The `14bb` header states a charge of -0.017. The patch is in fact neutral:
the donor residue changes by -0.122 and the acceptor residue by +0.122. The
header value is parsed and not used; every patched chain must be neutral to
1e-8.

### Operations after patching

The streams set `AUTOGENERATE ANGLES DIHEDRALS DRUDE PATCH`. The builder
therefore performs these steps after applying `14bb` along the chain:

1. Atoms keep the stream order. A patch `ATOM` that names an existing atom
   replaces it in place; new atoms are appended to their residue.
2. `14bb` lists `dele atom 2DO1`, the Drude of the deleted O1. Drude particles
   are generated after patching, so this deletion is accepted only when the
   same patch application deletes that polarizable parent.
3. One Drude particle (`D` + parent name, type `DRUD`, mass 0.4) follows each
   polarizable atom. Its charge is `-sqrt(2 k_D |alpha| / CCELEC)`; the parent
   keeps the remaining charge and mass. `k_D` is read from the stream's
   `X DRUD` bond (500 kcal/mol/A^2).
4. Angles and proper dihedrals are generated over real atoms only. Lone pairs
   and Drude particles never enter one. The streams define no angle or
   dihedral parameters for `LPD` or `DRUD`, which is consistent with this.
5. Every real-atom bond, angle, and dihedral must have a parameter in the
   streams, or the build is refused.
6. Anisotropy constants follow CHARMM: `A33 = 3 - A11 - A22`,
   `K33 = k_D/A33 - k_D`, `K11 = k_D/A11 - k_D - K33`, and likewise `K22`.
7. Thole pairs are the 1-2 and 1-3 polarizable pairs, with
   `a_ij = a_i + a_j`. NAMD derives them from the PSF; `thole_screened_pairs`
   provides the list for an engine that needs it.

### NAMD input requirements

Two requirements were found by running NAMD 3.0.3 on a native finite bundle
during review. Both are enforced by the writer and by unit tests.

1. Bond order. Every parent--Drude bond must follow every ordinary and
   lone-pair bond in `!NBOND`, across the whole system. With the Drude bonds
   mixed in, NAMD stops with `UNABLE TO FIND BOND PARAMETERS FOR CD31HA DRUD`.
   With them last it reports the lone-pair and Drude bonds as ignored and
   loads the structure.
2. `rigidBonds`. With `drude on` NAMD assumes the SWM4 water model and refuses
   to start unless `rigidBonds` is `water` or `all`. The input uses
   `rigidBonds water` with `useSettle on`, which leaves every cellulose bond
   flexible.

`tests/unit/test_drude.py::test_namd_reads_and_runs_the_native_finite_bundle`
runs a short real NAMD job when `CELLULOSE_FORGE_DRUDE_TOPPAR` is set and
`namd3` is available.

### Position restraints

Finite runnable systems accept the shared `--restrain-chains`,
`--restrain-layers` and `--restraint-force` interface. All real heavy atoms in
the selected chains are restrained by default; `--restrain-atom-names` may
limit them using canonical mapped-template names. Hydrogens, Drude particles
and off-centre lone pairs remain free. OpenMM represents the three Cartesian force
constants directly with a `CustomExternalForce`. NAMD uses its legacy
constraint-reference PDB; because that format stores one scalar coefficient,
the nonzero constants on all active axes must be equal. Other combinations
fail before output is written. The PDB contains a unit mask, and the exact
converted coefficient is written as `constraintScaling`; this avoids rounding
small forces in the fixed-width PDB field. NAMD's legacy potential is
`U = k d^2`, so conversion from the public `U = 1/2 K d^2` convention includes
the factor of one half. Construction-only periodic systems cannot apply
restraints because no runnable engine input exists.

### Coulomb constant

CHARMM derives the Drude charge with `CCELEC = 332.0716`; psfgen uses NAMD's
332.0636. The two give Drude charges that differ by about 1.2e-5 relative.
The writer uses the CHARMM value, because the acceptance reference is a
CHARMM-generated PSF, and records the value in the PSF title. In NAMD this
makes each effective polarizability about 2.4e-5 smaller than the stream
value. This choice should be revisited when the CHARMM reference is available.

### Chain boundaries

Both open and cyclic chains are built. Cyclic chains apply `14bb` from the
last residue back to the first and are used for chain-periodic crystals,
slabs, and fibrils. Residue 1 is the reducing end and links run O4 -> C1, as
in the rest of CelluloseForge.

Known limitation: in a cyclic chain, the three lone pairs on the closing
linkage oxygen have one host in the neighbouring periodic image. The writer
places them from nearest images. OpenMM virtual sites do not apply periodic
boundaries, and OpenMM moved exactly those three sites per chain (by up to
0.67 A) when it rebuilt them. NAMD 3.0.3 was tested on a periodic DP8 I-beta
slab: it reported large distances between the cross-boundary lone-pair hosts
and stopped before energy evaluation with a low-global-exclusion-count error.
Periodic-chain output is therefore construction-only: the command line refuses
periodic Drude chains unless `--drude-construction-only` is given, and then
writes the PSF and coordinates without a NAMD run input. Periodic execution
needs periodic-aware lone-pair handling in the engine or a separately validated
topology representation.

## Evidence so far

Checked on the official streams:

- DP8 finite: 268 explicit sites, 275 bonds, 97 lone pairs, 41 anisotropy
  records, 89 Drude particles, 329 angles, 528 dihedrals.
- DP8 periodic: 264, 272, 96, 40, 88, 328, 528.
- Against a psfgen 2.0 DP2 PSF (built with NAMD's Coulomb constant): all 93
  atom records per chain agree in type, charge, mass, particle kind,
  polarizability, and Thole factor; all 94 bonds agree; anisotropy constants
  agree. The psfgen file lacks the 4 angles and 10 dihedrals that span the
  glycosidic link, because psfgen does not regenerate them after the patch.
  The native topology includes them.
- OpenMM 8.6.1 reads a native finite I-beta slab PSF with the two streams and
  assigns every parameter. It finds the same Drude, lone-pair, and anisotropy
  counts and the same number of Thole pairs (179 per DP6 chain), and its
  rebuilt lone-pair positions match the written ones to 1e-14 A.
- NAMD 3.0.3 refuses a periodic DP8 I-beta slab before energy evaluation after
  reporting large distances between the closing-linkage lone-pair hosts. This
  is the evidence for keeping periodic chains construction-only.

psfgen and OpenMM are consumers of the format, not independent generators of
the topology. They do not replace the CHARMM reference below.

## Required CHARMM references

The native writer is not accepted until it matches CHARMM-generated
references, produced once and committed as fixtures:

- systems: BGLC monomer, DP2, DP3 (the smallest chain with a residue that is
  both donor and acceptor), DP4 finite, and cyclic DP2 and DP4;
- for each: the full PSF and coordinates including Drude particles and lone
  pairs, and CHARMM energies by term and forces on those coordinates;
- comparison by record, keyed by segment, residue, and atom name: types,
  charges, masses, polarizabilities, Thole factors, bonds, angles, dihedrals,
  lone-pair records, and anisotropy records.

Three points can only be settled by that reference: the atom order CHARMM
writes for patch-added atoms, the Coulomb-constant choice above, and the
exclusion of lone pairs and Drude particles from generated angles and
dihedrals.

## Geometry scope

- periodic bulk crystals: cellulose I-alpha, I-beta, II, and III-I
  (periodic chains required, so construction-only for now);
- finite slabs for all four allomorphs, with a NAMD run input;
- finite fibrils for the existing literature-defined I-alpha, I-beta, and
  cellulose II cross-sections, with a NAMD run input;
- chain-periodic slabs and fibrils, construction-only for now;
- no new cellulose III-I fibril preset without a published cross-section.

## Nonbonded protocol

The streams specify `vswitch` with a 10--12 A switch and `e14fac 1.0`. The
NAMD input therefore uses potential switching (`vdwForceSwitching off`),
`switchdist 10`, `cutoff 12`, and unscaled 1-4 interactions.

The long-range Lennard-Jones correction currently follows the project default
used for the other non-CHARMM36 force fields: on for homogeneous periodic
crystals, off for slabs and fibrils, where an isotropic correction is not
valid across the vacuum. This is a placeholder policy, not a validated one.
It has not been checked against a CHARMM or NAMD reference or against the
parameterization papers, and the build report says so.

`drudeBondLen` is left at NAMD's documented 0.25 A default with
`drudeHardWall on`. This confirms the input follows NAMD's hard-wall default;
it does not replace production validation of the resulting Drude excursions.

## Engine scope

| Engine | Status | Requirement before enabling |
| --- | --- | --- |
| NAMD | finite regression and matched OpenMM single points pass; periodic construction-only | production validation; engine support for periodic lone-pair hosts |
| OpenMM | finite regression and matched NAMD single points pass; periodic construction-only | production validation; a solution for periodic-chain lone pairs |
| GROMACS | blocked upstream | general anisotropic polarization and Thole screening in mainline |
| LAMMPS | blocked upstream | anisotropic polarization and general lone-pair sites |
| Amber | unsupported | no verified CHARMM Drude execution route |

An engine is not exposed merely because it has a generic oscillator or shell
model. It must pass term-by-term single-point comparison on identical
coordinates against the NAMD/CHARMM Drude reference. NAMD and OpenMM are
exposed for finite chains; periodic output remains construction-only.

OpenMM and NAMD use the same 10--12 A switching interval, but their switching
functions are not identical. Their Lennard-Jones energies therefore require an
explicit convention correction or a documented tolerance in any future
cross-engine comparison.

## Energy schema

Schema 1 of the energy report has nine terms and no slot for polarization.
The matched validator therefore preserves engine-native bookkeeping instead
of forcing unlike terms into the same category. NAMD includes Drude spring and
anisotropy energy in `BOND` and Thole screening in `ELECT`; its `DRUDEBOND`
column is a thermostat temperature, not an energy. OpenMM includes spring,
anisotropy, and screened-pair contributions in `DrudeForce`.

`validate_drude_cross_engine.py` evaluates identical native PSF parameters,
full-precision binary coordinates, and boxes with validation-only matched
nonbonded settings. A zero-length parent--Drude vector is not used because
NAMD can initialize it nondeterministically even for `run 0`. Every Drude is
instead displaced deterministically by 0.01 and 0.02 A on cycling Cartesian
axes. The validator compares ordinary bonded terms, NAMD `BOND` against
OpenMM covalent bond plus spring/anisotropy energy, NAMD nonbonded energy
against OpenMM nonbonded plus Thole energy, total potential, and the 0.01 to
0.02 A polarization response.

Spring and anisotropy responses agree within 2.5 x 10^-5 relative. The total
polarization-response difference is 3.9 x 10^-5 to 2.5 x 10^-4 relative for
the three slabs. The remaining difference is small but currently
uncharacterized; it is recorded explicitly and gated against the response
itself at 1 x 10^-3, not hidden by scaling the tolerance to the much larger
absolute potential. These checks validate engine translation, not the force
field's accuracy for crystalline cellulose.

## Validation gates

1. Parser fixtures for every selected directive and fail-closed unknowns
   (done, `tests/unit/test_drude_topology.py`).
2. Topology equality against the CHARMM references listed above.
3. Single-point energies by term, in two states: every particle at identical
   positions, including Drude particles and lone pairs; and Drude particles
   relaxed with atoms fixed, comparing energy and total dipole. Per-atom
   forces are compared as well as energies.
4. Minimization, then at least 100 ps of NVT and, for bulk crystals, NPT.
   Report warm and Drude temperatures separately (Drude near 1 K), the
   maximum Drude displacement, the number of hard-wall events (expected
   zero), and cell stability. A 1 ps run is a smoke test only.
5. Geometry matrix across every supported allomorph, shape, and chain
   boundary.
6. Clean-tree provenance for the final validation results.

Force-field accuracy and material properties are outside this implementation
gate and require a separate scientific benchmark.
