# ADR 0009: experimental AMOEBA cellulose backend

## Status

Experimental construction milestone. Finite slabs of cellulose I-beta, II and
III-I and literature-defined fibrils of cellulose I-alpha, I-beta and II
can be written as AMOEBA/OpenMM bundles using Tinker-format topology and
parameters read by OpenMM 8.6.1. Force-field accuracy, reference-engine energy
agreement, water, and
production dynamics have not been validated.

## Source and redistribution

The parameter source is the public Welborn-group repository
`Newman_D-Glucose_2026`, commit
`a2adb2a3b21f4df9696fcd13cf0cc21d36ae6c71`, associated with the 2026 AMOEBA
D-glucose and 1->4-linked oligomer work. The repository has no explicit licence
file as audited on 2026-11-20. CelluloseForge therefore does not include,
download, or redistribute its files. Users supply their own checkout.

The initial backend accepts only the cleaned published beta parameter file
`beta_glucuse_1-4.prm` and the same-DP beta oligomer topology. Every accepted
file is pinned by SHA-256. A missing, modified, or unsupported file fails before
output is created.

## Initial scientific scope

- polymer: native cellulose;
- allomorphs: I-beta, II and III-I for slabs; I-alpha, I-beta and II
  for literature-defined fibrils;
- shapes: finite slabs and finite crystallographic fibrils;
- chain lengths: DP 2, 4, 6, 8, and 10, matching the supplied references;
- engine: OpenMM;
- solvent: none.

Periodic chains are excluded because no cyclic AMOEBA topology or multipole
frame reference is supplied. A III-I fibril is excluded because CelluloseForge
has no literature-defined cross-section for that allomorph. Water is excluded
from the construction milestone so that solute topology and engine translation
can be validated first.

I-alpha slabs are temporarily refused. The shared slab builder can produce
close lateral-image contacts for some chain counts because the I-alpha
translation carries a component along the chain direction. The I-alpha
36-chain fibril has transverse vacuum and no periodic seam, so it remains in
scope. The shared geometry correction is intentionally separate from this
force-field backend.

The same published beta-(1->4)-D-glucose atom types, multipoles and bonded
parameters are transferred without refitting to each crystalline allomorph.
This is a construction hypothesis, not evidence that AMOEBA reproduces their
relative stability or crystal properties; those claims require the validation
gates below.

## Construction

The cleaned beta parameter file defines beta-D-glucose and the cellulose
beta-(1->4) linkage, including permanent multipoles, atomic polarizabilities,
mutual induction, buffered 14-7 van der Waals, stretch-bend and torsion terms.

For each chain, CelluloseForge reads the audited published topology with the
requested DP. Its atom types, atom order and connectivity are retained exactly.
The coordinates are replaced by CelluloseForge's crystallography-first finite
chain coordinates through the shared polymer mapper. Repeated chains use
offset copies of the same graph. Terminal O1/HO1 and HO4 coordinates come from
the existing finite-chain mapping contract.

The first implementation requires an orthorhombic simulation box. This covers
the supported finite slab and fibril orientations without rotating the
crystallographic coordinates.

## Engine strategy

OpenMM reads the Tinker-format XYZ and parameter file using
`openmm.app.TinkerFiles`; CelluloseForge does not reproduce the AMOEBA
functional form itself. Tinker is a future reference engine, not a currently
supported output target.
The generated runner lets OpenMM select the fastest available platform rather
than forcing the slow Reference implementation.
The generated dry slab and fibril smoke test is NVT. It intentionally has no
barostat: isotropic pressure coupling would compress the vacuum surrounding
these non-bulk geometries.

The shared chain/layer selectors restrain all real heavy atoms by default and
may be intersected with canonical mapped-template names through
`--restrain-atom-names`. The resulting indices feed an OpenMM
`CustomExternalForce`. Hydrogens remain free, and independent X/Y/Z force
constants are supported. The exact selected chain IDs, layer mapping, force
constants, selection mode, names and resolved count are recorded in
`build_report.json`.

The vacuum rule disables the isotropic vdW long-range correction. OpenMM's
`TinkerFiles` enables the corresponding `AmoebaVdwForce` flag by default, so
the generated runner explicitly turns it off and refuses an unexpected
vdW-force count. This differs from the published oligomer key because the
present systems contain transverse vacuum.

OpenMM 8.6.1 constructs the cleaned published file with these forces:

- `CustomBondForce`;
- `CustomAngleForce`;
- `CustomCompoundBondForce`;
- `PeriodicTorsionForce`;
- `AmoebaVdwForce`;
- `AmoebaMultipoleForce`.

The larger combined `a_glc.prm`/`b_glc.prm` files are intentionally not used.
In OpenMM 8.6.1 their unrelated angle-torsion records trigger an error when the
corresponding in-plane-angle force is absent from the glucose system. The
cleaned files supplied specifically for 1->4 glucose oligomers load without
editing and are the defensible input for this backend.

## Evidence at this milestone

The initial I-beta systems built from the audited external checkout were:

- I-beta (110), 3 layers x 6 chains, finite DP6: 2,322 atoms;
- I-beta 18-chain fibril, finite DP6: 2,322 atoms.

OpenMM 8.6.1 created a complete 2,322-particle AMOEBA System from each bundle.
The slab produced a finite initial potential energy on the CPU platform. These
are construction checks, not force-field validation.

During development, the allomorph expansion constructed a DP2 I-alpha (010)
slab before the shared-builder seam-contact defect was identified. I-alpha
slabs are now refused and are not part of the supported live matrix. The
current matrix builds I-beta (110), II (110) and III-I (110) slabs, plus the
supported I-alpha, I-beta and II fibrils.
OpenMM constructed every System, found exactly one vdW force with its vacuum
long-range correction disabled, and returned a finite potential energy for
each. This verifies parameter coverage and engine construction, but not the
transferred model's scientific accuracy.

## Required validation before expanding scope

1. Install Tinker9 and compare Tinker and OpenMM single-point energies term by
   term for the supplied DP2-DP10 references.
2. Repeat the comparison after substituting crystal coordinates.
3. Check minimization and short NVT stability for dry slabs and fibrils.
4. Add matching AMOEBA water only after dry-system agreement.
5. Compare fixed-charge CHARMM36, CHARMM Drude and AMOEBA using observables
   defined before production simulations.
6. Request redistribution/licensing clarification before considering bundled
   parameters or automated installation.

## Non-goals

This decision does not make AMOEBA a universal topology hub, convert AMOEBA to
GROMACS/NAMD/LAMMPS/Amber, or mix AMOEBA parameters with CHARMM Drude. The two
polarizable models remain independent backends with independent validation.
