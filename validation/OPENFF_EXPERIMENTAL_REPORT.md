# Experimental OpenFF validation report

## Status

This is a finite- and periodic-chain feasibility study on
`feature/openff-experimental`. It
does not add `--force-field openff` to the public CLI. OpenFF support should
not be claimed in the manuscript or README until the implementation is
independently reviewed, deliberately merged, and given longer stability tests.

## Reproducible environment

The optional environment is separate from the normal CelluloseForge
installation:

```bash
conda env create -f environment-openff.yml
conda activate cellulose-forge-openff
```

It pins Python 3.12, OpenFF Toolkit 0.19.0, Interchange 0.5.5, OpenFF NAGL
0.6.1, OpenFF NAGL models 2026.9.0, OpenMM 8.6.1, AmberTools 26, and ParmEd.
GROMACS is intentionally not installed by this
file and is required externally for the cross-engine matrix. AmberTools and
`tleap` remain necessary for the current coordinate-mapping scaffold, even
when NAGL supplies the partial charges.

## Scope tested

- Force field: OpenFF Sage `openff-2.3.0.offxml`.
- Charge methods:
  - graph-based NAGL model `openff-gnn-am1bcc-1.0.0.pt` (recommended
    experimental route);
  - AmberTools AM1-BCC as a conformer-dependent comparison.
- Chemistry: finite beta-D-cellobiose (DP2) and beta-D-cellohexaose (DP6).
- Coordinates: CelluloseForge crystallographic cellulose coordinates reordered
  through the existing finite-chain mapping.
- Slabs: two layers, four chains per layer, finite DP6 chains, and 12 A of
  longitudinal padding for I-alpha (010), I-beta (110), II (110), and III-I
  (110). Each system contains 8 chains and 1,032 atoms.
- Periodic slabs: the same four reference faces and eight-chain layout, using
  cyclic DP6 chains with 1,008 atoms per system. The current periodic matrix
  uses NAGL charges; periodic AM1-BCC has not been tested.
- Engines: OpenMM and GROMACS. Amber export is not treated as equivalent
  because Amber cannot reproduce Sage's 8 A switching function and substitutes
  a hard cutoff.

Homogeneous-fluid dispersion correction is disabled in both engines for these
vacuum slabs. The comparison uses an OpenMM Ewald tolerance of 1e-6 and
GROMACS `ewald-rtol = 1e-6`, `pme-order = 8`, `fourierspacing = 0.05` nm, and
`DispCorr = no`. A PME scan showed that the earlier order-4 interpolation,
rather than unconverged OpenMM electrostatics, caused the approximately
1 kJ/mol offset. Eighth-order interpolation reduces the largest finite-chain
residual to 2.78e-5 kJ mol-1 atom-1. Across the periodic NAGL matrix, the
largest residual is 3.24e-5 kJ mol-1 atom-1. A 5.0e-5 electrostatics tolerance
allows reasonable headroom for mixed-precision and platform-dependent
rounding.

## Molecular graph and charge provenance

- DP2 graph: PubChem CID 10712 (beta-cellobiose), 45 atoms and 46 bonds.
- DP6 graph: PubChem CID 440948 (beta-D-cellohexaose), 129 atoms and 134 bonds.
- Both graphs contain explicit stereochemistry and hydrogens and both charge
  methods produce a neutral molecule to numerical precision.
- NAGL charges depend on the molecular graph, not a conformer, and complete in
  seconds for DP6.
- OpenFF's AM1-BCC call was made with `use_conformers=None`; OpenFF generated
  an internal charge conformer rather than using the retained crystallographic
  coordinates. DP6 charging took approximately 13 minutes on the test machine.
- Every probe JSON records the graph, crystallographic-coordinate, and charged
  cache SHA-256 hashes plus the toolkit and charging backend. NAGL probes also
  record OpenFF NAGL 0.6.1 and model package 2026.9.0.
- Equivalent crystal chains reuse one cached charged DP6 template; they are
  never charged independently.

## Geometry and stereochemistry checks

Before minimization, every molecular bond in all 32 chains is checked (1,072
bonds per slab). Across the matrix, X-H bonds span 0.9603--0.9966 A and heavy
bonds span 1.3776--1.5581 A. Every chain is reconstructed as a conformer and
compared with the explicit stereochemical reference graph using RDKit/OpenFF.
Existing chiral tags are cleared before stereochemistry is perceived from the
3D coordinates. All 32 comparisons pass for both charge matrices, and each
case also mirrors its first chain as a negative control and correctly rejects
that inverted geometry.

These checks validate the known molecular bond graph. They are not a general
nonbonded-clash detector and do not independently discover a missing bond that
was absent from both the reference graph and the generated topology.

## Four-allomorph results

Both the NAGL and AM1-BCC matrices pass all four cases. GROMACS preprocessing
uses `-maxwarn 0`; every minimization and 1 ps NVT run reaches its final frame
without LINCS, NaN, or non-finite markers. OpenMM also completes minimization
and 1 ps NVT.

| Case | NAGL maximum OpenMM bond change | AM1-BCC maximum OpenMM bond change |
|---|---:|---:|
| I-alpha (010) | 0.098 A | 0.106 A |
| I-beta (110) | 0.109 A | 0.104 A |
| II (110) | 0.094 A | 0.105 A |
| III-I (110) | 0.117 A | 0.109 A |

The largest absolute OpenMM/GROMACS discrepancy per atom is:

| Energy term | NAGL (kJ mol-1 atom-1) | AM1-BCC (kJ mol-1 atom-1) | Allowed |
|---|---:|---:|---:|
| Bonds | 1.43e-5 | 1.43e-5 | 2.0e-5 |
| Angles | 3.94e-6 | 3.94e-6 | 2.0e-5 |
| Torsions | 8.95e-7 | 8.95e-7 | 2.0e-5 |
| vdW | 1.49e-6 | 1.49e-6 | 1.0e-5 |
| vdW 1-4 | 2.04e-6 | 2.15e-6 | 2.0e-5 |
| Electrostatics | 2.78e-5 | 2.67e-5 | 5.0e-5 |
| Electrostatics 1-4 | 1.36e-5 | 9.97e-6 | 2.0e-5 |

All seven required energy categories must be present in both engines; a
missing term is an error rather than being silently omitted.

## Periodic-chain results

The NAGL periodic matrix passes all four allomorphs. The cyclic template has
126 atoms and 132 bonds. Its stereochemical labels are derived independently
from the internal residue of the PubChem finite DP6 reference, rather than
from the crystal coordinates under test. All 960 C1--C5 stereocentres pass,
and local C1 inversion, local C3 inversion, and whole-chain mirroring are
rejected in every case.

Each system has eight O4--C1 bonds crossing the periodic boundary. Minimum-
image checks cover those boundary bonds and all 1,056 bonds before and after
OpenMM minimization and 1 ps NVT. The largest change in any bond is 0.122 A
for I-alpha, 0.104 A for I-beta, 0.108 A for II, and 0.097 A for III-I. Final
boundary-crossing O4--C1 bonds span 1.361--1.475 A across the matrix. The
largest OpenMM/GROMACS point-energy difference is 3.24e-5 kJ mol-1 atom-1,
from electrostatics; all seven required categories pass their tolerances.

All OpenMM minimization and NVT checks use the CPU platform with one thread,
and the platform and properties are recorded in every case. This removes
platform selection and CPU thread scheduling as sources of regenerated-value
drift.

The current periodic dynamics check is OpenMM-only. GROMACS is used for the
independent point-energy comparison, but periodic GROMACS minimization and NVT
remain to be added before making a broad periodic-support claim.

## Automated checks

- Normal CelluloseForge suite on the final branch state: 269 passed, 11
  skipped.
- Both four-allomorph real-engine matrices pass.
- The four-allomorph periodic NAGL matrix passes.
- Ruff and `git diff --check` are run before handoff.
- Existing public behavior is unchanged: this branch adds only an optional
  environment, validation scripts, and result artifacts.

## Known limitations and design decisions

1. Periodic support is experimental and currently uses only NAGL charges.
   Periodic bulk crystals, fibrils, GROMACS dynamics, and AM1-BCC remain to be
   tested; the present evidence covers four reference slab faces.
2. OpenFF's supported `Topology.from_pdb(..., unique_molecules=...)` route does
   not recognize the multi-residue GLYCAM `ROH/4GB/0GB` representation. The
   probe therefore uses deprecated `Molecule.from_pdb_and_smiles`. Production
   code should construct the explicit OpenFF graph directly and attach the
   crystallographic coordinates by a tested atom map.
3. NAGL avoids conformer-dependent charges and is therefore the cleaner
   experimental default. AM1-BCC remains useful comparative evidence but its
   internally generated charge conformer must be disclosed.
4. One-picosecond trajectories are smoke tests, not evidence of long-time
   crystal stability or reproduction of experimental cell parameters.
5. Amber output is not validated because of the switching-function mismatch.
6. No public CLI, normal dependency group, build-report schema, CI job, or
   release documentation has been implemented.

## Files for review

- `environment-openff.yml`
- `validation/probe_openff_cellobiose.py`
- `validation/probe_openff_periodic.py`
- `validation/validate_openff_finite_slab.py`
- `validation/validate_openff_periodic_slab.py`
- `validation/results/openff_cellodextrin_dp2_probe.json`
- `validation/results/openff_cellodextrin_dp6_probe.json`
- `validation/results/openff_cellodextrin_dp6_nagl_probe.json`
- `validation/results/openff_cellodextrin_dp6_periodic_nagl_probe.json`
- `validation/results/openff_finite_slab_matrix.json`
- `validation/results/openff_finite_slab_nagl_matrix.json`
- `validation/results/openff_periodic_slab_nagl_matrix.json`

The corresponding charged-template JSON files are cached evidence and include
the exact partial charges used by each matrix.

## Questions for independent review

1. Is the corrected coordinate-derived stereochemistry check, including its
   mirrored-chain negative control, sufficient evidence for every allomorph?
2. Is NAGL an acceptable default for experimental finite-chain support, with
   AM1-BCC retained only as a documented alternative?
3. Is the measured 5.0e-5 kJ mol-1 atom-1 PME tolerance appropriate for the
   documented eighth-order, 0.05 nm GROMACS calculation?
4. Should a production implementation construct the OpenFF graph directly,
   rather than retain the deprecated PDB-plus-SMILES loader?
5. Should experimental exposure wait for periodic bulk and fibril validation,
   GROMACS periodic dynamics, and longer stability runs?
