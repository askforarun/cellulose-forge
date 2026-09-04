# CelluloseForge Project Plan

## Project objective

Develop a modern Python package that constructs crystalline cellulose structures and produces simulation-ready GROMACS systems parameterized with GLYCAM.

Version 1 will support:

- Cellulose I-alpha, I-beta, II, and III-I.
- Parallelepiped crystallites.
- Surface-oriented slabs exposing user-selected crystallographic faces.
- Fibrils and nanocrystallites with selectable, explicitly reported bounding faces.
- Supported 36-chain elementary fibril presets.
- Supported monolayers.
- Finite cellulose chains.
- Covalently periodic cellulose chains.
- GLYCAM parameters represented correctly for GROMACS.
- GROMACS `.gro`, `.top`, and `.itp` output.
- Optional PDB coordinates and machine-readable build reports.

Version 1 will not include functionalized cellulose, CHARMM36 or OPLS parameterization, native CHARMM/PSF output, or native Amber `prmtop`/`rst7` output.

The original Cellulose-Builder and doGlycans packages will be retained as GPLv3-licensed reference implementations and validation material. New Python code will remain separated from the legacy implementations.

## Guiding principles

1. Crystallographic geometry must remain independent of force-field assignment.
2. Experimentally derived coordinates must not be replaced by idealized force-field coordinates.
3. Molecular connectivity must be explicit before a force-field adapter is applied.
4. GLYCAM assignment must consume the force-field-independent crystallographic and molecular model.
5. An output is not considered supported until `gmx grompp` and the relevant validation tests pass.
6. Scientific provenance, generated coordinates, applied patches, and force-field versions must be recorded.

---

## Phase 0 - Scientific specification

### Purpose

Define the scientific meaning and boundaries of Version 1 before implementation.

### Tasks

- Document supported allomorphs and crystal shapes.
- Define how replication arguments map to crystallographic axes.
- Define finite reducing and nonreducing chain ends.
- Define covalent periodicity along the cellulose-chain direction.
- Define transverse periodic-compatibility options.
- Define origin and center chains for cellulose I-beta and II.
- Define cellulose II antiparallel-chain connectivity.
- Select the precise GLYCAM parameter release.
- Define how external force-field data is discovered and versioned.
- Define GROMACS output contents and naming conventions.
- Establish coordinate, charge, and energy tolerances.
- Record features excluded from Version 1.

### Deliverables

- `docs/specification.md`
- `docs/scientific-conventions.md`
- `docs/force-field-policy.md`
- `docs/output-formats.md`

### Completion criteria

- Every proposed CLI option has an unambiguous scientific meaning.
- Finite ends and periodic bonds are fully specified.
- Force-field versions and provenance policies are documented.
- Version 1 scope is frozen.

---

## Phase 1 - Repository and package foundation

### Purpose

Establish a maintainable open-source Python project suitable for eventual JOSS submission.

### Proposed structure

```text
cellulose-forge/
|-- src/
|   `-- cellulose_forge/
|-- tests/
|   |-- unit/
|   |-- integration/
|   `-- reference/
|-- data/
|   |-- crystallography/
|   `-- schemas/
|-- legacy/
|   |-- cellulose-builder-2013/
|   `-- doglycans-2017/
|-- validation/
|-- examples/
|-- docs/
|-- paper/
|-- pyproject.toml
|-- LICENSE
|-- README.md
|-- CITATION.cff
`-- CONTRIBUTING.md
```

### Tasks

- Adopt GPLv3 and retain all legacy copyright and license notices.
- Create a `src`-layout Python package.
- Set the minimum supported Python version.
- Configure pytest, Ruff, type checking, and documentation generation.
- Add continuous integration on Linux and macOS.
- Establish semantic versioning and a changelog.
- Add contribution, support, and governance documentation.
- Record architectural decisions in the public repository.
- Maintain an AI-assistance log for the eventual JOSS disclosure.

### Completion criteria

- The package installs in a clean environment.
- Tests and quality checks run through documented commands.
- Continuous integration passes.
- Legacy code is isolated and clearly attributed.

---

## Phase 2 - Crystallographic data model

### Purpose

Represent the published crystallographic definitions as validated, traceable data.

### Core models

```text
UnitCell
SpaceGroupOperation
AsymmetricUnit
AtomSite
AllomorphDefinition
CrystalBasis
MillerIndex
CrystalPlaneDefinition
FacetDefinition
FibrilPreset
```

`SpaceGroupOperation` represents a crystallographic transformation of fractional coordinates:

```text
x' = R x + t
```

where `R` is a rotation matrix and `t` is a fractional translation vector.

`MillerIndex` stores an exact signed, nonzero integer triplet `(h, k, l)`. Indices must not be reduced automatically: for example, `(2, 0, 0)` and `(1, 0, 0)` have parallel normals but different interplanar spacings. Likewise, `(h, k, l)` and `(-h, -k, -l)` identify opposite oriented faces. Parallelism, orientation, and exact plane identity must therefore be represented as separate comparisons.

Fibril presets are declarative crystallographic data. They identify bounding faces, chain arrangements, applicable allomorphs, labels, and literature provenance; they do not contain force-field atom types or generated Cartesian coordinates. A slab requires no surface preset: its exposed face is supplied directly as a Miller index.

### Tasks

- Import unit-cell lengths and angles.
- Import asymmetric-unit fractional coordinates.
- Preserve atom ordering, labels, elements, and source provenance.
- Encode the P1 and P2(1) symmetry operations used by the supported allomorphs.
- Store allomorph definitions in a human-readable structured format.
- Define a validated signed Miller-index representation and unambiguous exact serialization without reducing index magnitude or changing sign.
- Define crystallographic-plane records containing a Miller index, human-readable label, and provenance.
- Define a simple slab request schema containing the allomorph, Miller index, layers or thickness, lateral dimensions, and optional vacuum spacing.
- Define fibril and nanocrystallite presets as ordered sets of bounding facets and chain arrangements rather than hard-coded coordinate masks.
- Transcribe the I-alpha, I-beta, and nanocrystallite face sets used by the published reference models, including negative indices and source citations.
- Distinguish a crystallographic face family from an oriented face without exposing a separate surface-termination option.
- Add schemas for user-supplied arbitrary faces and multi-facet cross-sections.
- Validate coordinate dimensions and atom labels.
- Validate that Miller indices are integral and not `(0, 0, 0)`, preset identifiers are unique, referenced allomorphs exist, and facet lists contain no unintended duplicates.
- Round-trip all crystallographic and surface definitions through their human-readable storage format without losing signs, ordering, units, or provenance.
- Record the primary crystallographic publication for every dataset.
- Compare imported data with the legacy input files.

### Tests

- Miller-index validation, magnitude/sign preservation, equality, parallelism, and serialization tests.
- Tests distinguishing exact plane indices, parallel orientations, and opposite oriented faces, including `(1, 1, 0)`, `(2, 2, 0)`, and `(-1, -1, 0)`.
- Plane-label and unit round-trip tests.
- Schema tests for Miller-index-driven slab requests and ordered multi-facet cross-sections.
- Reference-preset tests asserting the complete ordered face lists for I-alpha, I-beta, and the published nanocrystallite models.
- Provenance and unknown-allomorph validation tests.

### Completion criteria

- All four allomorph definitions load and validate.
- Every atom site has a name, element, and fractional coordinate.
- Calculated unit-cell volumes and vectors match reference values.
- Arbitrary valid signed Miller indices load, validate, compare, and serialize without losing magnitude or face orientation.
- All published fibril and nanocrystallite presets load with complete face lists and traceable citations.
- Invalid zero indices, malformed facet definitions, duplicate preset identifiers, and unknown allomorph references fail with clear errors.
- No force-field concepts are present in the crystallographic models.

---

## Phase 3 - Crystallographic geometry engine

### Purpose

Generate complete Cartesian bulk crystals, surface slabs, fibrils, and nanocrystallites without assigning molecular mechanics parameters.

### Geometry models

```text
BulkCrystal
Slab
FibrilCrossSection
GeneratedSurface
```

Phase III consumes the signed face and preset definitions established in Phase II. The geometry engine must derive each plane normal and interplanar spacing from the reciprocal lattice, so face support is not restricted to orthogonal unit cells or to a hard-coded list of orientations.

### Shape-selection contract

- A slab is selected with `--shape slab --face H K L`; the user supplies one Miller index.
- A fibril is selected with `--shape fibril --fibril-model NAME`; the named model supplies the complete cross-section, chain arrangement, and bounding-face set.
- The initial fibril model is `36-chain`, reproduced from the literature independently for each allomorph for which adequate crystallographic evidence is available.
- `36-chain` denotes the published discrete chain arrangement; it must not be interpreted or generated as an ideal regular hexagon.
- `--cellobiose-units N` controls fibril length using the original Cellulose-Builder convention. One cellobiose unit corresponds to two anhydroglucose residues along each chain.
- Users do not need to enter the individual faces of a named fibril model. Every face used by the model must be listed in the build report.

### Implementation order

1. Cellulose I-beta.
2. Cellulose II.
3. Cellulose III-I.
4. Cellulose I-alpha.

Cellulose I-beta is the first vertical target because it exercises P2(1) symmetry and origin/center chains and has strong legacy references.

### Tasks

- Implement general triclinic cell-vector construction.
- Convert fractional coordinates to Cartesian coordinates.
- Apply space-group symmetry operations.
- Normalize or wrap fractional coordinates when appropriate.
- Replicate unit cells along crystallographic axes.
- Compute reciprocal-lattice vectors, plane normals, and interplanar spacings for general `(h k l)` faces.
- Generate a slab from a user-supplied Miller index, number of molecular layers or thickness, lateral dimensions, and optional vacuum spacing.
- Use a documented, deterministic whole-layer cutting convention internally; Version 1 will not expose a separate termination parameter.
- Generate fibrils and nanocrystallites from named models by applying their complete bounding-face sets and chain arrangements.
- Avoid split residues and dangling fragments at cut boundaries by applying a documented whole-chain/whole-residue inclusion policy.
- Provide published fibril cross-section models; arbitrary individual Miller-index input remains available through slab generation.
- Include reference presets covering the faces illustrated for cellulose I-alpha, I-beta, and cellulose nanocrystallites, including opposite signed faces where required to close the cross-section:
  - I-alpha: `(1, 1, 0)`, `(0, 1, 0)`, `(1, 0, 0)`, `(-1, 1, 0)`, `(1, -1, 0)`, and `(-1, -1, 0)`.
  - I-beta: `(2, 0, 0)`, `(1, 1, 0)`, `(1, -1, 0)`, `(0, 1, 0)`, `(0, -1, 0)`, and the additional opposite faces needed to close the cross-section.
  - Nanocrystallite reference cross-section: `(0, 2, 0)`, `(1, 2, 0)`, `(1, 0, 0)`, `(1, -1, 0)`, `(0, -2, 0)`, `(-1, 0, 0)`, and `(-1, 1, 0)`.
- Express face labels canonically in machine-readable build reports while accepting conventional overbar/signed input notation.
- Preserve deterministic atom ordering.
- Generate output box vectors.
- Write diagnostic XYZ and PDB coordinates.
- Compare generated structures with the legacy builder.

### Tests

- Unit-cell vector and volume tests.
- Fractional-to-Cartesian round-trip tests.
- Space-group operation tests.
- Symmetry-equivalent atom tests.
- Atom-count tests.
- Translation-invariance tests.
- Reciprocal-lattice plane-normal and `d(hkl)` tests for orthogonal and non-orthogonal cells.
- Slab orientation, layer count, thickness, vacuum, exposed-face, and deterministic-boundary tests.
- Tests for both members of each opposite-face pair, including negative Miller indices.
- Fibril cross-section tests for every reference face and all published presets.
- Boundary tests ensuring that cutting does not silently create partial residues or chains.
- Determinism tests for equivalent face specifications and repeated builds.
- Coordinate comparisons with legacy reference structures.

### Completion criteria

- All four allomorphs reproduce their reference geometry within documented tolerances.
- Monoclinic and triclinic cell vectors are correct.
- A slab can be generated for any valid user-supplied `(h k l)` face for every supported allomorph.
- Fibril and nanocrystallite cross-sections reproduce the documented reference face sets and dimensions within tolerance.
- Every generated surface reports its requested Miller index, Cartesian normal, interplanar spacing, realized thickness, area, and applied deterministic boundary convention.
- No supported surface operation leaves an unreported partial residue or chemically ambiguous cut.
- Geometry generation requires no VMD, Octave, Perl, sed, psfgen, or Open Babel.

---

## Phase 4 - Molecular graph and cellulose chemistry

### Purpose

Convert crystallographic atoms into explicit glucose residues, cellulose chains, and chemical bonds.

### Core models

```text
Atom
Bond
Residue
PolymerChain
PeriodicBond
MolecularSystem
```

### Tasks

- Assign atoms to anhydroglucose residues.
- Assign residues to cellulose chains.
- Number atoms, residues, and chains deterministically.
- Identify origin and center chains where applicable.
- Represent parallel and antiparallel chains explicitly.
- Add intraglucose covalent bonds.
- Add beta(1->4) glycosidic bonds.
- Implement finite reducing and nonreducing ends.
- Implement glycosidic bonds across periodic boundaries.
- Track whether each coordinate is experimental, symmetry-generated, translated, or geometrically generated.

### Special cases

- Antiparallel origin and center chains in cellulose II.
- The different crystallographic-axis arrangement in cellulose I-alpha.
- Missing hydroxyl hydrogen coordinates in cellulose II.
- Terminal H and OH groups for finite chains.
- Removal of terminal caps for covalently periodic chains.

### Completion criteria

- Every atom belongs to one residue and one chain.
- Every chain has an explicit and correct residue order.
- No atom has an invalid valence.
- Finite chains have the expected chemical formula.
- Periodic chains contain the correct cross-boundary glycosidic bonds.

---

## Phase 5 - GLYCAM-to-GROMACS backend

### Purpose

Assign GLYCAM carbohydrate topology and produce a complete GROMACS system without copying the doGlycans implementation.

### Tasks

- Select and document the GLYCAM parameter release.
- Define terminal and internal beta-D-glucose residue naming.
- Assign GLYCAM atom names, atom types, masses, and partial charges.
- Implement beta(1->4) linkage behavior.
- Implement finite reducing and nonreducing ends.
- Implement covalently periodic chains.
- Implement GLYCAM electrostatic and Lennard-Jones 1-4 scaling.
- Convert GLYCAM bonded terms into correct GROMACS functional forms.
- Write `.gro`, `.itp`, `.top`, and optional PDB files.
- Verify complete parameter coverage.
- Compare isolated finite chains against independent doGlycans or AmberTools references.

### Tests

- GLYCAM residue-name and atom-name tests.
- Residue and system charge tests.
- Linkage topology tests.
- 1-4 pair and scaling tests.
- Proper and improper torsion tests.
- Parameter-coverage tests.
- Reference-chain comparison tests.
- `gmx grompp` integration tests.
- Small-system energy minimization tests.

### Completion criteria

- GROMACS accepts all supported GLYCAM structures.
- No topology-related missing-parameter warnings remain.
- Charges, connectivity, and scaling match documented GLYCAM conventions.
- An isolated cellulose chain agrees with an independent reference.

---

## Phase 6 - Crystal shapes and periodic options

### Purpose

Implement the full set of published cellulose structures after basic geometry and force-field assignment are stable.

### Tasks

- General parallelepiped crystallites.
- I-alpha, I-beta, and II 36-chain elementary fibrils.
- I-alpha monolayers.
- I-beta origin and center monolayers.
- Cellulose II origin and center monolayers.
- Transverse periodic-compatibility options.
- Covalent periodicity along the chain direction.
- Surface and Miller-index metadata.
- Diagnostic chain-selection and cross-section outputs.

### Tests

- Expected chain-count tests.
- Fibril cross-section tests.
- Origin/center membership tests.
- Monolayer-selection tests.
- Periodic edge-matching tests.
- Shape-specific basis-vector tests.

### Completion criteria

- Every advertised shape matches its published or legacy reference.
- Selections and atom ordering are deterministic.
- Requested periodic edges match within numerical tolerance.

---

## Phase 7 - Scientific validation campaign

### Purpose

Demonstrate chemical and physical correctness rather than only valid file syntax.

### Structural validation

- Unit-cell parameters and basis vectors.
- Fractional and Cartesian coordinates.
- Atom, residue, and chain counts.
- Bond lengths, bond angles, and glycosidic torsions.
- Glucose-ring conformations.
- Parallel and antiparallel chain directionality.

### Topological validation

- Atom names, types, masses, and charges.
- Bonded interaction counts.
- Exclusions and 1-4 scaling.
- Proper and improper dihedrals.
- Periodic covalent bonds.
- Complete force-field parameter coverage.

### Simulation validation

- Successful `gmx grompp` preprocessing.
- Restrained and unrestrained energy minimization.
- Periodic crystal minimization.
- Short NVT simulations.
- Preservation of connectivity and stereochemistry.
- Agreement with documented GLYCAM conventions and independent reference systems.

### Minimum validation matrix

| Allomorph | Finite chains | Periodic chains | GLYCAM |
|---|---:|---:|---:|
| I-alpha | Yes | Yes | Yes |
| I-beta | Yes | Yes | Yes |
| II | Yes | Yes | Yes |
| III-I | Yes | Yes | Yes |

### Completion criteria

- All validation procedures are reproducible from documented commands.
- Reference outputs are versioned and archived with checksums.
- Expected structural differences between finite and periodic models are distinguished from implementation errors.

---

## Phase 8 - CLI, documentation, and usability

### Purpose

Turn the validated library into a usable research tool.

### Proposed CLI

Version 1 uses GLYCAM exclusively, so no force-field selection option is required. The exact GLYCAM release and parameter provenance are always recorded in the build report.

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape crystal \
  --cells 4 4 5 \
  --chain-ends finite \
  --output ibeta_4x4x5
```

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape fibril \
  --fibril-model 36-chain \
  --cellobiose-units 20 \
  --output ibeta_36chain_fibril
```

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape slab \
  --face 1 1 0 \
  --layers 8 \
  --lateral-cells 6 6 \
  --vacuum 3.0 \
  --output ibeta_110_slab
```

```bash
cellulose-forge validate ibeta_4x4x5/system.top
```

### Expected output

```text
ibeta_4x4x5/
|-- structure.gro
|-- structure.pdb
|-- cellulose.itp
|-- system.top
|-- cell.json
`-- build_report.json
```

### Documentation

- Installation and dependency management.
- Quick-start tutorial.
- Scientific and crystallographic conventions.
- Allomorph and crystal-shape descriptions.
- Force-field selection and provenance.
- Finite and periodic-chain behavior.
- Output-file descriptions.
- Reproducible GROMACS examples.
- Python API reference.
- Troubleshooting and limitations.

### Completion criteria

- A new user can install the package and build a validated example without reading source code.
- Every CLI option and output field is documented.
- Documentation examples are tested automatically.

---

## Phase 9 - Version 1 release and JOSS preparation

### Purpose

Release stable research software and prepare it for peer review.

### Tasks

- Publish version `1.0.0`.
- Archive the release and validation dataset on Zenodo.
- Finalize `CITATION.cff`.
- Add the JOSS `paper.md` and bibliography.
- Document the statement of need and research significance.
- Compare capabilities with Cellulose-Builder, doGlycans, and related tools.
- Publish reproducible validation results and benchmarks.
- Complete the required AI-use disclosure.
- Ensure an adequate public development history.
- Invite external testing and record user feedback.

### Completion criteria

- The public release is stable and archived with a DOI.
- Continuous integration and documentation builds pass.
- The JOSS checklist is satisfied.
- At least one independent user has reproduced a documented workflow if possible.

---

## Milestones

1. **M0 - Specification:** scientific specification and repository foundation complete.
2. **M1 - I-beta geometry:** a small cellulose I-beta crystal is generated entirely in Python.
3. **M2 - Allomorph geometry:** all four allomorphs reproduce reference coordinates.
4. **M3 - Chemistry:** explicit residues, chains, termini, and beta(1->4) bonds are correct.
5. **M4 - GLYCAM I-beta:** an I-beta GLYCAM system passes `gmx grompp`.
6. **M5 - GLYCAM complete:** every supported allomorph passes GLYCAM validation.
7. **M6 - Shapes and periodicity:** slabs, fibrils, monolayers, and periodic chains are supported.
8. **M7 - Validation:** the complete scientific validation matrix passes.
9. **M8 - Release:** documentation, Version 1.0, archival, and JOSS materials are ready.

## First vertical prototype

Before generalizing the package, implement one complete narrow workflow:

```text
Cellulose I-beta
  -> 2 x 2 x 1 crystallite
  -> explicit chains and glucose residues
  -> finite beta(1->4)-linked chains
  -> GLYCAM assignment
  -> GROMACS GRO/TOP/ITP
  -> gmx grompp validation
  -> energy minimization
```

This prototype should expose the major crystallographic, topology, and GROMACS-integration issues early. Generalization to other allomorphs begins only after the prototype is validated.

## Deferred work

The following features are explicitly deferred until after Version 1:

- Functionalized cellulose.
- Antechamber integration.
- CHARMM36 parameterization.
- OPLS parameterization.
- Native CHARMM PSF output.
- Native Amber `prmtop` and `rst7` output.
- Protein, lipid, and arbitrary branched-glycan construction.
- A general doGlycans-compatible sequence language.
