# CelluloseForge Project Plan

## Project objective

Develop a modern Python package that constructs crystalline cellulose structures and produces simulation-ready GROMACS systems parameterized with either CHARMM36 or GLYCAM.

Version 1 will support:

- Cellulose I-alpha, I-beta, II, and III-I.
- Parallelepiped crystallites.
- Supported 36-chain elementary fibrils.
- Supported monolayers.
- Finite cellulose chains.
- Covalently periodic cellulose chains.
- CHARMM36 parameters represented correctly for GROMACS.
- GLYCAM parameters represented correctly for GROMACS.
- GROMACS `.gro`, `.top`, and `.itp` output.
- Optional PDB coordinates and machine-readable build reports.

Version 1 will not include functionalized cellulose, OPLS parameterization, native CHARMM/PSF output, or native Amber `prmtop`/`rst7` output.

The original Cellulose-Builder and doGlycans packages will be retained as GPLv3-licensed reference implementations and validation material. New Python code will remain separated from the legacy implementations.

## Guiding principles

1. Crystallographic geometry must remain independent of force-field assignment.
2. Experimentally derived coordinates must not be replaced by idealized force-field coordinates.
3. Molecular connectivity must be explicit before a force-field adapter is applied.
4. CHARMM36 and GLYCAM must share the same crystallographic and molecular model.
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
- Select precise CHARMM36 and GLYCAM parameter versions.
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
```

`SpaceGroupOperation` represents a crystallographic transformation of fractional coordinates:

```text
x' = R x + t
```

where `R` is a rotation matrix and `t` is a fractional translation vector.

### Tasks

- Import unit-cell lengths and angles.
- Import asymmetric-unit fractional coordinates.
- Preserve atom ordering, labels, elements, and source provenance.
- Encode the P1 and P2(1) symmetry operations used by the supported allomorphs.
- Store allomorph definitions in a human-readable structured format.
- Validate coordinate dimensions and atom labels.
- Record the primary crystallographic publication for every dataset.
- Compare imported data with the legacy input files.

### Completion criteria

- All four allomorph definitions load and validate.
- Every atom site has a name, element, and fractional coordinate.
- Calculated unit-cell volumes and vectors match reference values.
- No force-field concepts are present in the crystallographic models.

---

## Phase 3 - Crystallographic geometry engine

### Purpose

Generate complete Cartesian crystal geometries without assigning molecular mechanics parameters.

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
- Coordinate comparisons with legacy reference structures.

### Completion criteria

- All four allomorphs reproduce their reference geometry within documented tolerances.
- Monoclinic and triclinic cell vectors are correct.
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

## Phase 5 - CHARMM36-to-GROMACS backend

### Purpose

Assign CHARMM36 carbohydrate topology and produce a complete GROMACS system without requiring psfgen.

### Tasks

- Represent the CHARMM36 `BGLC` residue template as structured data.
- Implement the CHARMM `14bb` beta(1->4) glycosidic patch.
- Assign CHARMM atom names, atom types, masses, and partial charges.
- Apply finite-terminal treatment.
- Apply periodic glycosidic linkages.
- Generate bonds, angles, proper dihedrals, and explicit impropers.
- Implement CHARMM-compatible GROMACS defaults and exclusions.
- Write `.gro`, `.itp`, `.top`, and optional PDB files.
- Verify complete parameter coverage against the selected CHARMM36 files.
- Record residue templates, patches, and force-field versions in the build report.

### Tests

- `14bb` deletion and atom-modification tests.
- Patch charge-conservation tests.
- Atom-name and atom-type assignment tests.
- Connectivity and interaction-count tests.
- Parameter-coverage tests.
- Comparisons with legacy psfgen-generated topology.
- `gmx grompp` integration tests.
- Small-system energy minimization tests.

### Completion criteria

- GROMACS accepts all supported CHARMM36 structures.
- Every required interaction has a parameter.
- Total charge and connectivity are correct.
- Small reference systems minimize without loss of cellulose stereochemistry.

---

## Phase 6 - GLYCAM-to-GROMACS backend

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

## Phase 7 - Crystal shapes and periodic options

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

## Phase 8 - Scientific validation campaign

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
- Documented comparison of CHARMM36 and GLYCAM behavior.

### Minimum validation matrix

| Allomorph | Finite chains | Periodic chains | CHARMM36 | GLYCAM |
|---|---:|---:|---:|---:|
| I-alpha | Yes | Yes | Yes | Yes |
| I-beta | Yes | Yes | Yes | Yes |
| II | Yes | Yes | Yes | Yes |
| III-I | Yes | Yes | Yes | Yes |

### Completion criteria

- All validation procedures are reproducible from documented commands.
- Reference outputs are versioned and archived with checksums.
- Expected force-field differences are distinguished from implementation errors.

---

## Phase 9 - CLI, documentation, and usability

### Purpose

Turn the validated library into a usable research tool.

### Proposed CLI

```bash
cellulose-forge build \
  --allomorph I-beta \
  --shape crystal \
  --cells 4 4 5 \
  --force-field charmm36 \
  --chain-ends finite \
  --output ibeta_4x4x5
```

```bash
cellulose-forge build \
  --allomorph II \
  --shape fibril \
  --cellobiose 20 \
  --force-field glycam \
  --output cellulose_ii_fibril
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

## Phase 10 - Version 1 release and JOSS preparation

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
5. **M4 - CHARMM36 I-beta:** an I-beta CHARMM36 system passes `gmx grompp`.
6. **M5 - CHARMM36 complete:** every supported allomorph passes CHARMM36 validation.
7. **M6 - GLYCAM I-beta:** an I-beta GLYCAM system passes `gmx grompp`.
8. **M7 - GLYCAM complete:** every supported allomorph passes GLYCAM validation.
9. **M8 - Shapes and periodicity:** fibrils, monolayers, and periodic chains are supported.
10. **M9 - Validation:** the complete scientific validation matrix passes.
11. **M10 - Release:** documentation, Version 1.0, archival, and JOSS materials are ready.

## First vertical prototype

Before generalizing the package, implement one complete narrow workflow:

```text
Cellulose I-beta
  -> 2 x 2 x 1 crystallite
  -> explicit chains and glucose residues
  -> finite beta(1->4)-linked chains
  -> CHARMM36 assignment
  -> GROMACS GRO/TOP/ITP
  -> gmx grompp validation
  -> energy minimization
```

This prototype should expose the major crystallographic, topology, and GROMACS-integration issues early. Generalization to other allomorphs and GLYCAM begins only after the prototype is validated.

## Deferred work

The following features are explicitly deferred until after Version 1:

- Functionalized cellulose.
- Antechamber integration.
- OPLS parameterization.
- Native CHARMM PSF output.
- Native Amber `prmtop` and `rst7` output.
- Protein, lipid, and arbitrary branched-glycan construction.
- A general CHARMM RTF/STR parser.
- A general doGlycans-compatible sequence language.
