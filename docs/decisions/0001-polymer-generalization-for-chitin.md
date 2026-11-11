# 0001: Minimum generalization for chitin crystals

- Status: accepted; unit-cell and bulk construction implemented
- Branch: `feature/chitin`
- Scope: coordinate-only α- and β-chitin unit cells and periodic bulk crystals.
  Force-field output, slabs, fibrils, and amylopectin are out of scope here.

## Context

CelluloseForge already separates crystallographic data from geometry and
force-field mapping (see the architecture figure in the README). The audit
below records where cellulose is nevertheless assumed, so chitin can be added
without changing any cellulose output.

## Audit of cellulose-specific assumptions

| Area | Location | Assumption | Needed for chitin coordinates? |
|---|---|---|---|
| Crystal data schema | `crystallography.py` (`AsymmetricUnit`, `AllomorphDefinition`, JSON) | Already generic: atom labels, residue atom count, chain labels and directions, symmetry operations, chain axis, citation and DOI | No change except a `polymer` field and richer provenance |
| Allomorph lookup | `cli.py` (`--allomorph` choices and the name-to-file map) | Only the four cellulose allomorph names exist | Yes: add `--polymer` and per-polymer allomorph choices |
| Chirality check | `geometry.validate_cellulose_chirality` | Checks only C1 of β-D-glucose using O5, C2, H1 | Yes: generalize to every stereocentre declared in the polymer definition (chitin C1–C5) |
| Degree of polymerization | `geometry.py` (cellobiose wording and even-DP rule) | Two-residue crystallographic repeat named "cellobiose" | Generalize to "crystallographic repeat of *n* residues"; chitin's repeat is also a disaccharide (chitobiose) |
| Slab lattice and fibril presets | `geometry.py` (allomorph-name branches, `_FIBRIL_LAYER_SIZES`, 36-chain presets) | Cellulose-specific lateral repeats and published cellulose fibrils | Not in this phase; reject `--shape slab/fibril` for chitin until defined |
| PDB annotations | `pdb.py` (origin/center remark for I-beta) | Cellulose site naming | Make site labels data-driven |
| Periodic bond and clash checks | `clashes.py` (O4–C1 glycosidic bond) | β(1→4) linkage atoms are O4 and C1 | Chitin is also β(1→4) (O4–C1), so the check can be reused once the linkage atoms come from the definition |
| GLYCAM and CHARMM36 mapping | `glycam*.py`, `charmm36.py` (residues 4GB/0GB/ROH, BGLC, 14bb) | Glucose residues, cellulose termini and patches | Not in this phase; force-field builds must reject non-cellulose polymers |

## Decision

1. **Add a polymer layer to the data files without changing cellulose.**
   Each crystal JSON gains optional fields with cellulose defaults:
   `polymer` (default `cellulose`), `repeat_residues` (2),
   `linkage` (`{"donor": "C1", "acceptor": "O4", "type": "beta-1-4"}`),
   `stereocentres` (cellulose: `C1` only, matching today's check), and
   `site_labels`. Chitin files live in `data/crystallography/chitin/` as
   `alpha.json` and `beta.json`; cellulose files stay where they are.
2. **CLI:** add `--polymer {cellulose,chitin}` (default `cellulose`).
   `--allomorph` choices depend on the polymer. Existing commands produce
   byte-identical output.
3. **Chirality:** replace the single C1 test with a definition-driven check of
   every listed stereocentre. Cellulose keeps exactly its current convention and
   expected sign, so its build reports are unchanged.
4. **Shapes:** chitin supports `unit-cell` and `bulk` only. Slab, fibril, and
   any `--force-field` other than `none` fail with a clear message.
5. **Regression guard:** before any refactor, record SHA-256 hashes of the
   committed cellulose examples and a fixed matrix of coordinate builds, and
   assert them in a test, so "byte-identical" is checked rather than assumed.

## Data requirements (blocking)

Coordinates must come from primary sources, never from figures or secondary
diagrams. For each allomorph the package needs the unit cell, space group,
fractional coordinates of the asymmetric unit (including hydrogens, or a
documented procedure for placing them), chain orientation, experimental
temperature and hydration state, and the licence for redistributing the
coordinates.

Candidate primary sources to verify before use (not yet checked):

- α-chitin: high-resolution X-ray fibre diffraction refinement, e.g.
  Sikorski, Hori & Wada, *Biomacromolecules* (2009).
- β-chitin: anhydrous β-chitin at high resolution, e.g. Nishiyama, Noishiki &
  Wada, *Macromolecules* (2011).
- γ-chitin: excluded unless a complete crystallographic model exists.

The coordinate tables or CIF files from these papers' supplementary material
are required. Bibliographic details and DOIs must be confirmed against the
publisher records before they enter the data files.

## Validation for this phase

Unit-cell matrix and volume; symmetry-expanded atom counts; stereochemistry at
C1–C5 of every residue; glycosidic O4–C1 distances including across periodic
boundaries; minimum nonbonded distances including periodic images; agreement
of generated sites with the source coordinates; deterministic chain and residue
numbering; and the cellulose hash regression test.

## Consequences

Cellulose behaviour is untouched and enforced by the hash test. Chitin becomes
available for coordinate-only bulk builds once verified data are added. The
amylopectin work (graph topology, α linkages, branching) is not constrained by
this decision beyond reusing the polymer field.
