## v0.4.1 (2026-06-01)

### Fix

- __version__

## v0.4.0 (2026-06-01)

### Feat

- export ReferenceBuilder/ReferenceSpec/RepeatFeature; bump 0.3.0
- reference_and_documents() paired strategy
- reference-consistent documents() with violation labels
- references() strategy drawing ReferenceSpec with planted repeats
- ReferenceSpec.write bgzipped+faidx FASTA
- ReferenceBuilder/ReferenceSpec with tandem-repeat provenance
- general per-variant labels carried into GroundTruth

### Refactor

- share draw_ref_alt between Reference and ReferenceSpec

## v0.2.1 (2026-05-31)

### Fix

- declare hypothesis as a runtime dependency

## v0.2.0 (2026-05-31)

### Feat

- generate multiallelic records in documents() strategy
- percent-encode reserved chars in string values
- add field-value generation strategy for Number x Type matrix
- public API exports and README
- add Hypothesis strategies and coverage tables
- add reference-aware REF/ALT adapter
- add text + bgzip/index IO
- add VcfBuilder with eager validation
- add ground-truth deriver
- add VCF serializer
- add variant-class constructors and classify
- add Record, ContigDef, VcfDocument model
- add Genotype parse/render
- add curated reserved-field registry
- add FieldDef with validity invariants
- add Number=G genotype ordering
- add Number model and cardinality
- add VCF Type enum

### Fix

- serialize non-finite floats as missing; clearer error for unknown reserved field
