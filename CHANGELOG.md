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
