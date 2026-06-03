---
name: vcfixture
description: Use when writing or modifying tests that import `vcfixture` to generate VCF fixtures with a decoded ground-truth oracle — building documents with VcfBuilder, deriving GroundTruth numpy arrays, or fuzzing a parser with the Hypothesis strategies. For maintainers of consuming libraries (genoray, GenVarLoader) replacing hand-coded expected arrays.
---

# vcfixture

## Overview

vcfixture generates small spec-conformant VCFs (v4.1–4.5) **and returns the
decoded ground truth alongside**, so parser tests assert against a known oracle
instead of hand-maintained numpy literals. One immutable `VcfDocument` is the
hub: a builder (or Hypothesis strategies) produces it; `.render()`/`.write()`
serialize it; `.truth()` derives the oracle. The oracle and the bytes can never
diverge — they come from the same object.

## When to use

- Replacing hand-coded VCF fixtures + hand-derived expected arrays in a parser test.
- Asserting a parser's genotype matrix / INFO / FORMAT decode against a trusted oracle.
- Property-testing a parser across the Number×Type matrix or variant classes (use the strategies).

## The five things that bite you first

These are the API shapes that are NOT what you'd guess. Get them right and the rest follows.

| You'd guess | Actually |
|-------------|----------|
| `VcfBuilder(samples=...)`, then `.add_contig()` | `VcfBuilder(samples=[...], contigs=[("chr1", 248956422)])` — contigs are `(id, length)` pairs at construction; `length` may be `None`. There is no `add_contig`. |
| `.add_info(number=, type=)` | `.info("AF", Number.A, Type.FLOAT)` / `.fmt("DS", Number.A, Type.FLOAT)` — chained, return `self`. **Reserved IDs (AF, AC, AN, DP, AD, PL, GL, GT, SVLEN, SVCLAIM…) resolve Number/Type automatically when omitted: just `.info("AF")`.** Prefer omitting them — an explicit Number/Type is taken verbatim, NOT checked against the reserved definition, so a wrong pair is silently accepted. |
| `alt=["T"]` (string) | `alt=[Seq("T")]` — **ALT must be typed `Allele` instances**, never strings. `Seq`/`Sym`/`Star`/`Unspecified`/`Bnd`. |
| `Genotype(alleles=[0,1], phased=True)` | `gt=["0\|1", "1\|1"]` — genotypes are **strings** per sample (`\|`=phased, `/`=unphased, `.`=missing). |
| FORMAT via `formats={"DS": ...}` | FORMAT values are **keyword args**: `.record(..., DS=[[0.4], [1.9]])`. Shape is `[per-sample][per-value]` (inner list is per-ALT for Number=A/R). |

Outputs are **builder methods**, not module functions: `doc.render()` (text),
`doc.truth()` (`GroundTruth`), `doc.write(path, bgzip=True, index=True)`. There
is no `vcfixture.io.write` / `vcfixture.truth.derive` in the public API.

## Canonical example (end-to-end parser test)

```python
import numpy as np
from cyvcf2 import VCF  # your parser under test
from vcfixture import VcfBuilder, Number, Type, Seq

def test_parser_matches_oracle(tmp_path):
    doc = (
        VcfBuilder(samples=["s1", "s2"], contigs=[("chr1", 100_000)])
        .info("AF")                          # AF is reserved → Number/Type resolved for you
        .fmt("GT")
        .fmt("DS", Number.A, Type.FLOAT)
        .record(
            "chr1", 1000,
            ref="A", alt=[Seq("T")],         # typed Allele, not "T"
            gt=["0|1", "1|1"],               # one string per sample
            info={"AF": [0.25]},             # list even for Number=A length-1
            DS=[[0.4], [1.9]],               # FORMAT as kwarg, [sample][per-ALT]
        )
    )
    path = doc.write(tmp_path / "fix.vcf.gz", bgzip=True, index=True)
    t = doc.truth()

    vcf = VCF(str(path))
    rec = next(iter(vcf))

    # GroundTruth is indexed by RECORD first.
    np.testing.assert_array_equal(t.genotypes[0], [[0, 1], [1, 1]])  # (rec,sample,ploidy)
    assert t.phasing[0, 0] == True                                   # (rec,sample) bool
    assert t.info[0]["AF"] == [0.25]                                 # list-per-record of dicts
    assert t.format[0][0]["DS"] == [0.4]                             # [rec][sample][id]

    # assert your parser's decode against those oracle values...
    parsed_gt = np.array([g[:2] for g in rec.genotypes], dtype=np.int32)
    np.testing.assert_array_equal(parsed_gt, t.genotypes[0])
```

## GroundTruth shape (the other common mistake)

`GroundTruth` is **indexed by record first**, never a single flat dict:

- `genotypes` — `(records, samples, ploidy)` int ndarray; **`-1` = missing** (matches genoray/GVL).
- `phasing` — `(records, samples)` bool; `True` iff fully phased.
- `pos` — `(records,)` int64, 1-based. `ref` / `alts` — per-record strings.
- `info` — `list[dict]`, one `{id: value}` per record → `t.info[rec][id]`.
- `format` — `list[list[dict]]`, per record per sample → `t.format[rec][sample][id]`.
- `variant_class` — one label per record: `SNP MNP INS DEL DELINS SPANNING_DEL`, `MULTIALLELIC` (any >1 ALT), `SV_DEL`/`SV_DUP`/… (symbolic), `BND`.
- `alts_truth[rec][alt]` — `AlleleTruth` with `.kind`, `.is_sequence`, `.sv_type`, `.svlen`, `.sv_end`.

## Allele types (ALT must be one of these)

`Seq("A")` sequence · `Star()` spanning deletion `*` · `Unspecified()` `<*>` ·
`Bnd.parse("T[chr2:5[")` breakend · `Sym.deletion()`/`.insertion()`/
`.duplication()`/`.inversion()`/`.cnv()` symbolic (`Sym.duplication("TANDEM")`
for subtypes). Symbolic + breakend ALTs need a single-base REF padding base;
symbolic alleles need `SVLEN` in INFO; `<DEL>`/`<DUP>` need `SVCLAIM` (`"D"`/`"J"`/`"DJ"`).

## Spec version

Default is the latest (4.5). Pass `version=VcfVersion.V4_2` to `VcfBuilder` to
pin a version; the builder rejects fields/encodings not valid for it (e.g.
`SVCLAIM` before 4.4, and `SVLEN`'s sign convention flips at 4.3/4.4).

## Property-testing a parser (Hypothesis strategies)

`from vcfixture import strategies`. These feed the *same* builder, so shrinking
works end to end. Draw a doc, round-trip it through your parser, assert against
its truth:

```python
from hypothesis import given
from vcfixture import strategies as vs

@given(vs.documents())                       # small VcfDocument
def test_roundtrip(doc):
    truth = doc.truth()
    # write doc.render() / doc.write(...), parse, assert vs truth

@given(vs.reference_and_documents())          # (ReferenceSpec, doc, GroundTruth) triple
def test_reference_consistent(triple):
    spec, doc, truth = triple
```

Key strategies: `documents()` (plain), `documents_with_fields()` (full
Number×Type matrix), `symbolic_documents()` (SV alleles),
`reference_and_documents()` (reference-consistent triple, for tools that run
`bcftools norm`/`consensus`). Pass `reference=` / `violations=` to `documents()`
for reference-aware mode. Drive variant-class / Number×Type coverage with
`@pytest.mark.parametrize` over `vs.ALL_VARIANT_CLASSES` /
`vs.ALL_NUMBER_TYPE_COMBOS`, reserving Hypothesis for values within a combo.

## Common mistakes

| Symptom | Fix |
|---------|-----|
| `TypeError`/validation passing `alt=["T"]` | Wrap in a typed allele: `alt=[Seq("T")]`. |
| `AttributeError: add_contig` / `add_info` | Contigs go in the constructor; declare fields with chained `.info()`/`.fmt()`. |
| Trying `vcfixture.io.write(...)` / `truth.derive(...)` | Use the builder methods `doc.write(...)` / `doc.truth()`. |
| `truth.info["AF"]` → KeyError/TypeError | Index by record first: `truth.info[0]["AF"]`. |
| Passing FORMAT in a `formats=` dict | FORMAT fields are keyword args on `.record()`. |
| Genotype missing decoded as `.` vs `-1` | Oracle uses `-1` for missing in `genotypes`; match that. |
| `ValueError: ... cardinality` | Value count must equal resolved Number (Number=R is per-allele *including REF*; A is per-ALT; G is per-genotype). |
```
