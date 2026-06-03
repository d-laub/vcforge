# Symbolic & non-sequence ALT alleles — design

- **Date:** 2026-06-02
- **Status:** Approved (pass 1 of a multi-pass program toward full VCF 4.5 SV coverage)
- **Spec reference:** `docs/reference/VCFv4.5.tex` (vendored VCF 4.5 LaTeX source)

## Context & motivation

vcfixture currently models only literal-sequence ALT alleles (SNP/MNP/INS/DEL/DELINS)
plus the bare `*` spanning deletion. It has no representation for **symbolic alleles**
(`<DEL>`, `<INS>`, …), the **unspecified allele** `<*>`, or **breakends** (`t[p[`).

Two concrete downstream failures motivate this:

- `gvl.write()` does not filter symbolic alleles and produces garbage haplotypes by
  splicing the raw symbolic-allele string into the sequence.
- genoray PR#51 passes raw `<*>` allele strings into vcfixture, which silently
  mishandles them today — tests built on it are not trustworthy.

The unifying observation: every ALT is either **literal DNA a tool may splice into a
haplotype**, or **something it must not** (`*`, `<…>`, `<*>`, a breakend). The consumers
do not realize symbolic alleles — they **detect and filter** them. So vcfixture's job
here is to (a) *emit* spec-correct VCFs containing these alleles, and (b) give ground
truth that **unambiguously flags, per allele, whether it is literal sequence**, plus the
SV geometry where defined. That per-allele flag is the real consumer contract.

## Goals

1. Represent and serialize symbolic SV alleles `<DEL> <INS> <DUP> <INV> <CNV>` (and
   colon subtypes), the unspecified allele `<*>`, and breakend replacement strings.
2. Model **SV geometry** for the angle-bracket SVs: per-allele `SVLEN` (Number=A),
   reference-consistent span, and derived inclusive end.
3. Emit ground truth that flags each ALT as sequence-or-not and carries resolved SV
   geometry.
4. Make invalid construction **impossible where the type system allows** and a
   **precise eager error** everywhere else — `build()` returns a spec-valid document or
   raises; it never silently emits invalid VCF.
5. Leave clean architectural seams for the deferred passes (below).

## Non-goals (deferred to later passes)

These are explicitly out of scope for pass 1. The architecture must not foreclose them.

- **Full breakend topology**: `MATEID`/`PARID` mate linking, `EVENT`/`EVENTTYPE`
  grouping, single breakends, explicit partners, telomeric virtual breakends,
  multi-mate adjacencies, orientation consistency across records. Pass 1 *recognizes,
  flags, and serializes* breakends but does not validate mate topology or generate
  consistent pairs.
- **CNV:TR tandem-repeat family**: `RN`/`RUS`/`RUL`/`RUC`/`RB`/`CIRUC`/`CIRB`/`RUB`
  list-of-list encoding. `<CNV:TR>` is *representable and flagged* but the repeat-unit
  machinery is not modeled.
- **SV FORMAT keys**: `CICN`/`CNQ`/`CNL`/`CNP`/`NQ`/`HAP`/`AHAP`, and the copy-number
  *semantics* of `CN`. `CN` is declarable + buildable in pass 1 and its
  `CN ⇒ equal-SVLEN` cross-allele guard is enforced, but copy-number interpretation,
  the rest of these keys, and fuzzer generation are deferred.
- **`FORMAT/LEN` + `<*>` gVCF reference blocks**: the block-of-reference-calls idiom.
- **IUPAC ambiguity symbolic alleles** (`<R>`, `<M>`, …): valid 4.5 symbolic alleles,
  not generated.
- Localized-allele `LA`/`LR`/`LG` (already out of project scope).

## Approach: a typed `Allele` union (hybrid)

ALTs become a **sealed, typed `Allele` union**, used both as the *public construction
vocabulary* at the builder boundary and as the *stored* representation on `Record`.
This mirrors the existing model precedent: `Record` already stores typed `Genotype`
objects in its sample dicts, not GT strings.

Why typed objects and not strings: the builder accepts only `Allele` (the maximally
type-safe option), so we always hold parsed objects at construction. Storing them
avoids a parse→render→re-parse round trip and gives every consumer
(`serialize`/`truth`/`build`) an **exhaustive `match`** over the union — pyrefly enforces
that adding a future member (e.g. a richer breakend, or `CnvTrAllele`) surfaces as a
type error at every site that must handle it. That exhaustiveness is the central seam
that makes the deferred passes safe.

### The union (`allele.py`, new top-level public module, sibling of `genotype.py`)

```
Allele =
    | SequenceAllele(bases: str)          # literal [ACGTN]+; sub-classifies SNP/MNP/INS/DEL/DELINS vs REF
    | SpanningDeletion                    # bare "*"
    | SymbolicAllele(first_type, subtypes)# <DEL>, <CNV:TR>, ...
    | UnspecifiedAllele                   # <*>
    | BreakendAllele(raw, [parsed fields])# t[p[ etc.
```

- **Smart constructors validate at construction ("parse, don't validate")**:
  `SequenceAllele` validates `[ACGTN]+`; `SymbolicAllele("DEL")` rejects an unknown
  *first type* (`{DEL,INS,DUP,INV,CNV}`) while preserving unknown *subtypes* verbatim
  (the spec permits tool-defined subtypes); `BreakendAllele.parse(s)` rejects strings
  that do not match the breakend grammar. Malformed alleles cannot be built, so they
  never reach `record()`.
- `BreakendAllele` **always stores `raw`** and round-trips verbatim; it parses the
  fields it can today (orientation, mate `chr:pos`, inserted sequence, single-vs-paired)
  and leaves room to grow the rest in the BND pass — no model change required then.
- Each member has `.render() -> str` (the wire form). `serialize` calls it; the wire
  format lives cohesively on the alleles.
- Ergonomic aliases keep fixtures terse: `Seq("ACG")`, `Sym.deletion()` /
  `Sym("DUP","TANDEM")`, `Bnd.parse(...)`, `Star()`, `Unspecified()`.
- Re-exported from `__init__.py` (public API).

`variants.py` keeps its sequence helpers; its SNP/MNP/… classification moves behind
`SequenceAllele`, and `record_class` defers to the union for SV/BND/unspecified classes.

## Validity strategy (the construction contract)

VCF validity rules split by what the type system can enforce:

- **Statically enforced** (illegal states unrepresentable): value types; closed sets via
  `Literal`/enums/typed constructors (SV first-type, `SVCLAIM ∈ {D,J,DJ}`, phasing);
  the `Allele` union itself (a typo'd `<DLE>` cannot be expressed).
- **Eager runtime errors** (value-dependent / stateful / grammar — *not* expressible as
  types in Python): per-allele `SVLEN` count (Number=A), GT index `< n_alt`, cardinality,
  REF `[ACGTN]+`, breakend grammar, reference-consistent span. Enforced in `build`/at
  construction with precise messages, never deferred to `serialize`.

> **Contract:** `build()` either returns a spec-valid `VcfDocument` or raises. The
> round-trip oracle test (serialize → pysam/cyvcf2 → assert against `GroundTruth`) is the
> backstop proving the *generated* space really is valid.

## Components

### `_spec/reserved.py` — SV field registry (additive)

Add reserved INFO: `SVLEN` (A, Integer), `SVTYPE` (1, String — legacy), `SVCLAIM`
(A, String), `END` (1, Integer — deprecated/back-compat), `CIPOS`/`CIEND`/`CILEN`
(`.`, Integer), `MATEID`/`PARID` (A, String), `IMPRECISE` (Flag). These are declarations
only; values flow through the existing INFO machinery. (Deferred-pass fields are added
in their respective passes — the registry is purely additive.)

### `model.py`

- `Record.alts` becomes `tuple[Allele, ...]`. `__repr__`/`n_alt` updated to render
  alleles. (Frozen `Allele` dataclasses are hashable — `Record` stays frozen/immutable.)
- New `AltDef(id, description)` frozen dataclass; `VcfDocument.alt_defs: tuple[AltDef, …]`
  carries explicit `##ALT` description overrides.

### `build.py`

- `record(alt: Sequence[Allele], …)` — typed input only.
- An exhaustive `match` over each allele applies the value-dependent rules: padding base
  required for `Symbolic`/`Breakend` (POS = preceding base) and *not* for `Unspecified`
  (`<*>` interval includes POS — the documented exception); per-allele `SVLEN` present &
  defined for INS/DUP/INV/DEL/CNV and missing for breakends/`<*>` (negative normalized to
  absolute); `SVCLAIM` per-type rules (DEL/DUP require D/J/DJ; CNV D-or-missing; INS/INV/BND
  J-or-missing); the `CN ⇒ equal-SVLEN` guard across `<CNV>/<DEL>/<DUP>` alleles;
  reference-consistent span in reference-aware mode; existing GT-index and cardinality
  checks (GT may index symbolic alleles).
- New `.alt(id, description)` to override a `##ALT` description.
- **Auto-`##ALT` headers**: the builder collects every symbolic-allele ID actually used
  and emits the corresponding `##ALT` lines (default description, or the `.alt()` override).
  This deletes the "symbolic ALT used but undeclared" error class structurally.

### `serialize.py`

- `_render_record` renders each allele via `.render()`. Symbolic/breakend strings are
  **not** percent-encoded (they are not INFO values).
- Emit auto-collected `##ALT=<ID=…,Description="…">` meta lines (unioned with `.alt()`
  overrides). INFO rendering (`SVLEN`/`SVCLAIM`/`END`) unchanged.

### `truth.py` — the consumer contract

New frozen dataclass, derived by matching each `Allele` + resolving INFO:

```
AlleleTruth(
    kind: str,            # SNP|MNP|INS|DEL|DELINS|SPANNING_DEL|SYMBOLIC|UNSPECIFIED|BND
    is_sequence: bool,    # True iff literal DNA a tool may splice — the headline flag
    sv_type: str | None,  # "DEL"/"CNV"/"DUP:TANDEM"... for symbolic; else None
    svlen: int | None,    # resolved per-allele (abs); None where undefined
    sv_end: int | None,   # 1-based inclusive end = POS + svlen for DEL/DUP/INV/CNV; None otherwise
)
```

`GroundTruth` gains:
- `alts_truth: list[list[AlleleTruth]]` — per record, per ALT. ALT index `i` maps to
  genotype allele index `i+1`, so consumers cross-reference the existing `genotypes`
  matrix against `is_sequence` to know which *called* alleles to filter.
- `is_sequence_mask: list[np.ndarray]` — ragged per-record convenience for the common
  "drop any non-sequence allele" assertion.

`variant_class` vocabulary extends to `SV_DEL`/`SV_INS`/`SV_DUP`/`SV_INV`/`CNV`/`BND`/
`UNSPECIFIED`, with `MULTIALLELIC` still covering mixed records (`A,<DEL>,<*>`). The
genotype matrix, phasing, and INFO/FORMAT truth are unchanged — symbolic alleles remain
indexed alleles; this only adds the per-allele classification layer.

The `<*>` case lands cleanly: `kind=UNSPECIFIED`, `is_sequence=False`, `svlen=None`, POS
included in the reference interval — a consumer test asserts "index k unspecified → filter"
without string special-casing.

### `strategies.py`

Feed the *same* builder (house rule):
- Add `SV_DEL/SV_INS/SV_DUP/SV_INV/CNV` to the variant-class vocabulary. Reference-aware
  generation draws POS + single preceding REF base + a fitting `SVLEN` so `POS+SVLEN`
  stays within the contig, emitting matching `SVLEN`/`SVCLAIM`. Reference-free mode draws
  an arbitrary REF base with consistent `SVLEN`.
- `<*>` generation as an extra ALT (realistic gVCF-ish shape) and standalone.
- **Breakends are not fuzzed** (recognize+flag only) — covered by a small fixed table of
  explicit example records (paired bracket forms, single breakend, inserted-sequence form)
  driven through the builder.

## Extensibility seams (door left open)

- **Sealed union + exhaustive `match`**: adding a future allele member (richer breakend,
  `CnvTrAllele`, IUPAC) is a type error at every handling site until addressed.
- **`BreakendAllele.raw`**: full-topology validation (MATEID/PARID/EVENT, partners,
  telomeres, single breakends) is a future *document-level validator* over already-stored
  raw strings — no model change.
- **Additive reserved registry**: deferred-pass fields (CNV:TR family, SV FORMAT keys,
  LEN) slot in without touching existing entries.
- **`AlleleTruth` is a dataclass**: future fields (copy number, repeat units) extend it
  without breaking the existing per-allele contract.
- **`AltDef`/auto-`##ALT`**: already generalizes to IUPAC and any future symbolic ID.

## Testing

- Extend the oracle self-validation: serialize → parse via pysam/cyvcf2 → assert the
  third-party decode of ALT, `SVLEN`/`END`, and per-allele kind matches `AlleleTruth`
  (`is_sequence`, `svlen`, `sv_end`).
- Hypothesis property test asserting the round-trip for arbitrary symbolic-SV documents.
- Explicit breakend fixture table (build → serialize → parse → assert flagged BND,
  `is_sequence=False`, verbatim round-trip).
- `test_benchmark.py` budget extended to the new draw paths.
- `src/` stays `pyrefly check` strict-clean; the union's exhaustiveness is part of that.

## Migration impact

`record(alt=…)` now requires `Allele` objects, so existing fixtures/tests wrap literal
bases in `Seq(...)` (and `*` becomes `Star()`). Tests are out of the type-check scope but
the call-site change is mechanical and broad; it is part of this pass.
