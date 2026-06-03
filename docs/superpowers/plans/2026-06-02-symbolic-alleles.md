# Symbolic & Non-Sequence ALT Alleles Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Represent, serialize, validate, and flag symbolic SV alleles (`<DEL>/<INS>/<DUP>/<INV>/<CNV>` + subtypes), the unspecified allele `<*>`, and (recognize+flag) breakends, with per-allele ground truth that tells consumers which ALTs are literal sequence.

**Architecture:** A new sealed, typed `Allele` union becomes both the public construction vocabulary and the stored `Record.alts` type (mirroring how `Genotype` is already stored). Geometry (`SVLEN`) lives in INFO; `truth.py` derives a per-allele `AlleleTruth` carrying an `is_sequence` flag + resolved span. Builder validates the value-dependent rules eagerly; the round-trip oracle (serialize → cyvcf2 → assert) backstops the generated space.

**Tech Stack:** Python 3.10–3.13, `uv`, pytest, Hypothesis, numpy, pysam/cyvcf2, pyrefly (strict on `src/`), ruff. Spec reference: `docs/reference/VCFv4.5.tex`.

**Spec:** `docs/superpowers/specs/2026-06-02-symbolic-alleles-design.md`

---

## File Structure

| File | Responsibility | Action |
|---|---|---|
| `src/vcfixture/allele.py` | Typed `Allele` union, smart constructors, `classify_allele`, `.render()` | **Create** |
| `src/vcfixture/_spec/reserved.py` | SV reserved INFO/FORMAT field registry | Modify |
| `src/vcfixture/variants.py` | `record_class` over `Allele`; sequence helpers unchanged | Modify |
| `src/vcfixture/model.py` | `Record.alts: tuple[Allele,...]`; `AltDef`; `VcfDocument.alt_defs` | Modify |
| `src/vcfixture/serialize.py` | Render alleles; emit `##ALT` lines | Modify |
| `src/vcfixture/truth.py` | `AlleleTruth`, `GroundTruth.alts_truth` + `is_sequence_mask` | Modify |
| `src/vcfixture/build.py` | `record(alt: Sequence[Allele])`; eager SV validation; `.alt()` | Modify |
| `src/vcfixture/strategies.py` | Wrap alts as `Allele`; generate SV + `<*>` records | Modify |
| `src/vcfixture/__init__.py` | Re-export `Allele` vocabulary | Modify |
| `tests/test_allele.py` | Unit tests for the union + classifier | **Create** |
| `tests/test_symbolic_*.py` | SV build/serialize/truth/roundtrip + breakend fixtures | **Create** |
| existing `tests/*` constructing `Record`/`record(alt=...)` | Migrate to `Allele` objects | Modify |

**Conventions to follow:** frozen dataclasses subclass `CompactRepr` and define `@override __repr__` (see `genotype.py`). `from __future__ import annotations` at top of every module. Run everything via `uv run`. Commit messages use Conventional Commits. Do NOT touch `pyproject.toml` version or `CHANGELOG.md`.

---

## Task 1: The `Allele` union (`allele.py`)

Pure new, isolated module. Suite stays green.

**Files:**
- Create: `src/vcfixture/allele.py`
- Test: `tests/test_allele.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_allele.py
import pytest

from vcfixture.allele import (
    BreakendAllele,
    SequenceAllele,
    SpanningDeletion,
    SymbolicAllele,
    UnspecifiedAllele,
    classify_allele,
)


def test_sequence_allele_renders_and_validates():
    assert SequenceAllele("ACGT").render() == "ACGT"
    with pytest.raises(ValueError, match="bases"):
        SequenceAllele("AC-GT")


def test_symbolic_allele_type_and_render():
    assert SymbolicAllele("DEL").render() == "<DEL>"
    assert SymbolicAllele("DUP", ("TANDEM",)).render() == "<DUP:TANDEM>"
    assert SymbolicAllele("DUP", ("TANDEM",)).type_str == "DUP:TANDEM"


def test_symbolic_rejects_unknown_first_type_keeps_unknown_subtype():
    with pytest.raises(ValueError, match="first type"):
        SymbolicAllele("DLE")
    # unknown SUBTYPE is preserved, not rejected
    assert SymbolicAllele("DEL", ("FOO",)).render() == "<DEL:FOO>"


def test_spanning_and_unspecified_render():
    assert SpanningDeletion().render() == "*"
    assert UnspecifiedAllele().render() == "<*>"


def test_breakend_parse_paired_and_single():
    assert BreakendAllele.parse("T[chr2:5[").render() == "T[chr2:5["
    assert BreakendAllele.parse("]chr2:5]T").render() == "]chr2:5]T"
    bnd = BreakendAllele.parse(".TGCA")
    assert bnd.single is True and bnd.render() == ".TGCA"
    with pytest.raises(ValueError, match="breakend"):
        BreakendAllele.parse("T[chr2:5]")  # mismatched brackets


def test_classify_allele_dispatch():
    assert isinstance(classify_allele("ACGT"), SequenceAllele)
    assert isinstance(classify_allele("*"), SpanningDeletion)
    assert isinstance(classify_allele("<*>"), UnspecifiedAllele)
    assert classify_allele("<DUP:TANDEM>") == SymbolicAllele("DUP", ("TANDEM",))
    assert isinstance(classify_allele("C[2:321682["), BreakendAllele)
    assert isinstance(classify_allele("G."), BreakendAllele)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_allele.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'vcfixture.allele'`

- [ ] **Step 3: Implement `allele.py`**

```python
# src/vcfixture/allele.py
from __future__ import annotations

import re
from dataclasses import dataclass

from typing_extensions import assert_never

from ._repr import CompactRepr, override

_SEQ_RE = re.compile(r"^[ACGTNacgtn]+$")
_SV_FIRST_TYPES = frozenset({"DEL", "INS", "DUP", "INV", "CNV"})
# t[p[ , t]p] , [p[t , ]p]t  — both brackets identical (\2 backref); mate is chr:pos
_BND_PAIRED_RE = re.compile(r"^([ACGTNacgtn]*)([\[\]])([^\[\]]+:\d+)\2([ACGTNacgtn]*)$")
_BND_SINGLE_RE = re.compile(r"^(?:\.[ACGTNacgtn]+|[ACGTNacgtn]+\.)$")


@dataclass(frozen=True)
class SequenceAllele(CompactRepr):
    bases: str

    def __post_init__(self) -> None:
        if not _SEQ_RE.match(self.bases):
            raise ValueError(f"sequence allele bases must be [ACGTN]+, got {self.bases!r}")

    def render(self) -> str:
        return self.bases

    @override
    def __repr__(self) -> str:
        return f"Seq({self.bases})"


@dataclass(frozen=True)
class SpanningDeletion(CompactRepr):
    def render(self) -> str:
        return "*"

    @override
    def __repr__(self) -> str:
        return "Star()"


@dataclass(frozen=True)
class SymbolicAllele(CompactRepr):
    first_type: str
    subtypes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.first_type not in _SV_FIRST_TYPES:
            raise ValueError(
                f"symbolic SV first type must be one of {sorted(_SV_FIRST_TYPES)}, "
                f"got {self.first_type!r}"
            )

    @property
    def type_str(self) -> str:
        return ":".join((self.first_type, *self.subtypes))

    def render(self) -> str:
        return f"<{self.type_str}>"

    @override
    def __repr__(self) -> str:
        return f"Sym(<{self.type_str}>)"

    @classmethod
    def deletion(cls, *subtypes: str) -> SymbolicAllele:
        return cls("DEL", subtypes)

    @classmethod
    def insertion(cls, *subtypes: str) -> SymbolicAllele:
        return cls("INS", subtypes)

    @classmethod
    def duplication(cls, *subtypes: str) -> SymbolicAllele:
        return cls("DUP", subtypes)

    @classmethod
    def inversion(cls, *subtypes: str) -> SymbolicAllele:
        return cls("INV", subtypes)

    @classmethod
    def cnv(cls, *subtypes: str) -> SymbolicAllele:
        return cls("CNV", subtypes)


@dataclass(frozen=True)
class UnspecifiedAllele(CompactRepr):
    def render(self) -> str:
        return "<*>"

    @override
    def __repr__(self) -> str:
        return "Unspecified()"


@dataclass(frozen=True)
class BreakendAllele(CompactRepr):
    raw: str
    single: bool = False

    def render(self) -> str:
        return self.raw

    @override
    def __repr__(self) -> str:
        return f"Bnd({self.raw})"

    @classmethod
    def parse(cls, s: str) -> BreakendAllele:
        if _BND_SINGLE_RE.match(s):
            return cls(raw=s, single=True)
        if _BND_PAIRED_RE.match(s):
            return cls(raw=s, single=False)
        raise ValueError(f"not a valid breakend replacement string: {s!r}")


Allele = (
    SequenceAllele
    | SpanningDeletion
    | SymbolicAllele
    | UnspecifiedAllele
    | BreakendAllele
)

# Ergonomic aliases for terse fixtures.
Seq = SequenceAllele
Sym = SymbolicAllele
Star = SpanningDeletion
Unspecified = UnspecifiedAllele
Bnd = BreakendAllele


def classify_allele(alt: str) -> Allele:
    """Parse a raw ALT string into a typed Allele (syntactic dispatch)."""
    if alt == "*":
        return SpanningDeletion()
    if alt == "<*>":
        return UnspecifiedAllele()
    if alt.startswith("<") and alt.endswith(">"):
        parts = alt[1:-1].split(":")
        return SymbolicAllele(parts[0], tuple(parts[1:]))
    if "[" in alt or "]" in alt:
        return BreakendAllele.parse(alt)
    if len(alt) > 1 and (alt.startswith(".") or alt.endswith(".")):
        return BreakendAllele.parse(alt)
    return SequenceAllele(alt)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_allele.py -q && uv run pyrefly check`
Expected: PASS; pyrefly clean.

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/allele.py tests/test_allele.py
git commit -m "feat: typed Allele union with classifier and smart constructors"
```

---

## Task 2: SV reserved fields (`_spec/reserved.py`)

Additive. Suite stays green.

**Files:**
- Modify: `src/vcfixture/_spec/reserved.py`
- Test: `tests/test_reserved.py`

- [ ] **Step 1: Write the failing test (append to `tests/test_reserved.py`)**

```python
def test_sv_reserved_info_fields():
    from vcfixture._spec.number import NumberKind
    from vcfixture._spec.reserved import reserved
    from vcfixture._spec.types import Type

    svlen = reserved("SVLEN", "INFO")
    assert svlen.number.kind is NumberKind.A and svlen.type is Type.INTEGER
    assert reserved("SVCLAIM", "INFO").number.kind is NumberKind.A
    assert reserved("END", "INFO").type is Type.INTEGER
    assert reserved("MATEID", "INFO").type is Type.STRING
    assert reserved("IMPRECISE", "INFO").type is Type.FLAG
    assert reserved("CN", "FORMAT").type is Type.FLOAT
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_reserved.py::test_sv_reserved_info_fields -q`
Expected: FAIL with `KeyError: 'SVLEN'`

- [ ] **Step 3: Add the fields to `_spec/reserved.py`**

Add these entries to the `_INFO` dict (after `"H2"`):

```python
    "END": FieldDef("END", Number.ONE, Type.INTEGER, "End position (deprecated)", "INFO"),
    "SVTYPE": FieldDef("SVTYPE", Number.ONE, Type.STRING, "Type of structural variant", "INFO"),
    "SVLEN": FieldDef("SVLEN", Number.A, Type.INTEGER, "Length of structural variant", "INFO"),
    "SVCLAIM": FieldDef("SVCLAIM", Number.A, Type.STRING, "Structural variant claim", "INFO"),
    "CIPOS": FieldDef("CIPOS", Number.DOT, Type.INTEGER, "Confidence interval around POS", "INFO"),
    "CIEND": FieldDef("CIEND", Number.DOT, Type.INTEGER, "Confidence interval around END", "INFO"),
    "CILEN": FieldDef("CILEN", Number.DOT, Type.INTEGER, "Confidence interval around SVLEN", "INFO"),
    "MATEID": FieldDef("MATEID", Number.A, Type.STRING, "ID of mate breakend", "INFO"),
    "PARID": FieldDef("PARID", Number.A, Type.STRING, "ID of partner breakend", "INFO"),
    "IMPRECISE": FieldDef("IMPRECISE", Number.FLAG, Type.FLAG, "Imprecise structural variant", "INFO"),
```

Add to the `_FORMAT` dict (after `"PS"`):

```python
    "CN": FieldDef("CN", Number.ONE, Type.FLOAT, "Copy number", "FORMAT"),
    "LEN": FieldDef("LEN", Number.ONE, Type.INTEGER, "Length of <*> reference block", "FORMAT"),
```

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest tests/test_reserved.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/_spec/reserved.py tests/test_reserved.py
git commit -m "feat: add reserved structural-variant INFO/FORMAT fields"
```

---

## Task 3: Migrate the model to store `Allele` objects

This is the cross-cutting type change. `Record.alts` becomes `tuple[Allele, ...]`; the builder accepts only `Allele`; serialize/truth/variants render alleles. **No new behavior** — every document renders byte-identically. The task is green only at the end (all call sites + directly-constructing tests migrate together). Work through the steps in order; run the full suite at the end.

**Files:**
- Modify: `src/vcfixture/model.py`, `src/vcfixture/variants.py`, `src/vcfixture/serialize.py`, `src/vcfixture/truth.py`, `src/vcfixture/build.py`, `src/vcfixture/strategies.py`
- Modify (tests): `tests/test_model.py`, `tests/test_truth.py`, `tests/test_serialize.py`, `tests/test_repr.py`, `tests/test_build.py`, `tests/test_io.py`, `tests/test_labels.py`, `tests/test_benchmark.py`, `tests/test_genoray_parity.py`, `tests/test_variants.py`

- [ ] **Step 1: `model.py` — type and repr**

Change the `alts` field annotation (`model.py:40`) and `__repr__` (`model.py:50`):

```python
    alts: tuple[Allele, ...]  # typed alleles; may include SpanningDeletion()/symbolic/...
```

```python
        alts = ",".join(a.render() for a in self.alts) if self.alts else "."
```

Add the import near the other model imports:

```python
from .allele import Allele
```

- [ ] **Step 2: `variants.py` — `record_class` over `Allele`**

Replace `record_class` (`variants.py:43-46`) and add imports:

```python
from typing_extensions import assert_never

from .allele import (
    Allele,
    BreakendAllele,
    SequenceAllele,
    SpanningDeletion,
    SymbolicAllele,
    UnspecifiedAllele,
)


def record_class(ref: str, alts: tuple[Allele, ...]) -> str:
    if len(alts) != 1:
        return "MULTIALLELIC"
    a = alts[0]
    if isinstance(a, SequenceAllele):
        return classify(ref, a.bases)
    if isinstance(a, SpanningDeletion):
        return "SPANNING_DEL"
    if isinstance(a, UnspecifiedAllele):
        return "UNSPECIFIED"
    if isinstance(a, BreakendAllele):
        return "BND"
    if isinstance(a, SymbolicAllele):
        return a.first_type if a.first_type == "CNV" else f"SV_{a.first_type}"
    assert_never(a)
```

(Leave `classify` and the `snp`/`mnp`/… helpers unchanged.)

- [ ] **Step 3: `serialize.py` — render alleles**

Change `serialize.py:78`:

```python
    alt = ",".join(a.render() for a in rec.alts) if rec.alts else "."
```

- [ ] **Step 4: `truth.py` — render alts list**

Change `truth.py:45` (the `alts.append(...)` line) to render strings (preserving `GroundTruth.alts: list[list[str]]`):

```python
        alts.append([a.render() for a in rec.alts])
```

(`record_class(rec.ref, rec.alts)` on the next line now receives the typed tuple — no change needed.)

- [ ] **Step 5: `build.py` — accept `Allele` only**

Add import:

```python
from .allele import Allele
```

Change the `alt` parameter type (`build.py:82`):

```python
        alt: Sequence[Allele],
```

The body already does `alts = tuple(alt)` and stores `alts=alts` — no further change this task. GT-index check uses `n_alt = len(alts)` and is unchanged (symbolic alleles are still indexed).

- [ ] **Step 6: `strategies.py` — wrap produced alts as alleles**

The strategies build raw ALT strings then call `b.record(alt=...)`. Wrap at the two record-building call sites so they pass `Allele` objects.

Add import:

```python
from .allele import Allele, classify_allele
```

Reference-free body (`strategies.py:411`) — change `alt=alts` to:

```python
        b.record("chr1", pos, ref=ref, alt=[classify_allele(a) for a in alts], gt=gts)
```

Reference-aware body (`strategies.py:324-331`, the `b.record(...)` call) — change `alt=alts` to:

```python
            alt=[classify_allele(a) for a in alts],
```

Matrix strategy (`strategies.py:172-183`, `documents_with_fields`) — same change to its `alt=alts`:

```python
            alt=[classify_allele(a) for a in alts],
```

- [ ] **Step 7: Migrate directly-constructing tests**

Apply these exact edits. In each, add `from vcfixture.allele import Seq, Star` at the top and wrap literal ALTs.

`tests/test_model.py`: `alts=("A",)` → `alts=(Seq("A"),)`; in `test_ploidy_varies_uses_max` the positional `("T",)` → `(Seq("T"),)`; assertion `doc.records[0].alts == ("A",)` → `== (Seq("A"),)`.

`tests/test_truth.py`: `alts=("A",)` → `alts=(Seq("A"),)`. (The `t.alts == [["A"]]` assertion stays — truth renders to strings.)

`tests/test_serialize.py`: `alts=("A",)` → `alts=(Seq("A"),)`.

`tests/test_repr.py`: `_make_record(alts=("T", "G"))` default → `alts=(Seq("T"), Seq("G"))`; the `alts=()` case stays `()`.

`tests/test_build.py`: every `alt=["A"]`/`alt=["T"]` → `alt=[Seq("A")]`/`alt=[Seq("T")]`.

`tests/test_io.py`, `tests/test_labels.py`, `tests/test_benchmark.py`: `alt=["T"]`/`alt=["C"]` → `alt=[Seq("T")]`/`alt=[Seq("C")]`.

`tests/test_genoray_parity.py`: each `alt=["A"]`/`alt=["C"]` → `alt=[Seq("A")]`/`alt=[Seq("C")]`.

`tests/test_variants.py` `test_record_class_multiallelic`:

```python
def test_record_class_multiallelic():
    from vcfixture.allele import Seq

    assert v.record_class("G", (Seq("A"), Seq("C"))) == "MULTIALLELIC"
    assert v.record_class("G", (Seq("A"),)) == "SNP"
```

- [ ] **Step 8: Run the full suite + type check**

Run: `uv run pytest -q && uv run pyrefly check`
Expected: ALL PASS; pyrefly strict-clean. (Documents render byte-identically; the round-trip tests prove no behavior change.)

- [ ] **Step 9: Commit**

```bash
git add src/vcfixture tests
git commit -m "refactor: store typed Allele objects in Record.alts"
```

---

## Task 4: Per-allele ground truth (`AlleleTruth`)

Additive to `GroundTruth`. Derives `is_sequence` for every ALT and SV geometry where present.

**Files:**
- Modify: `src/vcfixture/truth.py`
- Test: `tests/test_symbolic_truth.py` (create)

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_symbolic_truth.py
import numpy as np

from vcfixture.allele import Seq, Star, Sym, Unspecified
from vcfixture.genotype import Genotype
from vcfixture.model import ContigDef, Record, VcfDocument


def _doc(alts, info):
    rec = Record(
        chrom="chr1", pos=100, ids=None, ref="G", alts=alts, qual=None,
        filters=None, info=info, fmt_keys=("GT",),
        samples=({"GT": Genotype.parse("0/1")},),
    )
    return VcfDocument("VCFv4.5", (), (), (), (ContigDef("chr1", 1000),), ("s1",), (rec,))


def test_symbolic_del_geometry_and_flag():
    t = _doc((Sym.deletion(),), {"SVLEN": [50]}).truth()
    at = t.alts_truth[0][0]
    assert at.kind == "SYMBOLIC" and at.is_sequence is False
    assert at.sv_type == "DEL" and at.svlen == 50 and at.sv_end == 150
    assert t.variant_class == ["SV_DEL"]


def test_insertion_has_no_end():
    t = _doc((Sym.insertion(),), {"SVLEN": [30]}).truth()
    at = t.alts_truth[0][0]
    assert at.svlen == 30 and at.sv_end is None


def test_negative_svlen_normalized_absolute():
    t = _doc((Sym.deletion(),), {"SVLEN": [-50]}).truth()
    assert t.alts_truth[0][0].svlen == 50


def test_unspecified_and_spanning_and_mixed_mask():
    t = _doc((Seq("A"), Star(), Unspecified()), {}).truth()
    kinds = [a.kind for a in t.alts_truth[0]]
    assert kinds == ["SNP", "SPANNING_DEL", "UNSPECIFIED"]
    np.testing.assert_array_equal(t.is_sequence_mask[0], [True, False, False])
    assert t.variant_class == ["MULTIALLELIC"]
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run pytest tests/test_symbolic_truth.py -q`
Expected: FAIL — `AttributeError: 'GroundTruth' object has no attribute 'alts_truth'`

- [ ] **Step 3: Implement in `truth.py`**

Add imports:

```python
from typing_extensions import assert_never

from .allele import (
    Allele,
    BreakendAllele,
    SequenceAllele,
    SpanningDeletion,
    SymbolicAllele,
    UnspecifiedAllele,
)
from .variants import classify, record_class
```

Add the dataclass (after the imports, before `GroundTruth`):

```python
@dataclass(frozen=True)
class AlleleTruth:
    kind: str  # SNP|MNP|INS|DEL|DELINS|SPANNING_DEL|SYMBOLIC|UNSPECIFIED|BND
    is_sequence: bool  # True iff literal DNA a tool may splice
    sv_type: str | None  # e.g. "DEL"/"DUP:TANDEM" for symbolic; else None
    svlen: int | None  # resolved per-allele (absolute); None where undefined
    sv_end: int | None  # 1-based inclusive end = POS + svlen for DEL/DUP/INV/CNV; else None


_SV_SPANNING = frozenset({"DEL", "DUP", "INV", "CNV"})


def _allele_truth(ref: str, pos: int, allele: Allele, svlen_val: object) -> AlleleTruth:
    if isinstance(allele, SequenceAllele):
        return AlleleTruth(classify(ref, allele.bases), True, None, None, None)
    if isinstance(allele, SpanningDeletion):
        return AlleleTruth("SPANNING_DEL", False, None, None, None)
    if isinstance(allele, UnspecifiedAllele):
        return AlleleTruth("UNSPECIFIED", False, None, None, None)
    if isinstance(allele, BreakendAllele):
        return AlleleTruth("BND", False, None, None, None)
    if isinstance(allele, SymbolicAllele):
        svlen = None if svlen_val is None else abs(int(svlen_val))
        end = (
            pos + svlen
            if svlen is not None and allele.first_type in _SV_SPANNING
            else None
        )
        return AlleleTruth("SYMBOLIC", False, allele.type_str, svlen, end)
    assert_never(allele)
```

Add two fields to `GroundTruth` (after `variant_class`):

```python
    alts_truth: list[list[AlleleTruth]]  # per record, per ALT
    is_sequence_mask: list[np.ndarray]  # per record: bool array over ALTs
```

In `derive_truth`, add accumulators next to the others:

```python
    alts_truth: list[list[AlleleTruth]] = []
    seq_mask: list[np.ndarray] = []
```

Inside the record loop (after `vclass.append(...)`), add:

```python
        svlen_list = rec.info.get("SVLEN")
        per_alt: list[AlleleTruth] = []
        for ai, allele in enumerate(rec.alts):
            sv = (
                svlen_list[ai]
                if isinstance(svlen_list, (list, tuple)) and ai < len(svlen_list)
                else None
            )
            per_alt.append(_allele_truth(rec.ref, rec.pos, allele, sv))
        alts_truth.append(per_alt)
        seq_mask.append(np.array([a.is_sequence for a in per_alt], dtype=bool))
```

Add to the `GroundTruth(...)` return:

```python
        alts_truth=alts_truth,
        is_sequence_mask=seq_mask,
```

(Remove the now-redundant local `from .variants import record_class` if it was a function-local import; it is now module-level.)

- [ ] **Step 4: Run to verify pass + type check**

Run: `uv run pytest tests/test_symbolic_truth.py -q && uv run pyrefly check`
Expected: PASS; pyrefly clean.

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/truth.py tests/test_symbolic_truth.py
git commit -m "feat: per-allele AlleleTruth with is_sequence flag and SV geometry"
```

---

## Task 5: Builder validation for symbolic alleles

Eager, precise errors for the value-dependent rules.

**Files:**
- Modify: `src/vcfixture/build.py`
- Test: `tests/test_symbolic_build.py` (create)

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_symbolic_build.py
import pytest

from vcfixture._spec.number import Number
from vcfixture._spec.types import Type
from vcfixture.allele import Seq, Sym, Unspecified
from vcfixture.build import VcfBuilder


def _b():
    return (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 1000)])
        .fmt("GT")
        .info("SVLEN", Number.A, Type.INTEGER)
        .info("SVCLAIM", Number.A, Type.STRING)
    )


def test_valid_symbolic_del_builds():
    doc = _b().record(
        "chr1", 100, ref="G", alt=[Sym.deletion()],
        gt=["0/1"], info={"SVLEN": [50], "SVCLAIM": ["DJ"]},
    ).build()
    assert doc.records[0].alts[0].render() == "<DEL>"


def test_symbolic_requires_single_base_ref():
    with pytest.raises(ValueError, match="padding base"):
        _b().record("chr1", 100, ref="GA", alt=[Sym.deletion()],
                    gt=["0/1"], info={"SVLEN": [50], "SVCLAIM": ["DJ"]})


def test_symbolic_sv_requires_svlen():
    with pytest.raises(ValueError, match="SVLEN"):
        _b().record("chr1", 100, ref="G", alt=[Sym.deletion()],
                    gt=["0/1"], info={"SVCLAIM": ["DJ"]})


def test_svclaim_del_must_be_djd():
    with pytest.raises(ValueError, match="SVCLAIM"):
        _b().record("chr1", 100, ref="G", alt=[Sym.deletion()],
                    gt=["0/1"], info={"SVLEN": [50], "SVCLAIM": ["X"]})


def test_unspecified_allele_does_not_require_padding_or_svlen():
    doc = _b().record("chr1", 100, ref="G", alt=[Unspecified()], gt=["0/1"]).build()
    assert doc.records[0].alts[0].render() == "<*>"


def test_breakend_svlen_must_be_missing():
    from vcfixture.allele import Bnd

    with pytest.raises(ValueError, match="missing"):
        _b().record("chr1", 100, ref="G", alt=[Bnd.parse("G[chr2:9[")],
                    gt=["0/1"], info={"SVLEN": [10]})
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run pytest tests/test_symbolic_build.py -q`
Expected: FAIL (records build without raising).

- [ ] **Step 3: Implement validation in `build.py`**

Add imports:

```python
from .allele import (
    Allele,
    BreakendAllele,
    SequenceAllele,
    SpanningDeletion,
    SymbolicAllele,
    UnspecifiedAllele,
)
```

Add a module-level constant and a static validator method on `VcfBuilder`:

```python
_SVCLAIM_RULES = {
    "DEL": {"D", "J", "DJ"},
    "DUP": {"D", "J", "DJ"},
    "CNV": {"D"},
    "INS": {"J"},
    "INV": {"J"},
}
_SVLEN_DEFINED = frozenset({"DEL", "INS", "DUP", "INV", "CNV"})
```

```python
    @staticmethod
    def _validate_alleles(
        ref: str, alts: tuple[Allele, ...], info: Mapping[str, object] | None
    ) -> None:
        info = info or {}
        svlen = info.get("SVLEN")
        svclaim = info.get("SVCLAIM")
        needs_padding = any(
            isinstance(a, (SymbolicAllele, BreakendAllele)) for a in alts
        )
        if needs_padding and len(ref) != 1:
            raise ValueError(
                "symbolic/breakend ALT requires a single preceding REF padding base, "
                f"got REF={ref!r}"
            )
        for i, a in enumerate(alts):
            sv = svlen[i] if isinstance(svlen, (list, tuple)) and i < len(svlen) else None
            cl = svclaim[i] if isinstance(svclaim, (list, tuple)) and i < len(svclaim) else None
            if isinstance(a, SymbolicAllele) and a.first_type in _SVLEN_DEFINED:
                if sv is None:
                    raise ValueError(f"SVLEN required for symbolic allele {a.render()}")
                allowed = _SVCLAIM_RULES[a.first_type]
                if cl is not None and cl not in allowed:
                    raise ValueError(
                        f"SVCLAIM {cl!r} invalid for {a.render()}; allowed {sorted(allowed)}"
                    )
                if a.first_type in {"DEL", "DUP"} and cl is None:
                    raise ValueError(f"SVCLAIM required for {a.render()} (D/J/DJ)")
            elif isinstance(a, (BreakendAllele, UnspecifiedAllele, SpanningDeletion)):
                if sv is not None:
                    raise ValueError(f"SVLEN must be missing for {a.render()}")
        # FORMAT CN guard: all <CNV>/<DEL>/<DUP> alleles share one SVLEN.
        cn_types = [
            (i, a) for i, a in enumerate(alts)
            if isinstance(a, SymbolicAllele) and a.first_type in {"CNV", "DEL", "DUP"}
        ]
        # (CN presence is checked by the caller; this guard only applies if CN is set
        #  on any sample — see _has_cn below.)
        return None
```

Call it from `record()`, right after `alts = tuple(alt)` (so it fires before sample assembly). Note: the CN guard needs sample data, so fold the CN check into `record()` after samples are built:

```python
        alts = tuple(alt)
        n_alt = len(alts)
        self._validate_alleles(ref, alts, info)
```

After the FORMAT fields are assembled (after the `for key, per_sample in fmt_fields.items()` loop), add the CN equal-SVLEN guard:

```python
        if "CN" in fmt_keys:
            cn_svlens = {
                (info or {}).get("SVLEN", [None])[i]
                for i, a in enumerate(alts)
                if isinstance(a, SymbolicAllele) and a.first_type in {"CNV", "DEL", "DUP"}
            }
            if len(cn_svlens) > 1:
                raise ValueError("FORMAT CN requires equal SVLEN across <CNV>/<DEL>/<DUP> alleles")
```

(`Mapping` is already imported in `build.py`.) Negative-SVLEN normalization is a *truth* concern (Task 4) — the builder stores values as given; the spec says readers take the absolute value.

- [ ] **Step 4: Run to verify pass + type check**

Run: `uv run pytest tests/test_symbolic_build.py -q && uv run pyrefly check`
Expected: PASS; pyrefly clean.

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/build.py tests/test_symbolic_build.py
git commit -m "feat: eager builder validation for symbolic/breakend alleles"
```

---

## Task 6: `##ALT` header emission (`AltDef` + serialize)

**Files:**
- Modify: `src/vcfixture/model.py`, `src/vcfixture/build.py`, `src/vcfixture/serialize.py`
- Test: `tests/test_symbolic_serialize.py` (create)

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_symbolic_serialize.py
from vcfixture._spec.number import Number
from vcfixture._spec.types import Type
from vcfixture.allele import Bnd, Sym
from vcfixture.build import VcfBuilder


def _doc():
    return (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 1000)])
        .fmt("GT")
        .info("SVLEN", Number.A, Type.INTEGER)
        .info("SVCLAIM", Number.A, Type.STRING)
        .alt("DEL", "Custom deletion description")
        .record("chr1", 100, ref="G", alt=[Sym.deletion()],
                gt=["0/1"], info={"SVLEN": [50], "SVCLAIM": ["DJ"]})
        .record("chr1", 200, ref="G", alt=[Sym.duplication("TANDEM")],
                gt=["0/1"], info={"SVLEN": [20], "SVCLAIM": ["DJ"]})
        .build()
    )


def test_alt_header_lines_emitted_and_deduped():
    text = _doc().render()
    lines = text.splitlines()
    assert '##ALT=<ID=DEL,Description="Custom deletion description">' in lines
    assert '##ALT=<ID=DUP:TANDEM,Description="DUP:TANDEM structural variant">' in lines
    assert sum(line.startswith("##ALT=<ID=DEL,") for line in lines) == 1


def test_symbolic_and_breakend_alts_render_unencoded():
    text = _doc().render()
    data = [l for l in text.splitlines() if l.startswith("chr1\t100")][0]
    assert data.split("\t")[4] == "<DEL>"
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run pytest tests/test_symbolic_serialize.py -q`
Expected: FAIL — `VcfBuilder` has no `.alt`.

- [ ] **Step 3: Add `AltDef` + `alt_defs` (model.py)**

Add the dataclass (after `ContigDef`):

```python
@dataclass(frozen=True)
class AltDef(CompactRepr):
    id: str
    description: str

    @override
    def __repr__(self) -> str:
        return f"AltDef({self.id})"

    def header_line(self) -> str:
        return f'##ALT=<ID={self.id},Description="{self.description}">'
```

Add a field to `VcfDocument` with a default so existing positional constructors keep working (place it **last**):

```python
    alt_defs: tuple[AltDef, ...] = ()
```

- [ ] **Step 4: Builder `.alt()` + wiring (build.py)**

Add import: `from .model import AltDef, ContigDef, Record, VcfDocument`.
Init: add `self._alt_defs: dict[str, str] = {}` in `__init__`.
Add method:

```python
    def alt(self, id: str, description: str) -> VcfBuilder:
        self._alt_defs[id] = description
        return self
```

In `build()`, collect every symbolic-allele type used across records, union with explicit overrides, and pass `alt_defs`:

```python
        alt_ids: dict[str, str] = {}
        for rec in self._records:
            for a in rec.alts:
                if isinstance(a, SymbolicAllele):
                    alt_ids.setdefault(a.type_str, f"{a.type_str} structural variant")
        alt_ids.update(self._alt_defs)
        alt_defs = tuple(AltDef(i, d) for i, d in alt_ids.items())
```

Add `alt_defs=alt_defs` to the `VcfDocument(...)` return.

- [ ] **Step 5: Emit `##ALT` lines (serialize.py)**

In `render_document`, after the FILTER loop and before/after FORMAT (order is not significant to parsers), add:

```python
    for ad in doc.alt_defs:
        lines.append(ad.header_line())
```

- [ ] **Step 6: Run to verify pass + type check + full suite**

Run: `uv run pytest tests/test_symbolic_serialize.py -q && uv run pytest -q && uv run pyrefly check`
Expected: ALL PASS; pyrefly clean.

- [ ] **Step 7: Commit**

```bash
git add src/vcfixture/model.py src/vcfixture/build.py src/vcfixture/serialize.py tests/test_symbolic_serialize.py
git commit -m "feat: emit ##ALT header lines for symbolic alleles"
```

---

## Task 7: Strategies generate symbolic SV records + round-trip oracle

**Files:**
- Modify: `src/vcfixture/strategies.py`
- Test: `tests/test_symbolic_roundtrip.py` (create)

- [ ] **Step 1: Write the failing test (Hypothesis + cyvcf2 oracle)**

```python
# tests/test_symbolic_roundtrip.py
import tempfile
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings

from vcfixture import strategies as S

cyvcf2 = pytest.importorskip("cyvcf2")


@settings(max_examples=50, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(S.symbolic_documents())
def test_symbolic_alts_round_trip_through_cyvcf2(doc):
    truth = doc.truth()
    d = tempfile.mkdtemp()
    path = doc.write(Path(d) / "x.vcf.gz", bgzip=True, index=True)
    vf = cyvcf2.VCF(str(path))
    for ri, variant in enumerate(vf):
        for ai, alt in enumerate(variant.ALT):
            at = truth.alts_truth[ri][ai]
            assert alt == doc.records[ri].alts[ai].render()
            if at.kind == "SYMBOLIC":
                assert at.is_sequence is False
                got = variant.INFO.get("SVLEN")
                got_i = got[ai] if isinstance(got, (list, tuple)) else got
                assert abs(int(got_i)) == at.svlen
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run pytest tests/test_symbolic_roundtrip.py -q`
Expected: FAIL — `strategies` has no `symbolic_documents`.

- [ ] **Step 3: Add `symbolic_documents` strategy (strategies.py)**

Add import: `from .allele import SymbolicAllele, UnspecifiedAllele`.

```python
_SYMBOLIC_SV = ["DEL", "INS", "DUP", "INV", "CNV"]


@st.composite
def symbolic_documents(draw: DrawFn, max_samples: int = 3, max_records: int = 4):
    """Documents whose records carry symbolic SV alleles and <*>, with consistent
    per-allele SVLEN/SVCLAIM. Reference-free (arbitrary single REF base)."""
    n_samples = draw(st.integers(1, max_samples))
    samples = [f"s{i}" for i in range(n_samples)]
    ploidy = draw(st.integers(1, 2))
    b = (
        VcfBuilder(samples=samples, contigs=[("chr1", 1_000_000)])
        .fmt("GT")
        .info("SVLEN", Number.A, Type.INTEGER)
        .info("SVCLAIM", Number.A, Type.STRING)
    )
    n_rec = draw(st.integers(1, max_records))
    pos = 1000
    for _ in range(n_rec):
        ref = draw(st.sampled_from(_BASES))
        kind = draw(st.sampled_from([*_SYMBOLIC_SV, "UNSPEC"]))
        if kind == "UNSPEC":
            alts = [UnspecifiedAllele()]
            info: dict[str, object] = {}
        else:
            svlen = draw(st.integers(1, 1000))
            claim = "DJ" if kind in {"DEL", "DUP"} else ("D" if kind == "CNV" else "J")
            alts = [SymbolicAllele(kind)]
            info = {"SVLEN": [svlen], "SVCLAIM": [claim]}
        gts = [draw(genotypes(ploidy, n_alt=1)) for _ in samples]
        b.record("chr1", pos, ref=ref, alt=alts, gt=gts, info=info)
        pos += draw(st.integers(2000, 5000))
    return b.build()
```

- [ ] **Step 4: Run to verify pass**

Run: `uv run pytest tests/test_symbolic_roundtrip.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/vcfixture/strategies.py tests/test_symbolic_roundtrip.py
git commit -m "feat: Hypothesis strategy for symbolic SV documents + cyvcf2 oracle"
```

---

## Task 8: Breakend recognize-and-flag fixture tests

Breakends are not fuzzed; cover them with an explicit table driven through the builder.

**Files:**
- Test: `tests/test_breakends.py` (create)

- [ ] **Step 1: Write the test**

```python
# tests/test_breakends.py
import tempfile
from pathlib import Path

import pytest

from vcfixture._spec.number import Number
from vcfixture._spec.types import Type
from vcfixture.allele import Bnd
from vcfixture.build import VcfBuilder

cyvcf2 = pytest.importorskip("cyvcf2")

CASES = ["T[chr2:5[", "]chr2:5]T", "[chr2:5[T", "T]chr2:5]", ".TGCA", "TGCA."]


@pytest.mark.parametrize("bnd", CASES)
def test_breakend_flagged_and_round_trips(bnd):
    doc = (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 1000), ("chr2", 1000)])
        .fmt("GT")
        .info("MATEID", Number.A, Type.STRING)
        .record("chr1", 100, ref="T", alt=[Bnd.parse(bnd)], gt=["0/1"])
        .build()
    )
    t = doc.truth()
    assert t.alts_truth[0][0].kind == "BND"
    assert t.alts_truth[0][0].is_sequence is False
    assert t.variant_class == ["BND"]

    d = tempfile.mkdtemp()
    path = doc.write(Path(d) / "b.vcf.gz", bgzip=True, index=True)
    vf = cyvcf2.VCF(str(path))
    variant = next(iter(vf))
    assert variant.ALT[0] == bnd
```

- [ ] **Step 2: Run to verify pass (no new src code expected)**

Run: `uv run pytest tests/test_breakends.py -q`
Expected: PASS (everything needed exists from Tasks 1/3/4). If a case fails to round-trip through cyvcf2, note it — pysam may require a `##ALT`/`SVTYPE`; adjust the fixture to add `info={"MATEID": ["m1"]}` and a `SVTYPE` declaration rather than weakening the assertion.

- [ ] **Step 3: Commit**

```bash
git add tests/test_breakends.py
git commit -m "test: breakend alleles are flagged non-sequence and round-trip verbatim"
```

---

## Task 9: Public API, benchmark, and final checks

**Files:**
- Modify: `src/vcfixture/__init__.py`, `tests/test_public_api.py`, `tests/test_benchmark.py`

- [ ] **Step 1: Write the failing API test (append to `tests/test_public_api.py`)**

```python
def test_allele_vocabulary_exported():
    import vcfixture

    for name in ["Seq", "Sym", "Star", "Unspecified", "Bnd", "Allele"]:
        assert hasattr(vcfixture, name), name
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run pytest tests/test_public_api.py -q`
Expected: FAIL

- [ ] **Step 3: Re-export from `__init__.py`**

Add import:

```python
from .allele import (
    Allele,
    Bnd,
    BreakendAllele,
    Seq,
    SequenceAllele,
    SpanningDeletion,
    Star,
    Sym,
    SymbolicAllele,
    Unspecified,
    UnspecifiedAllele,
)
```

Add these names to `__all__`:

```python
    "Allele",
    "Seq",
    "Sym",
    "Star",
    "Unspecified",
    "Bnd",
    "SequenceAllele",
    "SymbolicAllele",
    "SpanningDeletion",
    "UnspecifiedAllele",
    "BreakendAllele",
```

- [ ] **Step 4: Extend the benchmark to the new draw path (tests/test_benchmark.py)**

Add a benchmark case that draws + serializes + derives truth for `symbolic_documents()`, mirroring the existing budgeted pattern in the file:

```python
def test_symbolic_draw_serialize_truth_within_budget(benchmark):
    from vcfixture import strategies as S

    @benchmark
    def _run():
        doc = S.symbolic_documents().example()
        doc.render()
        doc.truth()
```

(If `test_benchmark.py` does not use `pytest-benchmark` but a manual timer, follow that file's existing structure instead — replicate its budget assertion against the symbolic path.)

- [ ] **Step 5: Full verification**

Run:
```bash
uv run pytest -q
uv run ruff check && uv run ruff format --check
uv run pyrefly check
```
Expected: ALL PASS; lint clean; format clean; pyrefly strict-clean on `src/`.

- [ ] **Step 6: Clean-room import sanity (mirrors CI wheel smoke test)**

Run: `uv run python -c "import vcfixture; print(vcfixture.Sym.deletion().render())"`
Expected: prints `<DEL>` (confirms no `src/` import reaches a dev-only dependency).

- [ ] **Step 7: Commit**

```bash
git add src/vcfixture/__init__.py tests/test_public_api.py tests/test_benchmark.py
git commit -m "feat: export Allele construction vocabulary; benchmark symbolic path"
```

---

## Self-Review

**Spec coverage:**
- Symbolic SV alleles `<DEL/INS/DUP/INV/CNV>` + subtypes → Tasks 1, 5, 7. ✓
- `<*>` special semantics (no padding, no SVLEN, POS-inclusive) → Tasks 1, 4 (`UNSPECIFIED`), 5 (`test_unspecified...`). ✓
- Breakends recognize+flag, not fuzzed → Tasks 1 (`BreakendAllele.parse`), 8 (fixture table). ✓
- SV geometry (per-allele SVLEN, span, abs value) → Task 4 (`AlleleTruth.svlen/sv_end`, negative normalized). ✓
- Reserved fields (SVLEN/SVTYPE/SVCLAIM/END/CIPOS/CIEND/CILEN/MATEID/PARID/IMPRECISE/CN/LEN) → Task 2. ✓
- Typed `Allele` union as public vocab + stored model type → Tasks 1, 3, 9. ✓
- Validity contract (static where possible; eager runtime otherwise; never silent) → Tasks 1 (constructors), 5 (builder), 5 tests assert raises. ✓
- Ground-truth contract (`is_sequence` per allele + mask; variant_class vocab) → Task 4. ✓
- `##ALT` auto-emission + `.alt()` override → Task 6. ✓
- Round-trip oracle extended to symbolic space → Task 7. ✓
- Migration of existing fixtures/tests → Task 3 Step 7. ✓
- Deferred items (full BND topology, CNV:TR family, SV FORMAT semantics, FORMAT/LEN, IUPAC) → not implemented; seams (`BreakendAllele.raw`, additive registry, exhaustive `match`, `AlleleTruth` dataclass) preserved. ✓

**Placeholder scan:** No "TBD"/"add validation"/"similar to". Every code step shows real code. Task 8 Step 2 flags a *contingency* (if cyvcf2 rejects a breakend without `SVTYPE`) with a concrete remedy, not a placeholder.

**Type consistency:** `classify_allele` (no `ref` arg) used consistently (Tasks 1, 3). `record_class(ref, tuple[Allele,...])` matches caller in `truth.py` (Task 3 Step 4). `AlleleTruth(kind, is_sequence, sv_type, svlen, sv_end)` field order identical in Tasks 4 def and tests. `type_str` property used in Tasks 1, 4, 6. `SymbolicAllele.first_type` checked against the same `{DEL,INS,DUP,INV,CNV}` set in `allele.py`, `variants.py`, `truth.py`, `build.py`. Aliases `Seq/Sym/Star/Unspecified/Bnd` defined in Task 1, used everywhere after. `VcfDocument.alt_defs` defaulted last (Task 6) so existing positional constructors in `test_model.py`/`test_truth.py` (migrated in Task 3) stay valid.
