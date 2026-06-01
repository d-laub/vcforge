# Compact model reprs

**Date:** 2026-06-01
**Status:** Approved, pending implementation

## Problem

A failing Hypothesis property test prints a "Falsifying example" that is ~2,000
lines long, dominated by a deeply nested dump of `Number` objects:

```
number=Number(kind=<NumberKind.FIXED: 'FIXED'>, count=1,
 ONE=Number(...), A=Number(kind=<NumberKind.A: 'A'>, count=None,
  ONE=Number(...), R=Number(...), G=Number(...), ...
```

This makes test failures effectively unreadable.

## Root cause

`Number` declares its singletons (`ONE`, `A`, `R`, `G`, `DOT`, `FLAG`) as
`ClassVar[Number]` annotations **in the class body** (`_spec/number.py:25-30`).
Because the module uses `from __future__ import annotations`, dataclasses records
these as `_FIELD_CLASSVAR` entries in `Number.__dataclass_fields__` — and those
entries still carry `init=True`.

- Standard `repr()` correctly excludes `_FIELD_CLASSVAR` entries → compact output:
  `Number(kind=<NumberKind.G: 'G'>, count=None)`.
- Hypothesis's vendored pretty-printer (`hypothesis/vendor/pretty.py`) field-walks
  `__dataclass_fields__` filtering **only on `v.init`**, ignoring `_field_type`.
  It therefore recurses into every singleton, and each singleton recurses into all
  the others — exponential fan-out, rendered until it bails. This single leak is
  what explodes the entire `VcfDocument` repr.

Verified: the six singletons appear in `__dataclass_fields__` with
`_field_type=_FIELD_CLASSVAR, init=True`.

## Goals

1. Stop the recursive explosion (mandatory).
2. Make model-object reprs compact and human-readable in **all** contexts: plain
   `repr()`, pytest assertion diffs, and Hypothesis Falsifying examples.

Non-goals: changing equality/hashing, changing serialized VCF output, adding any
runtime dependency, or making reprs `eval`-able (the compact forms are
intentionally lossy — full detail remains available via `.render()` and field
access).

## Design

### 1. Root-cause fix (`_spec/number.py`)

Move the six `ClassVar[Number]` annotations into an `if TYPE_CHECKING:` block:

```python
from typing import TYPE_CHECKING, ClassVar

@dataclass(frozen=True)
class Number:
    kind: NumberKind
    count: int | None = None

    if TYPE_CHECKING:
        ONE: ClassVar[Number]
        A: ClassVar[Number]
        R: ClassVar[Number]
        G: ClassVar[Number]
        DOT: ClassVar[Number]
        FLAG: ClassVar[Number]
```

The runtime assignments (`Number.ONE = Number(...)`, etc., lines 61-66) are
unchanged. Because `TYPE_CHECKING` is `False` at runtime, the annotations are not
collected into the class, so the singletons leave `__dataclass_fields__` entirely
— any field-walker (Hypothesis, `dataclasses.asdict`, third-party reprs) stops
recursing. Type checkers still see the `ClassVar` declarations.

Verified: with this change `Number.__dataclass_fields__` contains only
`kind`/`count`; `pretty(Number.G)` prints `Number(kind=<NumberKind.G: 'G'>,
count=None)`; and `pyrefly check` reports 0 errors on the pattern.

This change alone resolves the reported 2k-line blowup.

### 2. Pretty-print mixin (`_repr.py`, new private module)

```python
from typing import Any

class CompactRepr:
    """Route Hypothesis's pretty-printer (and IPython's) through __repr__."""
    def _repr_pretty_(self, p: Any, cycle: bool) -> None:
        p.text(repr(self))
```

`_repr_pretty_` is the documented IPython pretty protocol; Hypothesis's vendored
pretty-printer checks for it **before** the dataclass field-walk (verified by
reading `pretty.py`). It requires **no Hypothesis import**, so the clean-room
wheel smoke test stays green (`p` is typed `Any`; `p.text(...)` is valid on `Any`
under pyrefly strict).

**Important inheritance rule:** the mixin provides *only* `_repr_pretty_`. Each
model dataclass must define its compact `__repr__` **in its own class body**.
`@dataclass` uses `_set_new_attribute`, which does not overwrite a `__repr__`
already present in the class `__dict__` — so a body-level `__repr__` is respected.
A `__repr__` inherited from the mixin would NOT be (dataclass would generate one
on the subclass and shadow it), which is why `__repr__` lives per-class and only
`_repr_pretty_` lives on the mixin.

### 3. Compact `__repr__` per class

Each frozen dataclass in the model layer inherits `CompactRepr` and defines
`__repr__`:

| Class | Module | Compact form |
|---|---|---|
| `Number` | `_spec/number.py` | `Number(G)`, `Number(2)`, `Number(FLAG)` |
| `FieldDef` | `_spec/fielddef.py` | `FieldDef(GT FORMAT Number=1 Type=String)` |
| `Genotype` | `genotype.py` | `Genotype(0\|1)` |
| `ContigDef` | `model.py` | `ContigDef(chr1:200)` / `ContigDef(chr1)` if length is None |
| `RepeatFeature` | `reference.py` | `RepeatFeature(chr1@10 AT×3)` |
| `Record` | `model.py` | `Record(chr1:5 A>T,G ×2 [snp,multiallelic])` |
| `ReferenceSpec` | `reference.py` | `ReferenceSpec(contigs=[chr1:200bp], repeats=2)` |
| `VcfDocument` | `model.py` | `VcfDocument(VCFv4.5 samples=2 records=5 info=0 format=1)` |

**Format details:**

- **`Number`** — token: `FIXED` → `str(count)`; `FLAG` → `FLAG`; otherwise
  `kind.value` (`A`/`R`/`G`/`.`). Rendered as `Number(<token>)`.
- **`FieldDef`** — `FieldDef(<id> <kind> Number=<number.header_str()> Type=<type.value>)`.
  Description is omitted (often long, not identifying).
- **`Genotype`** — `Genotype(<render()>)`, e.g. `Genotype(0|1)`, `Genotype(./.)`.
- **`ContigDef`** — `ContigDef(<id>:<length>)`, or `ContigDef(<id>)` when
  `length is None`.
- **`RepeatFeature`** — `RepeatFeature(<contig>@<pos0> <motif>×<count>)`.
- **`Record`** — `Record(<chrom>:<pos> <ref>><alts joined by ','> ×<n samples>[ <labels>])`.
  Alts joined with `,`; if `alts` is empty, show `.` (i.e. `A>.`). Sample count is
  `len(samples)`. Labels: sorted, `[a,b,c]`; the whole ` [...]` segment is omitted
  when `labels` is empty. `ids`, `qual`, `filters`, and `info` are omitted from the
  compact form (available via fields / `.render()`).
- **`ReferenceSpec`** — `ReferenceSpec(contigs=[<id>:<len(seq)>bp, ...], repeats=<n>)`.
  The full sequence strings are dropped (this is the second-biggest noise source
  after `Number`).
- **`VcfDocument`** — `VcfDocument(<fileformat> samples=<n> records=<n>
  info=<len(info_defs)> format=<len(format_defs)>)`.

`Reference` and `ReferenceBuilder` are plain classes (not dataclasses) and are not
part of this change.

### 4. Tests (`tests/`, outside the type-check scope)

- **Regression guard:** render a representative `VcfDocument` (and a
  `ReferenceSpec`) through Hypothesis's vendored `pretty()` and assert the output
  (a) is under a sane length bound (e.g. < 2,000 chars) and (b) contains no
  `ONE=Number` / `_FIELD_CLASSVAR`-recursion marker. This directly pins the
  reported bug.
- **Custom-repr-is-used guard:** assert `repr(Number.G) == "Number(G)"` and that
  `repr` of a sample `VcfDocument`/`Record`/etc. matches the expected compact form
  — this also catches an accidental `@dataclass`-generated `__repr__` shadowing the
  hand-written one.
- **Per-class format assertions:** one assertion per class covering the formats in
  the table, including the edge cases (`ContigDef` with `length=None`, `Record`
  with empty `labels`, `Record` with empty `alts`).

## Files touched

- `src/vcfixture/_spec/number.py` — TYPE_CHECKING block + `CompactRepr` + `__repr__`.
- `src/vcfixture/_repr.py` — **new**, the `CompactRepr` mixin.
- `src/vcfixture/_spec/fielddef.py` — `CompactRepr` + `__repr__`.
- `src/vcfixture/genotype.py` — `CompactRepr` + `__repr__`.
- `src/vcfixture/model.py` — `CompactRepr` + `__repr__` on `ContigDef`, `Record`,
  `VcfDocument`.
- `src/vcfixture/reference.py` — `CompactRepr` + `__repr__` on `RepeatFeature`,
  `ReferenceSpec`.
- `tests/test_repr.py` — **new**, the tests above.

## Verification

- `uv run pytest` — new repr tests pass; existing suite unaffected.
- `uv run pyrefly check` — strict, 0 errors (the `_repr_pretty_` `Any` param and
  the TYPE_CHECKING ClassVars are both clean).
- `uv run ruff check && uv run ruff format`.
- Clean-room wheel smoke test (CI) — unaffected; no runtime import added.
