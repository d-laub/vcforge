"""Performance gate: drawing + serializing + deriving truth for a
representative document must stay cheap."""

import time

from vcfixture._spec.number import Number
from vcfixture._spec.types import Type
from vcfixture.allele import Seq, Star, Sym
from vcfixture.build import VcfBuilder


def _representative_doc():
    b = VcfBuilder(
        samples=[f"s{i}" for i in range(10)], contigs=[("chr1", 1_000_000)]
    ).fmt("GT")
    for k in range(50):
        gts = ["0|1" if (i + k) % 2 else "1/1" for i in range(10)]
        b.record("chr1", 1000 + k * 10, ref="A", alt=[Seq("T")], gt=gts)
    return b.build()


def _representative_symbolic_doc():
    # 3 samples, 4 SV records covering the main symbolic SV types + <*>.
    # Scale is modest (4 records × 3 samples) because SV validation is heavier
    # per-record than plain sequence alleles; the budget matches the sibling.
    b = (
        VcfBuilder(samples=["s0", "s1", "s2"], contigs=[("chr1", 10_000_000)])
        .fmt("GT")
        .info("SVLEN", Number.A, Type.INTEGER)
        .info("SVCLAIM", Number.A, Type.STRING)
    )
    b.record(
        "chr1",
        1000,
        ref="A",
        alt=[Sym.deletion()],
        gt=["0|1", "1/1", "0/0"],
        info={"SVLEN": [500], "SVCLAIM": ["DJ"]},
    )
    b.record(
        "chr1",
        5000,
        ref="G",
        alt=[Sym.duplication("TANDEM")],
        gt=["0/1", "0|0", "1|1"],
        info={"SVLEN": [1000], "SVCLAIM": ["D"]},
    )
    b.record(
        "chr1",
        10000,
        ref="C",
        alt=[Sym.insertion()],
        gt=["1/1", "0|1", "0/0"],
        info={"SVLEN": [250], "SVCLAIM": ["J"]},
    )
    # <*> (spanning deletion) needs no SVLEN/SVCLAIM
    b.record(
        "chr1",
        15000,
        ref="T",
        alt=[Star()],
        gt=["0/1", "0/1", "0/0"],
    )
    return b.build()


def test_build_serialize_truth_under_budget():
    t0 = time.perf_counter()
    for _ in range(100):
        doc = _representative_doc()
        doc.render()
        doc.truth()
    elapsed = time.perf_counter() - t0
    assert elapsed < 5.0, f"too slow: {elapsed:.2f}s for 100 iterations"


def test_symbolic_draw_serialize_truth_under_budget():
    # Measures render+truth on a directly-built representative symbolic doc —
    # no Hypothesis overhead. Budget matches the sibling: both do render+truth
    # at comparable scale (the symbolic doc is smaller in records but heavier
    # per record, so the overall wall time is similar).
    t0 = time.perf_counter()
    for _ in range(100):
        doc = _representative_symbolic_doc()
        doc.render()
        doc.truth()
    elapsed = time.perf_counter() - t0
    assert elapsed < 5.0, f"too slow: {elapsed:.2f}s for 100 iterations"
