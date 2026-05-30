from __future__ import annotations
from hypothesis import strategies as st
from .build import VcfBuilder
from .genotype import Genotype
from .variants import (snp, mnp, insertion, deletion, delins, spanning_deletion)
from ._spec.number import Number, NumberKind
from ._spec.types import Type

def _build_number_type_combos():
    numbers = [Number.ONE, Number.fixed(2), Number.A, Number.R, Number.G,
               Number.DOT]
    combos = []
    for kind in ("INFO", "FORMAT"):
        allowed = Type.info_allowed() if kind == "INFO" else Type.format_allowed()
        for n in numbers:
            for t in allowed:
                if t is Type.FLAG:
                    continue
                combos.append((n, t, kind))
        if kind == "INFO":
            combos.append((Number.FLAG, Type.FLAG, "INFO"))
    return combos

ALL_NUMBER_TYPE_COMBOS = _build_number_type_combos()
ALL_VARIANT_CLASSES = ["SNP", "MNP", "INS", "DEL", "DELINS", "SPANNING_DEL"]

_BASES = "ACGT"

@st.composite
def _ref_alt(draw, klass: str):
    b = draw(st.sampled_from(_BASES))
    b2 = draw(st.sampled_from(_BASES))
    if klass == "SNP":
        alt = _BASES[(_BASES.index(b) + 1 + draw(st.integers(0, 2))) % 4]
        return snp(b, alt)
    if klass == "MNP":
        return mnp(b + b2, _BASES[(_BASES.index(b) + 1) % 4]
                   + _BASES[(_BASES.index(b2) + 1) % 4])
    if klass == "INS":
        return insertion(b, draw(st.text(_BASES, min_size=1, max_size=3)))
    if klass == "DEL":
        return deletion(b, draw(st.text(_BASES, min_size=1, max_size=3)))
    if klass == "DELINS":
        return delins(b + b2, draw(st.text(_BASES, min_size=1, max_size=3)))
    return spanning_deletion(b)

@st.composite
def genotypes(draw, ploidy: int, n_alt: int, missing_rate: float = 0.1):
    alleles = []
    for _ in range(ploidy):
        if draw(st.floats(0, 1)) < missing_rate:
            alleles.append(".")
        else:
            alleles.append(str(draw(st.integers(0, n_alt))))
    phased = draw(st.booleans())
    sep = "|" if phased else "/"
    return sep.join(alleles)

@st.composite
def documents(draw, max_samples: int = 3, max_records: int = 4):
    n_samples = draw(st.integers(1, max_samples))
    samples = [f"s{i}" for i in range(n_samples)]
    ploidy = draw(st.integers(1, 2))
    b = VcfBuilder(samples=samples, contigs=[("chr1", 100000)]).fmt("GT")

    n_rec = draw(st.integers(1, max_records))
    pos = 1000
    for _ in range(n_rec):
        klass = draw(st.sampled_from(ALL_VARIANT_CLASSES))
        ref, alt = draw(_ref_alt(klass))
        gts = [draw(genotypes(ploidy, n_alt=1)) for _ in samples]
        b.record("chr1", pos, ref=ref, alt=[alt], gt=gts)
        pos += draw(st.integers(1, 50))
    return b.build()
