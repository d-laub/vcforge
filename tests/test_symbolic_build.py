from __future__ import annotations

import pytest

from vcfixture._spec.number import Number
from vcfixture._spec.types import Type
from vcfixture.allele import Sym, Unspecified
from vcfixture.build import VcfBuilder


def _b() -> VcfBuilder:
    return (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 1000)])
        .fmt("GT")
        .info("SVLEN", Number.A, Type.INTEGER)
        .info("SVCLAIM", Number.A, Type.STRING)
    )


def test_valid_symbolic_del_builds():
    doc = (
        _b()
        .record(
            "chr1",
            100,
            ref="G",
            alt=[Sym.deletion()],
            gt=["0/1"],
            info={"SVLEN": [50], "SVCLAIM": ["DJ"]},
        )
        .build()
    )
    assert doc.records[0].alts[0].render() == "<DEL>"


def test_symbolic_requires_single_base_ref():
    with pytest.raises(ValueError, match="padding base"):
        _b().record(
            "chr1",
            100,
            ref="GA",
            alt=[Sym.deletion()],
            gt=["0/1"],
            info={"SVLEN": [50], "SVCLAIM": ["DJ"]},
        )


def test_symbolic_sv_requires_svlen():
    with pytest.raises(ValueError, match="SVLEN"):
        _b().record(
            "chr1",
            100,
            ref="G",
            alt=[Sym.deletion()],
            gt=["0/1"],
            info={"SVCLAIM": ["DJ"]},
        )


def test_svclaim_del_must_be_djd():
    with pytest.raises(ValueError, match="invalid for"):
        _b().record(
            "chr1",
            100,
            ref="G",
            alt=[Sym.deletion()],
            gt=["0/1"],
            info={"SVLEN": [50], "SVCLAIM": ["X"]},
        )


def test_unspecified_allele_does_not_require_padding_or_svlen():
    doc = _b().record("chr1", 100, ref="G", alt=[Unspecified()], gt=["0/1"]).build()
    assert doc.records[0].alts[0].render() == "<*>"


def test_breakend_svlen_must_be_missing():
    from vcfixture.allele import Bnd

    with pytest.raises(ValueError, match="missing"):
        _b().record(
            "chr1",
            100,
            ref="G",
            alt=[Bnd.parse("G[chr2:9[")],
            gt=["0/1"],
            info={"SVLEN": [10]},
        )


# ---------------------------------------------------------------------------
# Fix 1: scalar SVLEN normalization — builder/oracle divergence fix
# ---------------------------------------------------------------------------


def test_scalar_svlen_accepted_for_single_alt():
    """Scalar SVLEN (not wrapped in a list) must be treated as index-0 value."""
    doc = (
        _b()
        .record(
            "chr1",
            100,
            ref="G",
            alt=[Sym.deletion()],
            gt=["0/1"],
            info={"SVLEN": 50, "SVCLAIM": ["DJ"]},
        )
        .build()
    )
    assert doc.records[0].alts[0].render() == "<DEL>"


# ---------------------------------------------------------------------------
# Fix 2: valid SVCLAIM values for INS, INV, CNV
# ---------------------------------------------------------------------------


def _b_ins() -> VcfBuilder:
    return (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 1000)])
        .fmt("GT")
        .info("SVLEN", Number.A, Type.INTEGER)
        .info("SVCLAIM", Number.A, Type.STRING)
    )


def test_ins_svclaim_j_accepted():
    doc = (
        _b_ins()
        .record(
            "chr1",
            100,
            ref="G",
            alt=[Sym.insertion()],
            gt=["0/1"],
            info={"SVLEN": [100], "SVCLAIM": ["J"]},
        )
        .build()
    )
    assert doc.records[0].alts[0].render() == "<INS>"


def test_inv_svclaim_j_accepted():
    doc = (
        _b_ins()
        .record(
            "chr1",
            100,
            ref="G",
            alt=[Sym.inversion()],
            gt=["0/1"],
            info={"SVLEN": [500], "SVCLAIM": ["J"]},
        )
        .build()
    )
    assert doc.records[0].alts[0].render() == "<INV>"


def test_cnv_svclaim_d_accepted():
    doc = (
        _b_ins()
        .record(
            "chr1",
            100,
            ref="G",
            alt=[Sym.cnv()],
            gt=["0/1"],
            info={"SVLEN": [200], "SVCLAIM": ["D"]},
        )
        .build()
    )
    assert doc.records[0].alts[0].render() == "<CNV>"


# ---------------------------------------------------------------------------
# Fix 3: CN FORMAT guard — equal vs. unequal SVLEN
# ---------------------------------------------------------------------------


def _b_cn() -> VcfBuilder:
    """Builder with GT + CN FORMAT and SVLEN + SVCLAIM INFO."""
    return (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 10_000)])
        .fmt("GT")
        .fmt("CN")  # reserved: Number=1, Type=FLOAT
        .info("SVLEN", Number.A, Type.INTEGER)
        .info("SVCLAIM", Number.A, Type.STRING)
    )


def test_cn_equal_svlen_two_alleles_passes():
    """Two CN-eligible alleles with EQUAL SVLEN must build without error."""
    doc = (
        _b_cn()
        .record(
            "chr1",
            100,
            ref="G",
            alt=[Sym.deletion(), Sym.duplication()],
            gt=["0/1"],
            info={"SVLEN": [200, 200], "SVCLAIM": ["DJ", "DJ"]},
            CN=[2.0],
        )
        .build()
    )
    assert len(doc.records[0].alts) == 2


def test_cn_unequal_svlen_two_alleles_raises():
    """Two CN-eligible alleles with DIFFERENT SVLENs must raise ValueError."""
    with pytest.raises(ValueError, match="equal SVLEN"):
        _b_cn().record(
            "chr1",
            100,
            ref="G",
            alt=[Sym.deletion(), Sym.duplication()],
            gt=["0/1"],
            info={"SVLEN": [100, 200], "SVCLAIM": ["DJ", "DJ"]},
            CN=[2.0],
        )


# ---------------------------------------------------------------------------
# Consistency guard: _SVCLAIM_RULES key set == _SV_FIRST_TYPES
# ---------------------------------------------------------------------------


def test_svclaim_rules_cover_all_symbolic_types():
    from vcfixture.allele import _SV_FIRST_TYPES
    from vcfixture.build import _SVCLAIM_RULES

    assert set(_SVCLAIM_RULES) == set(_SV_FIRST_TYPES)
