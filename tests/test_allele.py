from __future__ import annotations

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
