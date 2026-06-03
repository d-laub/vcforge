from vcfixture._spec.number import Number
from vcfixture._spec.types import Type
from vcfixture.allele import Seq, Sym
from vcfixture.build import VcfBuilder


def _doc():
    return (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 1000)])
        .fmt("GT")
        .info("SVLEN", Number.A, Type.INTEGER)
        .info("SVCLAIM", Number.A, Type.STRING)
        .alt("DEL", "Custom deletion description")
        .record(
            "chr1",
            100,
            ref="G",
            alt=[Sym.deletion()],
            gt=["0/1"],
            info={"SVLEN": [50], "SVCLAIM": ["DJ"]},
        )
        .record(
            "chr1",
            200,
            ref="G",
            alt=[Sym.duplication("TANDEM")],
            gt=["0/1"],
            info={"SVLEN": [20], "SVCLAIM": ["DJ"]},
        )
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
    data = [line for line in text.splitlines() if line.startswith("chr1\t100")][0]
    assert data.split("\t")[4] == "<DEL>"


def test_no_symbolic_alts_emits_no_alt_header():
    doc = (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 1000)])
        .fmt("GT")
        .record("chr1", 100, ref="G", alt=[Seq("A")], gt=["0/1"])
        .build()
    )
    assert not any(line.startswith("##ALT") for line in doc.render().splitlines())


def test_explicit_alt_def_without_matching_record():
    doc = (
        VcfBuilder(samples=["s1"], contigs=[("chr1", 1000)])
        .fmt("GT")
        .alt("INV", "Standalone inversion desc")
        .record("chr1", 100, ref="G", alt=[Seq("A")], gt=["0/1"])
        .build()
    )
    lines = doc.render().splitlines()
    assert '##ALT=<ID=INV,Description="Standalone inversion desc">' in lines
