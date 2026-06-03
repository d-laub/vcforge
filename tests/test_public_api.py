import vcfixture


def test_public_exports():
    for name in [
        "VcfBuilder",
        "Number",
        "Type",
        "Genotype",
        "GroundTruth",
        "Reference",
        "ReferenceBuilder",
        "ReferenceSpec",
        "RepeatFeature",
        "strategies",
    ]:
        assert hasattr(vcfixture, name), f"missing public export: {name}"


def test_allele_vocabulary_exported():
    for name in [
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
    ]:
        assert hasattr(vcfixture, name), f"missing public export: {name}"

    # Verify short aliases point at the correct classes
    assert vcfixture.Seq is vcfixture.SequenceAllele
    assert vcfixture.Sym is vcfixture.SymbolicAllele
    assert vcfixture.Star is vcfixture.SpanningDeletion
    assert vcfixture.Unspecified is vcfixture.UnspecifiedAllele
    assert vcfixture.Bnd is vcfixture.BreakendAllele
