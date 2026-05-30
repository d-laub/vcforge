import vcforge

def test_public_exports():
    for name in ["VcfBuilder", "Number", "Type", "Genotype",
                 "GroundTruth", "Reference", "strategies"]:
        assert hasattr(vcforge, name), f"missing public export: {name}"
