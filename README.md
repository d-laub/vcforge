# vcfixture

[![PyPI version](https://img.shields.io/pypi/v/vcfixture.svg)](https://pypi.org/project/vcfixture/)
[![Docs](https://img.shields.io/badge/docs-d--laub.github.io%2Fvcfixture-blue.svg)](https://d-laub.github.io/vcfixture/)

Generate small VCF (v4.x) test data — via an explicit builder or Hypothesis
strategies — with the decoded ground truth returned alongside, so parser tests
assert against a known oracle instead of hand-coded literals.

```python
from vcfixture import VcfBuilder, Number, Type, Seq

doc = (VcfBuilder(samples=["s1", "s2"], contigs=[("chr1", 100000)])
       .info("AF", Number.A, Type.FLOAT)
       .fmt("GT").fmt("DS", Number.A, Type.FLOAT)
       .record("chr1", 81262, ref="GAT", alt=[Seq("A")],  # ALTs are typed Allele instances
               gt=["0|1", "1|1"], info={"AF": [0.5]}, DS=[[1.0], [2.0]]))

text  = doc.render()                       # VCF text
truth = doc.truth()                        # GroundTruth (numpy genotypes, ...)
doc.write("x.vcf.gz", bgzip=True, index=True)
```

See the [full documentation](https://d-laub.github.io/vcfixture/) for guides and
API reference, and [CONTRIBUTING.md](CONTRIBUTING.md) for development setup (uv
is required).
