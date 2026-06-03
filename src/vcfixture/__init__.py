"""vcfixture — generate small VCF test data with decoded ground truth."""

from importlib import metadata

from . import strategies
from ._spec.number import Number
from ._spec.types import Type
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
from .build import VcfBuilder
from .genotype import Genotype
from .reference import Reference, ReferenceBuilder, ReferenceSpec, RepeatFeature
from .truth import GroundTruth

__version__ = metadata.version("vcfixture")
__all__ = [
    "Allele",
    "Bnd",
    "BreakendAllele",
    "Genotype",
    "GroundTruth",
    "Number",
    "Reference",
    "ReferenceBuilder",
    "ReferenceSpec",
    "RepeatFeature",
    "Seq",
    "SequenceAllele",
    "SpanningDeletion",
    "Star",
    "Sym",
    "SymbolicAllele",
    "Type",
    "Unspecified",
    "UnspecifiedAllele",
    "VcfBuilder",
    "strategies",
    "__version__",
]
