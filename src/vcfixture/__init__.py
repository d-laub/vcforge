"""vcfixture — generate small VCF test data with decoded ground truth."""

from importlib import metadata

from . import strategies
from ._spec.number import Number
from ._spec.types import Type
from .build import VcfBuilder
from .genotype import Genotype
from .reference import Reference, ReferenceBuilder, ReferenceSpec, RepeatFeature
from .truth import GroundTruth

__version__ = metadata.version("vcfixture")
__all__ = [
    "VcfBuilder",
    "Genotype",
    "Reference",
    "ReferenceBuilder",
    "ReferenceSpec",
    "RepeatFeature",
    "GroundTruth",
    "Number",
    "Type",
    "strategies",
    "__version__",
]
