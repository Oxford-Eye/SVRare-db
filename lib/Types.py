from typing import Type
import attr
from enum import Enum

class Output_format(Enum):
    """
    if want the value:
    Output_format.JSON.value
    """
    JSON = 'json'
    TSV = 'tsv'
    CSV = 'csv'

class Genotype(Enum):
    HOM = 'hom'
    HET = 'het'
    REF = 'ref'
    UNKNOWN = None

class SVtype(Enum):
    LOSS = 'LOSS'
    GAIN = 'GAIN'
    INV = 'INV'
    INS = 'INS'
    REF = 'REF'
    UNKNOWN = None

@attr.s(frozen=True)
class Cutoffs:
    internal: float = attr.ib()
    gnomad: float = attr.ib()
    # dbvar doesn't have freq. Largest number is 7327, so use 73 as 0.01, roughly.
    dbvar: int = attr.ib()
    decipher: float = attr.ib()

@attr.s(frozen=True)
class Params:
    cutoffs: Cutoffs = attr.ib()
    distance: float = attr.ib()
    output_format: Output_format = attr.ib()

