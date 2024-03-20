import attr
from pathlib import Path
import os
from collections import Counter
import numpy as np
from cyvcf2 import VCF
from typing import Type, Tuple, List
from lib import Params, Types, Interval_base

class SV(Interval_base):
    def __init__(self, chrom, start, end, svtype, FILTER, source, vcf_id, genotype):
        super().__init__(chrom, start, end)
        self.name = f"{chrom}-{start}-{end}-{svtype.value}"
        self.svtype: Types.SVtype = svtype
        self.FILTER: str = FILTER
        self.source: str = source
        self.vcf_id: str = vcf_id
        self.genotype: Types.Genotype = genotype

@attr.s(frozen=True)
class Patient:
    """
    Class to hold patient details on manta/canvas variants
    """
    name: str = attr.ib()
    canvas_file: str = attr.ib()
    manta_file: str = attr.ib()
    svtools_file: str = attr.ib()
    pbsv_file: str = attr.ib()
    bam_file: str = attr.ib()

    def get_SV(self) -> List[SV]:
        svs = self.parse_vcf('manta')
        svs.extend(self.parse_vcf('canvas'))
        return svs
    
    def parse_vcf(self, what:str) -> List[SV]:
        """
        Parse Manta/Canvas file
        
        TODO: raise warning if finding no INV from manta
        """
        vcf = []
        what = what.lower()
        if what == 'manta':
            if self.manta_file is None or not os.path.isfile(self.manta_file):
                return []
            vcf = VCF(self.manta_file)
        elif what == 'canvas':
            if self.canvas_file is None or not os.path.isfile(self.canvas_file):
                return []
            vcf = VCF(self.canvas_file)
        elif what == 'svtools':
            if self.svtools_file is None or not os.path.isfile(self.svtools_file):
                return []
            vcf =  VCF(self.svtools_file)
        elif what == 'pbsv':
            if self.pbsv_file is None or not os.path.isfile(self.pbsv_file):
                return []
            vcf =  VCF(self.pbsv_file)
        else:
            raise ValueError(f"can only parse manta, canvas, svtools, pbsv file. Given {what}")
        sample_ind = vcf.samples.index(self.name)
        svs = []
        for variants in vcf:
            if not variants.ALT:
                continue
            for variant_ind, variant in enumerate(variants.ALT):
                svtype = translate_svtype(variants, variant)
                genotype = translate_genotype(tuple(variants.genotypes[sample_ind]))
                if svtype not in (Types.SVtype.LOSS, Types.SVtype.GAIN, Types.SVtype.INV):
                    continue
                if genotype not in (Types.Genotype.HOM, Types.Genotype.HET):
                    continue
                svs.append(SV(
                    variants.CHROM, 
                    int(variants.POS), 
                    int(variants.INFO.get('END')), 
                    svtype, 
                    variants.FILTER,
                    what,
                    variants.ID,
                    genotype
                    ))
        return svs
    

def translate_svtype(variants, variant) -> Types.SVtype:
    """
    Translate sv type.
    Manta: extract svtype from ID
    Canvas: extract svtype from variant
    """
    ID = variants.ID
    if ID.startswith('Canvas'):
        content = variant[1:-1]
        if content == 'DUP':
            return Types.SVtype.GAIN
        count = int(content[2:])
        if count < 1:
            return Types.SVtype.LOSS
        if count > 1:
            return Types.SVtype.GAIN
        return Types.SVtype.UNKNOWN
    if ID.startswith('Manta'):
        svtype = ID.split(':')[0][5:]
        if svtype == 'DEL':
            return Types.SVtype.LOSS
        if svtype == 'INS':
            return Types.SVtype.INS
        if svtype == 'DUP':
            return Types.SVtype.GAIN
        if svtype == 'INV':
            return Types.SVtype.INV
        return Types.SVtype.UNKNOWN
    svtype = variants.INFO.get('SVTYPE', None)
    if svtype is not None:
        if svtype == 'DEL':
            return Types.SVtype.LOSS
        if svtype == 'INS':
            return Types.SVtype.INS
        if svtype == 'DUP':
            return Types.SVtype.GAIN
        if svtype == 'INV':
            return Types.SVtype.INV
        return Types.SVtype.UNKNOWN
    return Types.SVtype.UNKNOWN
    
def translate_genotype(genotype: Tuple[int, int, bool]) -> Types.Genotype:
    """
    Translate genotype
    """
    ct = Counter(genotype[:2])
    if ct[1] == 2:
        return Types.Genotype.HOM
    if ct[1] == 1:
        return Types.Genotype.HET
    if ct[0] == 2:
        return Types.Genotype.REF
    return Types.Genotype.UNKNOWN



        
