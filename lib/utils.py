import gzip
import os
from typing import List, Set
from lib import Params, Types, Interval_base, Interval_base


def read_cnv_file(infile):
    '''
    Read single sample CNV file
    '''
    header = [
        'chrom',
        'start',
        'id',
        'ref',
        'alt',
        'qual',
        'filter',
        'info',
        'format',
        'genotype',
    ]
    result = []
    with gzip.open(infile, 'rt') as inf:
        for line in inf:
            if line.startswith('#'):
                pass
            row = line.rstrip().split('\t')
            row_dict = dict(zip(header, row))
            row_dict['start'] = int(row_dict['start'])
            id_details = row_dict['id'].split(':')
            row_dict['type'] = id_details[1]
            row_dict['end'] = int(id_details[2].split('-')[1])
            if row_dict['type'] != 'REF':
                result.append(row_dict)
    return result

def _get_overlap(interval:Interval_base, sequence_type, tbx_gtf, gene_id:bool, protein_coding_gene:bool) -> Set[str]:
    result = set()
    chrom = interval.chrom.lstrip('chr')
    if chrom == 'M':
        chrom = 'MT'
    for line in tbx_gtf.fetch(chrom, interval.start, interval.end):
        row = line.rstrip().split('\t')
        if row[2] != sequence_type:
            continue
        info = {}
        for field in row[-1].split('; '):
            key, value = field.split(' ', 1)
            info[key] = value.split('"')[1]
        if protein_coding_gene and info.get('transcript_biotype', info['gene_biotype']) != 'protein_coding':
            continue
        if gene_id:
            result.add(info['gene_id'])
        else:
            result.add(info['gene_name'])
    return result

def get_protein_coding_disrupted_genes(interval:Interval_base, interval_type:Types.SVtype, tbx_gtf, gene_id = True, feature = 'exon', protein_coding_gene = True) -> List[str]:
    '''
    For CDS/exon and Genes
    INV: to be protein-coding disrupting, it has to cover at least one exon, 
        and at least one of its BND is within the said gene
    GAIN: has to be intragenic, and cover at least one exon.
    '''
    if feature not in ('exon', 'CDS'):
        raise ValueError(f"feature has to be iether exon or CDS, got {feature}")
    covered_feature = _get_overlap(interval, feature, tbx_gtf, gene_id, protein_coding_gene)
    if covered_feature:
        if interval_type == Types.SVtype.LOSS:
            return list(covered_feature)
        start_genes = _get_overlap(Interval_base(
            interval.chrom,
            interval.start,
            interval.start + 1,
            ), 'gene', tbx_gtf, gene_id, protein_coding_gene)
        end_genes = _get_overlap(Interval_base(
            interval.chrom,
            interval.end,
            interval.end + 1,
            ), 'gene', tbx_gtf, gene_id, protein_coding_gene)
        boundary_genes = set()
        if interval_type == Types.SVtype.INV:
            boundary_genes = start_genes | end_genes
        elif interval_type == Types.SVtype.GAIN:
            boundary_genes = start_genes & end_genes
        return list(covered_feature & boundary_genes)
    return []
