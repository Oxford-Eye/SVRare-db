import gzip
import os
import glob
import copy
import csv
import argparse
from lib import models, gnomad, decipher, dbvar, utils, Interval_base, Types
from lib.patient import Patient, SV
import sqlalchemy as sa
from sqlalchemy.orm import Session
from typing import List
import pysam
import yaml

CHROMOSOMES = [str(i) for i in range(23)] + ['X', 'Y']
CHROMOSOMES += [f"chr{i}" for i in CHROMOSOMES]
CHROMOSOMES = set(CHROMOSOMES)

def get_duplicates(SVs, distance_cutoff):
    # mark duplicate on filtered, canvas call
    # since Canvas call tends to overestimate size
    bad_svs = set()
    SVs.sort(key=lambda x: (x.FILTER!='PASS', x.source == 'canvas'))
    for ind_i in range(len(SVs)-1):
        if SVs[ind_i].vcf_id in bad_svs:
            continue
        for ind_j in range(ind_i+1, len(SVs)):
            if SVs[ind_j].vcf_id in bad_svs:
                continue
            distance = Interval_base.get_distance(SVs[ind_i], SVs[ind_j])
            if distance <= distance_cutoff:
                bad_svs.add(SVs[ind_j].vcf_id)
    return bad_svs

def annotate(SVs: List[SV], freq_dbs: List, config) -> List:
    tbx_gtf = pysam.TabixFile(config['gene_tbx'])
    for sv in SVs:
        # external freqs. Note that only gnomad supports annotation for INV
        for freq_db in freq_dbs:
            if freq_db['name'] == 'dbvar':
                key = f"{freq_db['name']}_count"
                if sv.svtype == Types.SVtype.INV:
                    setattr(sv, key, None)
                else:
                    setattr(sv, key, (freq_db[sv.svtype.value].get_count(sv, 0.5)))
            else:
                key = f"{freq_db['name']}_freq"
                if sv.svtype == Types.SVtype.INV and freq_db['name'] != 'gnomad':
                    setattr(sv, key, None)
                else:
                    setattr(sv, key, (freq_db[sv.svtype.value].get_freq(sv, 0.5)))
        # genes
        genes = utils._get_overlap(sv, 'gene', tbx_gtf, True, False)
        setattr(sv, 'genes', genes)
        # cdss
        cdss = utils.get_protein_coding_disrupted_genes(sv, sv.svtype, tbx_gtf, feature='CDS')
        setattr(sv, 'cdss', cdss)
        # exons
        exons = utils.get_protein_coding_disrupted_genes(sv, sv.svtype, tbx_gtf, feature='exon')
        setattr(sv, 'exons', exons)
        
        
def main(config):
    engine = sa.create_engine(config['db'])
    session = Session(engine)
    freq_dbs = [
        {
            'name': 'gnomad',
            'LOSS': gnomad.Gnomad(config['gnomad'], 'LOSS'),
            'GAIN': gnomad.Gnomad(config['gnomad'], 'GAIN'),
            'INV': gnomad.Gnomad(config['gnomad'], 'INV'),
        },
        {
            'name': 'dbvar',
            'LOSS': dbvar.Dbvar(config['dbvar']['LOSS'], 'LOSS'),
            'GAIN': dbvar.Dbvar(config['dbvar']['GAIN'], 'GAIN'),
        },
        {
            'name': 'decipher',
            'LOSS': decipher.Decipher(config['decipher'], 'LOSS'),
            'GAIN': decipher.Decipher(config['decipher'], 'GAIN'),
        }
    ]
    
    # read patients
    patients = []
    with open(config['patients'], 'rt') as inf:
        csvreader = csv.reader(inf, delimiter='\t')
        header = []
        for row in csvreader:
            if not header:
                header = row
                continue
            patients.append(dict(zip(header, row)))
            
    # import patients
    entities = []
    for patient in patients:
        entities.append(models.Patient(
            name = patient['name'],
            family_id = patient['family_id'],
            manta_path = patient['manta_path'],
            canvas_path = patient['canvas_path'],
            bam_path = patient['bam_path'],
            is_solved = True if int(patient['is_solved']) != 0 else False,
            disease = patient['disease'],
            is_proband = True if int(patient['is_proband']) != 0 else False,
            relation_to_proband = patient['relation_to_proband'],
        ))
    
    session.add_all(entities)
    session.flush()
    
    # import patient_hpo
    entities = []
    for patient in patients:
        patient_id = None
        patient_id = session.query(models.Patient).filter(models.Patient.name == patient['name']).one().id
        if patient_id is None:
            raise ValueError(f"Patient not seen in db: {patient['name']}")
        patient['id'] = patient_id
        hpo_ids = ()
        if patient['HPO']:
            hpo_ids = (int(hpo.lstrip('HP:')) for hpo in patient['HPO'].split(','))
        for hpo_id in hpo_ids:
            entities.append(models.Patient_HPO(
                hpo_id = hpo_id,
                patient_id = patient_id
            ))
            
    session.add_all(entities)
    session.flush()
    
    
    # deal with SVs    
    for input_patient in patients:
        print(f"importing {input_patient['name']}")
        patient = Patient(
            name = input_patient['name'],
            canvas_file = input_patient['canvas_path'],
            manta_file = input_patient['manta_path'],
            bam_file = input_patient['bam_path'],
        )
        # get manta / canvas etc.
        SVs = patient.parse_vcf('manta') + patient.parse_vcf('canvas')
        # filter SVs on chromosomes
        SVs = list(filter(lambda x: x.chrom in CHROMOSOMES, SVs))
        # annotate SVs
        annotate(SVs, freq_dbs, config)
        # import SV
        for sv in SVs:
            # if sv already there?
            sv_id = None
            db_sv = session.query(models.SV).filter(models.SV.name == sv.name).first()
            if db_sv is not None:
                sv_id = db_sv.id
            else:
                sv_entity = models.SV(
                    name = sv.name,
                    chrom = sv.chrom,
                    start = sv.start,
                    end = sv.end,
                    sv_type = sv.svtype.value,
                    gnomad_freq = sv.gnomad_freq,
                    dbvar_count = sv.dbvar_count,
                    decipher_freq = sv.decipher_freq,
                )
                session.add(sv_entity)
                # session flush to get id
                session.flush()
                sv_id = sv_entity.id
                # SV_gene/CDS/exon
                for gene in sv.genes:
                    entity = models.SV_Gene(
                        sv_id = sv_id,
                        gene_id = int(gene.lstrip('ENSG')),
                    )
                    session.add(entity)
                for gene in sv.cdss:
                    entity = models.SV_CDS(
                        sv_id = sv_id,
                        gene_id = int(gene.lstrip('ENSG')),
                    )
                    session.add(entity)
                for gene in sv.exons:
                    entity = models.SV_Exon(
                        sv_id = sv_id,
                        gene_id = int(gene.lstrip('ENSG')),
                    )
                    session.add(entity)
                        
        # deal with duplicates
        groups = Interval_base.group(SVs)
        for group in groups:
            duplicate_svs = get_duplicates(group.intervals, config['params']['distance'])
            for sv in group.intervals:
                Patient_SV_entity = models.Patient_SV(
                    patient_id = input_patient['id'],
                    sv_id = sv_id,
                    genotype = sv.genotype.value,
                    vcf_id = sv.vcf_id,
                    source = sv.source,
                    filter = sv.FILTER,
                    is_duplicate = True if sv.vcf_id in duplicate_svs else False,
                )
                session.add(Patient_SV_entity)
        session.commit()
    session.commit()
        
        
if __name__ == '__main__':
    with open('config.yml', 'rt') as inf:
        config = yaml.safe_load(inf)
        main(config)
        