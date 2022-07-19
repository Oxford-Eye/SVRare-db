import sys
import os
import yaml
import sqlalchemy as sa
from sqlalchemy.orm import declarative_base

Base = declarative_base()

'''
No need for relations / back_polulate.
But in case it does, refer to https://docs.sqlalchemy.org/en/14/orm/basic_relationships.html#association-object
'''
class Patient(Base):
    __tablename__ = 'Patient'
    
    id = sa.Column(sa.Integer, primary_key=True, autoincrement=True)
    name = sa.Column(sa.String, index=True, unique=True)
    family_id = sa.Column(sa.String, index=True)
    manta_path = sa.Column(sa.String)
    canvas_path = sa.Column(sa.String)
    bam_path = sa.Column(sa.String)
    svtools_path = sa.Column(sa.String)
    is_solved = sa.Column(sa.Boolean)
    disease = sa.Column(sa.String)
    is_proband = sa.Column(sa.Boolean)
    relation_to_proband = sa.Column(sa.String)
    
class SV(Base):
    __tablename__ = 'SV'
    
    id = sa.Column(sa.Integer, primary_key=True, autoincrement=True)
    # name: chrom-start-end-sv_type
    name = sa.Column(sa.String, index=True)
    chrom = sa.Column(sa.String, index=True)
    start = sa.Column(sa.Integer, index=True)
    end = sa.Column(sa.Integer, index=True)
    sv_type = sa.Column(sa.String, index=True)
    N_carriers = sa.Column(sa.Integer, index=True)
    gnomad_freq = sa.Column(sa.Numeric, index=True)
    dbvar_count = sa.Column(sa.Integer, index=True)
    decipher_freq = sa.Column(sa.Numeric, index=True)
    
class Patient_SV(Base):
    __tablename__ = 'Patient_SV'
    
    id = sa.Column(sa.Integer, primary_key=True, autoincrement=True)
    patient_id = sa.Column(sa.Integer, sa.ForeignKey("Patient.id"), nullable=False)
    sv_id = sa.Column(sa.Integer, sa.ForeignKey("SV.id"), nullable=False)
    genotype = sa.Column(sa.String, nullable=False)
    vcf_id = sa.Column(sa.String, nullable=False)
    source = sa.Column(sa.String, nullable=False)
    filter = sa.Column(sa.String, index=True)
    is_duplicate = sa.Column(sa.Boolean)
    igv_real = sa.Column(sa.Boolean, index=True)
    validated_as_real = sa.Column(sa.Boolean, index=True)
    
class Gene(Base):
    __tablename__ = 'Gene'
    
    # id is ensembl id without the ENSG part
    id = sa.Column(sa.Integer, primary_key=True)
    symbol =  sa.Column(sa.String, index=True)
    pli = sa.Column(sa.Integer, index=True)
    prec = sa.Column(sa.Integer, index=True)
    oe_lof_upper = sa.Column(sa.Integer, index=True)
    chrom = sa.Column(sa.String, nullable=False)
    start = sa.Column(sa.Integer, index=True)
    end = sa.Column(sa.Integer, index=True)

class SV_Gene(Base):
    __tablename__ = 'SV_Gene'
    
    id = sa.Column(sa.Integer, primary_key=True, autoincrement=True)
    sv_id = sa.Column(sa.Integer, sa.ForeignKey("SV.id"), nullable=False)
    gene_id = sa.Column(sa.Integer, sa.ForeignKey("Gene.id"), nullable=False)

class SV_Exon(Base):
    __tablename__ = 'SV_exon'
    
    id = sa.Column(sa.Integer, primary_key=True, autoincrement=True)
    sv_id = sa.Column(sa.Integer, sa.ForeignKey("SV.id"), nullable=False)
    gene_id = sa.Column(sa.Integer, sa.ForeignKey("Gene.id"), nullable=False)

class SV_CDS(Base):
    __tablename__ = 'SV_cds'
    
    id = sa.Column(sa.Integer, primary_key=True, autoincrement=True)
    sv_id = sa.Column(sa.Integer, sa.ForeignKey("SV.id"), nullable=False)
    gene_id = sa.Column(sa.Integer, sa.ForeignKey("Gene.id"), nullable=False)

class HPO(Base):
    __tablename__ = 'HPO'
    
    # id is HPO id without the HP: part
    id = sa.Column(sa.Integer, primary_key=True)
    name = sa.Column(sa.String, index=True)
    definition = sa.Column(sa.String)
    comment = sa.Column(sa.String)

class HPO_HPO(Base):
    __tablename__ = 'HPO_HPO'
    
    id = sa.Column(sa.Integer, primary_key=True, autoincrement=True)
    hpo_id = sa.Column(sa.Integer, sa.ForeignKey("HPO.id"), nullable=False)
    parent_hpo_id = sa.Column(sa.Integer, sa.ForeignKey("HPO.id"), nullable=False)

class HPO_Gene(Base):
    __tablename__ = 'HPO_Gene'
    
    id = sa.Column(sa.Integer, primary_key=True, autoincrement=True)
    hpo_id = sa.Column(sa.Integer, sa.ForeignKey("HPO.id"), nullable=False)
    gene_id = sa.Column(sa.Integer, sa.ForeignKey("Gene.id"), nullable=False)

class Patient_HPO(Base):
    __tablename__ = 'Patient_HPO'
    
    id = sa.Column(sa.Integer, primary_key=True, autoincrement=True)
    hpo_id = sa.Column(sa.Integer, sa.ForeignKey("HPO.id"), nullable=True)
    patient_id = sa.Column(sa.Integer, sa.ForeignKey("Patient.id"), nullable=False)

