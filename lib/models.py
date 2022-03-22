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
class Patient_table(Base):
    __tablename__ = 'patient'
    
    id = sa.Column(sa.Integer, primary_key=True)
    name = sa.Column(sa.String, index=True, unique=True)
    manta_vcf = sa.Column(sa.String)
    canvas_vcf = sa.Column(sa.String)
    bam = sa.Column(sa.String)
    is_solved = sa.Column(sa.Boolean)
    disease = sa.Column(sa.String)
    
    
class Interval_table(Base):
    __tablename__ = 'interval'
    
    id = sa.Column(sa.Integer, primary_key=True)
    sv_id = sa.Column(sa.String, index=True)
    chrom = sa.Column(sa.String, index=True)
    start = sa.Column(sa.Integer, index=True)
    end = sa.Column(sa.Integer, index=True)
    sv_type = sa.Column(sa.String, index=True)
    internal_freq = sa.Column(sa.Numeric, index=True)
    gnomad_freq = sa.Column(sa.Numeric, index=True)
    dbvar_count = sa.Column(sa.Integer, index=True)
    decipher_freq = sa.Column(sa.Numeric, index=True)
    
class Patient_Interval_table(Base):
    __tablename__ = 'patient_interval'
    
    id = sa.Column(sa.Integer, primary_key=True)
    patient_id = sa.Column(sa.Integer, sa.ForeignKey("patient.id"), nullable=False)
    interval_id = sa.Column(sa.Integer, sa.ForeignKey("interval.id"), nullable=False)
    source = sa.Column(sa.String, nullable=False)
    filter = sa.Column(sa.String, index=True)
    is_duplicate = sa.Column(sa.Boolean, index=True)
    
class Gene(Base):
    __tablename__ = 'gene'
    
    # id is ensembl id without the ENSG part
    id = sa.Column(sa.Integer, primary_key=True)
    symbol =  sa.Column(sa.String, index=True)
    chrom = sa.Column(sa.String, nullable=False)
    start = sa.Column(sa.Integer, index=True)
    end = sa.Column(sa.Integer, index=True)

class Interval_Gene(Base):
    __tablename__ = 'interval_gene'
    
    id = sa.Column(sa.Integer, primary_key=True)
    interval_id = sa.Column(sa.Integer, sa.ForeignKey("interval.id"), nullable=False)
    gene_id = sa.Column(sa.Integer, sa.ForeignKey("gene.id"), nullable=False)

class HPO(Base):
    __tablename__ = 'HPO'
    
    # id is HPO id without the HP: part
    id = sa.Column(sa.Integer, primary_key=True)
    name = sa.Column(sa.String, index=True)
    definition = sa.Column(sa.String)
    comment = sa.Column(sa.String)

class HPO_HPO(Base):
    __tablename__ = 'HPO_HPO'
    
    id = sa.Column(sa.Integer, primary_key=True)
    hpo_id = sa.Column(sa.Integer, sa.ForeignKey("HPO.id"), nullable=False)
    parent_hpo_id = sa.Column(sa.Integer, sa.ForeignKey("HPO.id"), nullable=False)

class HPO_Gene(Base):
    __tablename__ = 'HPO_gene'
    
    id = sa.Column(sa.Integer, primary_key=True)
    hpo_id = sa.Column(sa.Integer, sa.ForeignKey("HPO.id"), nullable=False)
    gene_id = sa.Column(sa.Integer, sa.ForeignKey("gene.id"), nullable=False)

class Patient_HPO(Base):
    __tablename__ = 'Patient_HPO'
    
    id = sa.Column(sa.Integer, primary_key=True)
    hpo_id = sa.Column(sa.Integer, sa.ForeignKey("HPO.id"), nullable=False)
    patient_id = sa.Column(sa.Integer, sa.ForeignKey("patient.id"), nullable=False)

