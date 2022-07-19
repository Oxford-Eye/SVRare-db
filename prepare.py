'''
Prepare for action
0. remove db
1. Download HPO terms and HPO genes
2. import HPO and genes
'''
import os
import sys
import yaml
import gzip
import urllib.request
from sqlalchemy import create_engine
from sqlalchemy.orm import Session, declarative_base
from lib import models

def main(config):
    
    engine = create_engine(config['db'])
    models.Base.metadata.drop_all(engine)
    models.Base.metadata.create_all(engine)
    
    # HPO
    with urllib.request.urlopen(config['hpo_obo_url']) as f:
        html = f.read().decode('utf-8')
        entities = [] # HPO models
        entries = [] # HPO dict, which has more useful fields (is_a)
        for block in html.split('\n\n'):
            if not block.startswith('[Term]\n'):
                continue
            
            entity = {
                'is_a': []
            }
            for line in block.split('\n'):
                
                if line == '[Term]':
                    continue
                key, val = line.split(': ', 1)
                if key == 'is_a':
                    entity['is_a'].append(int(val.split(' ')[0].lstrip('HP:')))
                elif key == 'id':
                    entity['id'] = int(val.lstrip('HP:'))
                else:
                    entity[key] = val
            entries.append(entity)        
            entities.append(
                models.HPO(
                    id = entity['id'],
                    name = entity['name'],
                    definition = entity.get('def', None),
                    comment = entity.get('comment', None)
                )
            )
        with Session(engine) as session:
            session.add_all(entities)
            session.commit()
    
        # HPO is_a
        entities = []
        for entry in entries:
            for parent in entry['is_a']:
                entities.append(
                    models.HPO_HPO(
                        hpo_id = entry['id'],
                        parent_hpo_id = parent
                    )
                )
        with Session(engine) as session:
            session.add_all(entities)
            session.commit()
    
    # gene
    # get gnomad gene constraints
    gnomad_constraints = {}
    with urllib.request.urlopen(config['gnomad_constraint_url']) as zfd:
        fd = gzip.GzipFile(fileobj=zfd, mode="r")
        header = []
        for line in fd:
            line = line.decode('utf8')
            row = line.rstrip().split('\t')
            if not header:
                header = row
                continue
            row_dict = dict(zip(header, row))
            gene_id = int(row_dict['gene_id'].lstrip('ENSG'))
            gnomad_constraints[gene_id] = {
                'pli': row_dict['pLI'],
                'prec': row_dict['pRec'],
                'oe_lof_upper': row_dict['oe_lof_upper'],
            }
    # store genes[symbol] = id for HPO_gene,  which doesn't have ensembl id
    genes = {}
    with gzip.open(config['gene_tbx'], 'rt') as inf:
        entities = []
        for line in inf:
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            if row[2] != 'gene':
                continue
            info = {}
            for info_field in row[-1].split('; '):
                key, val = info_field.split(' ')
                val = val.strip('"')
                info[key] = val
            
            gene_id = int(info['gene_id'].lstrip('ENSG'))
            # sometimes in GRCh38 there's no gene_name
            symbol = info.get('gene_name', info['gene_id'])
            genes[symbol] = gene_id
            constraints = gnomad_constraints.get(gene_id, {
                'pli': None,
                'prec': None,
                'oe_lof_upper': None
            })
            
            entities.append(models.Gene(
                id = gene_id,
                symbol = symbol,
                chrom = row[0],
                start = int(row[3]),
                end = int(row[4]),
                pli = constraints['pli'],
                prec = constraints['prec'],
                oe_lof_upper =  constraints['oe_lof_upper'],
            ))
        with Session(engine) as session:
            session.add_all(entities)
            session.commit()
    # HPO_gene
    with urllib.request.urlopen(config['hpo_gene_url']) as f:
        entities = []
        html = f.read().decode('utf-8')
        for line in html.split('\n'):
            if line.startswith('#') or not line.strip():
                continue
            row = line.split('\t')
            try:
                hpo_id = int(row[0].lstrip('HP:'))
            except ValueError:
                print(line)
                raise
            symbol = row[3]
            if symbol not in genes:
                continue
            entities.append(models.HPO_Gene(
                hpo_id = hpo_id,
                gene_id = genes[symbol],
            ))
            
        with Session(engine) as session:
            session.add_all(entities)
            session.commit()
            
if __name__ == '__main__':
    with open('config.yml', 'rt') as inf:
        config = yaml.safe_load(inf)
        sys.exit(main(config))
