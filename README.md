# SVRare Import
This repo is for importing structural variants, together with relevant annotations, into a sqlite database. The size of the resulting sqlite database is around 200M for 100 whole genome sequenced germlines.

## Install
`pip install -r requirements.txt`

## get public SV database
### gnomAD
1. Download the right build from https://gnomad.broadinstitute.org/downloads (search for SV 2.1 sites VCF)
2. Modify the file to expose SV end with following snippet
```python
import gzip

infile = 'gnomad_v2.1_sv.sites.vcf.gz'
outfile = 'gnomad_v2.1_sv.sites.converted.vcf'

with gzip.open(infile, 'rt') as inf, open(outfile, 'wt') as outf:
  for line in inf:
    if line.startswith('##'):
      outf.write(line)
      continue
    row = line.split('\t')
    if line.startswith('#'):
      row = row[:2] + ['END'] + row[2:]
      outf.write('\t'.join(row))
    else:
      end = ''
      for info in row[-1].split(';'):
        if info.startswith('END'):
          end = info.split('=')[1]
          break
      row = row[:2] + [end] + row[2:]
      outf.write('\t'.join(row))
```
3. bgzip and index the output with `bgzip gnomad_v2.1_sv.sites.converted.vcf && tabix -s 1 -b 2 -e 3 gnomad_v2.1_sv.sites.converted.vcf.gz`

## dbVAR
Download non-redundant LOSS and GAIN files from https://github.com/ncbi/dbvar/blob/master/Structural_Variant_Sets/Nonredundant_Structural_Variants/README.md  
Since the gz files are not bgzipped, you'll have to unzip and bgzip the files, then index them with
`tabix -s 1 -b 2 -e 3 xxx.tsv.gz`

## Decipher
Download from https://decipher.sanger.ac.uk/about/downloads/data
Note that the file has some unsorted records, and needs to be sorted. You can use the following snippet
```python
import gzip
infile = 'population_cnv.txt.gz'
outfile = 'population_cnv.sorted.txt'
header = None
data = []
with gzip.open(infile, 'rt') as inf:
  for line in inf:
    if header is None:
      header = line
      continue
    row = line.rstrip().split('\t')
    row[3] = int(row[3])
    row[2] = int(row[2])
    data.append(row)
data.sort(key=lambda x: (x[1], x[2], x[3]))

with open(outfile, 'wt') as outf:
  outf.write(header)
  for d in data:
    outf.write('\t'.join([str(i) for i in d]) + '\n)
```
Then bgzip and index using `bgzip population_cnv.sorted.txt && tabix -s 2 -b 3 -e 4 population_cnv.sorted.txt.gz`

## Prepare database
First of all, edit `config.yml` for locations of the files. Note that the `params` section is not relevant and will be removed soon.

In `patients.tsv` (a tab-delimited file), one can adopt the following format:

| family_id | name | is_proband | relation_to_proband | disease | HPO | bam_path | canvas_path | manta_path | is_solved |
|---------|-------|------|------|------|------|------|------|------|------|
| F_01 | OE_01 | 1 | | Retinal dystrophy | HP:0000556,HP:0000505 | /path/to/bam/OE_01.bqsr.bam | /path/to/canvas.vcf.gz | /path/to/manta.vcf.gz* | 0

\* Note that Manta does not convert breaking points to inversions by default. Please consult Manta's [manual](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#inversions).

Then one can run `python prepare.py` to download relevant data, construct the database.

## Import data to database

Simply run `python import_SV.py`
