#!/usr/bin/env python2.7
"""
    Make a table of all phenotypes.

    Not run with qsub because that's not necessary.
"""

import csv
import os
import os.path as osp
import re
import shutil
import sys
import time

projdir = os.getenv('HOME') + '/project/metab_gwas'

outdir = projdir + '/results/tabulate_phenotypes/2018-05-03'

biolog_input_file = projdir + '/results/process_phenotypes/2016-08-17_24hr/biolog.tsv'
nb_input_file = projdir + '/data/phenotype/2016-07-15/non_biolog.tsv'
bioclim_input_file = projdir + '/results/process_phenotypes/bioclim_30sec_2016-10-13/pca.tsv'
bioclim_input_file2 = projdir + '/data/strain/2016-07-27/data_bioclim.tsv'
composite_input_file = projdir + '/results/process_phenotypes/composite_2017-01-31/pheno.tsv'
growth_input_file = projdir + '/data/phenotype/growth_rate_20170210/growth_rate_gwasid.tsv'
plant_input_file = projdir + '/data/phenotype/plant_data/2017-06-30/adjusted_2017-07-12.tsv'

mel_153_list = projdir + '/data/strain/lists/153meliloti.txt'


if not os.access(osp.dirname(outdir), os.F_OK):
    os.mkdir(osp.dirname(outdir))
if not os.access(outdir, os.F_OK):
    os.mkdir(outdir)
if not os.access(outdir + '/script_copies', os.F_OK):
    os.mkdir(outdir + '/script_copies')

shutil.copyfile(sys.argv[0], outdir + '/script_copies/' +
                time.strftime('%Y-%m-%d-%H%M') + '-' +
                osp.basename(sys.argv[0]))

mel_153 = set()
with open(mel_153_list, 'rb') as ihnd:
    for line in ihnd:
        mel_153.add(line.strip())

data = {}

biolog_traits = []
with open(biolog_input_file, 'rb') as ihnd:
    rdr = csv.DictReader(ihnd, delimiter='\t')
    first = True
    for row in rdr:
        strain = row['GWAS_ID']
        # Drop uninteresting strains
        if strain not in mel_153:
            s = 'USDA' + strain
            if s not in mel_153:
                continue
            strain = s
        # Drop the old data from Leo:
        if row['dataset'] != '1':
            continue
        data[strain] = {}
        for pheno in rdr.fieldnames[4:]:
            colname = pheno
            pname = re.sub('^X2', '2', pheno)
            if colname not in row:
                raise Exception('could not find %s' % colname)
            data[strain][pname] = row[colname]
            if first:
                biolog_traits.append(pname)
        first = False

with open(bioclim_input_file, 'rb') as ihnd:
    rdr = csv.DictReader(ihnd, delimiter='\t')
    for row in rdr:
        strain = row['GWAS_ID']
        # Drop uninteresting strains
        if strain not in mel_153:
            s = 'USDA' + strain
            if s not in mel_153:
                continue
            strain = s
        if strain not in data:
            data[strain] = {}
        for pheno in rdr.fieldnames[2:]:
            colname = pheno
            if colname not in row:
                raise Exception('could not find %s' % colname)
            data[strain][pheno] = row[colname]

with open(bioclim_input_file2, 'rb') as ihnd:
    rdr = csv.DictReader(ihnd, delimiter='\t')
    for row in rdr:
        strain = row['GWAS_ID']
        # Drop uninteresting strains
        if strain not in mel_153:
            s = 'USDA' + strain
            if s not in mel_153:
                continue
            strain = s
        if strain not in data:
            data[strain] = {}
        bioclim_raw_traits = [x for x in rdr.fieldnames if x.startswith('bio_')]
        for pheno in bioclim_raw_traits:
            colname = pheno
            if colname not in row:
                raise Exception('could not find %s' % colname)
            data[strain][pheno] = row[colname]

with open(growth_input_file, 'rb') as ihnd:
    rdr = csv.DictReader(ihnd, delimiter='\t')
    for row in rdr:
        strain = row['GWAS_ID']
        # Drop uninteresting strains
        if strain not in mel_153:
            s = 'USDA' + strain
            if s not in mel_153:
                continue
            strain = s
        if strain not in data:
            data[strain] = {}
        for pheno in ['Growth_Rate']:
            colname = pheno
            if colname not in row:
                raise Exception('could not find %s' % colname)
            data[strain][pheno] = row[colname]


nb_traits = []
with open(nb_input_file, 'rb') as ihnd:
    rdr = csv.DictReader(ihnd, delimiter='\t')
    first = True
    for row in rdr:
        strain = row['GWAS_ID']
        # Drop uninteresting strains
        if strain not in mel_153:
            s = 'USDA' + strain
            if s not in mel_153:
                continue
            strain = s
        if strain not in data:
            data[strain] = {}
        for pheno in rdr.fieldnames[2:]:
            colname = pheno
            if colname not in row:
                raise Exception('could not find %s' % colname)
            data[strain][pheno] = row[colname]
            if first:
                nb_traits.append(pheno)
        first = False

with open(composite_input_file, 'rb') as ihnd:
    rdr = csv.DictReader(ihnd, delimiter='\t')
    for row in rdr:
        strain = row['GWAS_ID']
        # Drop uninteresting strains
        if strain not in mel_153:
            s = 'USDA' + strain
            if s not in mel_153:
                continue
            strain = s
        if strain not in data:
            data[strain] = {}
        for pheno in rdr.fieldnames[2:]:
            colname = pheno
            if colname not in row:
                raise Exception('could not find %s' % colname)
            data[strain][pheno] = row[colname]


plant_traits_categories = ['nodule', 'weight']
plant_traits = ['nodule_A17', 'weight_A17', 'nodule_R108', 'weight_R108']
with open(plant_input_file, 'rb') as ihnd:
    rdr = csv.DictReader(ihnd, delimiter='\t')
    for row in rdr:
        strain = row['strain']
        # Drop uninteresting strains
        if strain not in mel_153:
            s = 'USDA' + strain
            if s not in mel_153:
                continue
            strain = s
        if strain not in data:
            data[strain] = {}
        for pheno in plant_traits_categories:
            colname = pheno
            if colname not in row:
                raise Exception('could not find %s' % colname)
            data[strain][pheno + '_' + row['plant_genotype']] = row[colname]

traits = sorted(data[list(mel_153)[0]].keys())
with open(outdir + '/phenotypes.mel153.tsv', 'wb') as ohnd:
    ohnd.write('\t'.join(['strain'] + traits) + '\n')
    for strain in sorted(mel_153):
        if strain not in data:
            continue
        ohnd.write(strain)
        for trait in traits:
            if trait in data[strain]:
                ohnd.write('\t' + data[strain][trait])
            else:
                ohnd.write('\tNA')
        ohnd.write('\n')

with open(outdir + '/nb_biolog_traits.txt', 'wb') as ohnd:
    for t in biolog_traits:
        ohnd.write(t + '\tbiolog\n')
    for t in nb_traits:
        ohnd.write(t + '\tbinary\n')
