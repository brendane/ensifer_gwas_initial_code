#!/usr/bin/env python2
"""
    Process the 24 hr OD540's from the biolog plates.

    Same as 2016-07-29_24hr, but it includes five extra strains.

    INPUT

    1. Absorbance data sent by Mike on 2016-06-14, and then processed by
       me (I subtracted the water absorbance).
    2. Leo Cartegena's biolog data
    3. Data processed by Liana to get strain name information

    PROCESSING

    The input data have already had the water well absorbance subtracted.
    This script takes the mean of alpha_D_Glucose, D_Fructose, and
    Sucrose for each strain and divides all the well by that mean. If
    there is a negative absorbance for any of the three compounds,
    the strain is thrown out (only one strain has this issue, and the
    absorbance is negative for all three compounds). Negative absorbances
    in other substrates are changed to zero.

    Two strains are removed because there were experimental problems in
    the most recent experiment:
        USDA1002 ("1002" in this dataset)
        USDA207  ("207")

    There is no need to run this script with qsub.
"""

remove = {'207', '1002'}

import copy
import csv
import os
import os.path as osp
import shutil
import sys
import time

projdir = '/home/tiffinp/epstein1/project/metab_gwas/'
biologfile = osp.join(projdir, 'data', 'phenotype', 'biolog_2016-06-14',
                      '24hrs', 'data.tsv')
oldblgfile = osp.join(projdir, 'data', 'phenotype', 'old_biolog',
                      'processed', 'data_24hr.tsv')
newestfile = osp.join(projdir, 'data', 'phenotype', '2016-08-15',
                      '24hrs', 'data.tsv')
lianafile = osp.join(projdir, 'data', 'phenotype', '2016-07-15',
                     'biolog.tsv')
outdir = osp.join(projdir, 'results', 'process_phenotypes', '2016-08-17_24hr')
scriptdir = osp.join(outdir, 'script_copies')

# Create the output directories
try:
    os.mkdir(osp.dirname(outdir))
except:
    pass
try:
    os.mkdir(outdir)
    os.mkdir(scriptdir)
except:
    pass

# Make a time-stamped copy of this script
shutil.copyfile(sys.argv[0],
                osp.join(scriptdir,
                         time.strftime('%Y-%m-%d-%H%M') + '-' + osp.basename(sys.argv[0])))


# Process the data
str_id = {}
with open(lianafile, 'rb') as handle:
    rdr = csv.DictReader(handle, delimiter='\t')
    for row in rdr:
        str_id[row['GWAS_ID']] = row['Strain_ID']

header = False
with open(osp.join(outdir, 'biolog.tsv'), 'wb') as out:
    for blgfile in [oldblgfile, biologfile, newestfile]:
        with open(blgfile, 'rb') as handle:
            rdr = csv.DictReader(handle, delimiter='\t')
            colnames = copy.deepcopy(rdr.fieldnames)
            colnames.insert(1, 'GWAS_ID')
            colnames.insert(2, 'dataset')
            colnames[0] = 'Strain_ID'
            if not header:
                out.write('\t'.join(colnames) + '\n')
                header = True
            for row in rdr:
                if row['strain'] in remove and blgfile == biologfile:
                    continue
                glucose = float(row['alpha_D_Glucose'])
                sucrose = float(row['Sucrose'])
                fructose = float(row['D_Fructose'])
                if glucose < 0 or sucrose < 0 or fructose < 0:
                    continue
                std = (glucose + fructose + sucrose) / 3.
                for c in colnames:
                    if c == 'Strain_ID':
                        try:
                            out.write(str_id[row['strain']])
                        except KeyError:
                            if row['strain'] in {'1467', '1566', '1593'}:
                                out.write('USDA' + row['strain'])
                            else:
                                out.write(row['strain'])
                        continue
                    if c == 'GWAS_ID':
                        out.write('\t' + row['strain'])
                        continue
                    if c == 'dataset':
                        if blgfile == oldblgfile:
                            out.write('\t0')
                        else:
                            out.write('\t1')
                        continue
                    od = float(row[c]) / std
                    if od < 0: od = 0.
                    out.write('\t' + str(od))
                out.write('\n')
