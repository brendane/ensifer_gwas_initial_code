#!/usr/bin/env python2.7
"""
    Given a clstr file from CD-Hit and a directory with ipscan
    annotations (named <strain>_...), annotate each cluster.

    annotate_cdhit_from_ipscan.py --output <output file> 
    --ref-strain <strain> --ref-file <gff3 file for reference>
    <cluster file> <input directory>

    A number of features of this script are probably specific to my
    current dataset.

    Prefers PFAM and then TIGRFAM.
"""

import argparse
import csv
import gzip
import itertools
import re
import os

def get_tags(gffstring):
    ret = {}
    for x in gffstring.split(';'):
        k, v = x.split('=')
        ret[k] = v
    return ret

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--ref-strain')
parser.add_argument('--ref-file')
parser.add_argument('clustr')
parser.add_argument('annotdir')
args = parser.parse_args()

strains = set()
cluster_number = -1
clusters = {}
ofun = open
if args.clustr.endswith('.gz'):
    ofun = gzip.open
with ofun(args.clustr, 'rb') as ih:
    for line in ih:
        line = line.strip()
        if line.startswith('>'):
            cluster_number = int(re.sub('>Cluster ', '', line))
        else:
            gene = re.sub('>', '', line.split(' ')[1])
            gene = re.sub('\.\.\.$', '', gene)
            strain = gene[1:].split('.')[0]
            strains.add(strain)
            clusters[gene] = cluster_number

ref_annot = {}
ref_strain = args.ref_strain
if args.ref_file is not None:
    with open(args.ref_file, 'rb') as ih:
        rdr = csv.reader(ih, delimiter='\t')
        for row in rdr:
            tags = get_tags(row[8])
            a = ''
            try:
                g = tags['gene']
                a += g + ': '
            except KeyError:
                pass
            try:
                p = tags['product']
                a += p + ' '
            except KeyError:
                pass
            i = tags['ID']
            l = tags['locus_tag']
            a += '(ref_id=' + i + ')'
            try:
                c = clusters[ref_strain + '.' + l]
            except KeyError:
                ## This gene not included - maybe too long
                continue
            if c not in ref_annot:
                ref_annot[c] = []
            ref_annot[c].append(a)

denovo_annot = {}
for fname in os.listdir(args.annotdir):
    strain = re.sub('_.+', '', fname)
    with open(args.annotdir + '/' + fname, 'rb') as ih:
        rdr = csv.reader(ih, delimiter='\t')
        for lt, rows in itertools.groupby(rdr, lambda x: x[0]):
            l = re.sub('_.+', '', lt)
            try:
                c = clusters[strain + '.' + l]
            except KeyError:
                ## This gene not included - maybe too long
                continue
            if c not in denovo_annot:
                denovo_annot[c] = []
            for row in rows:
                annot_type = row[3]
                if annot_type == 'Pfam' or annot_type == 'TIGRFAM':
                    annot = annot_type + ': ' + row[5]
                    if annot not in denovo_annot[c]:
                        denovo_annot[c].append(annot)

with open(args.output, 'wb') as oh:
    for c in xrange(cluster_number + 1):
        a = ''
        if c in ref_annot:
            a += 'REF=' + ',,'.join(ref_annot[c]) + '; '
            if len(ref_annot[c]) > 1:
                print c
        if c in denovo_annot:
            a += 'DENOVO=' + ',,'.join(denovo_annot[c]) + '; '
        if len(a) == 0:
            a = 'no annotation found'
        oh.write('denovo\t.\tCDS\t' + str(c + 1) + '\t' +
                 str(c + 1) + '\t.\t.\t0\t' + 'ID=cluster-' + str(c) +
                 '; ' + a + '\n')
