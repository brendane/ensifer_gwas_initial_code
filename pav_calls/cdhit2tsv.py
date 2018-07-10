#!/usr/bin/env python2.7
"""
    Turn a clstr file from CD-Hit into a tsv file with copy number.

    cdhit2vcf.py <clstr file> <output file> <contig name>

    Assumes genes are named <strain>.<gene name>
"""

import collections
import gzip
import re
import sys

infile = sys.argv[1]
outfile = sys.argv[2]
contig = sys.argv[3]

def ddint():
    return collections.defaultdict(int)

strains = set()
clusters = collections.defaultdict(ddint)

cluster_number = -1
ofun = open
if infile.endswith('.gz'):
    ofun = gzip.open
with ofun(infile, 'rb') as ih:
    for line in ih:
        line = line.strip()
        if line.startswith('>'):
            cluster_number = int(re.sub('>Cluster ', '', line))
        else:
            gene = line.split(' ')[1]
            strain = gene[1:].split('.')[0]
            strains.add(strain)
            clusters[cluster_number][strain] += 1

strs = sorted(strains)

with open(outfile, 'wb') as oh:
    oh.write('chr\tpos\tcluster\t' + '\t'.join(strs) + '\n')
    cls = sorted(clusters)
    for c in cls:
        oh.write(contig + '\t' + str(c + 1) + '\t' + str(c))
        for s in strs:
            oh.write('\t' + str(clusters[c][s]))
        oh.write('\n')
