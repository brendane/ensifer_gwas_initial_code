#!/usr/bin/env python2.7
"""
    Turn a clstr file from CD-Hit into a VCF file.

    cdhit2vcf.py <clstr file> <output file> <contig name>

    Assumes genes are named <strain>.<gene name>
"""

import collections
import re
import sys

infile = sys.argv[1]
outfile = sys.argv[2]
contig = sys.argv[3]

strains = set()
clusters = collections.defaultdict(set)

cluster_number = -1
with open(infile, 'rb') as ih:
    for line in ih:
        line = line.strip()
        if line.startswith('>'):
            cluster_number = int(re.sub('>Cluster ', '', line))
        else:
            gene = line.split(' ')[1]
            strain = gene[1:].split('.')[0]
            strains.add(strain)
            clusters[cluster_number].add(strain)

strs = sorted(strains)

with open(outfile, 'wb') as oh:
    oh.write('##fileformat=VCFv4.1\n')
    oh.write('##contig=<ID=' + contig + ',length=' +
             str(len(clusters)) + '>\n')
    oh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    oh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')
    for s in strs:
        oh.write('\t' + s)
    oh.write('\n')
    cls = sorted(clusters)
    for c in cls:
        present = clusters[c]
        oh.write(contig + '\t' + str(c + 1) + '\t' +
                 'rdv-cluster-' + str(c) + '-' + contig + '-' 
                 + str(c+1) + '-A\t' + 'T\tA\t0\t.\t.\tGT')
        for s in strs:
            if s in present:
                oh.write('\t0')
            else:
                oh.write('\t1')
        oh.write('\n')
