#!/usr/bin/env python2
"""
    Transform a VCF file into a FASTA file using just the variants.

    vcf2tsv.py --output <output file> <input file>

    OPTIONS

        --subset        A bed file with positions to include, default is
                        to include everything
        --min-gt        Minimum genotyping rate (0)
        --max-pos       Maximum position in genome (for getting rid of RDVs)
        --exclude       Exclude these replicons

    Ignores heterozygous sites.
"""

#==============================================================================#

import argparse
import itertools
import re

import vcf

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--min-gt', default=0.0, type=float)
parser.add_argument('--subset')
parser.add_argument('--max-pos', default=10E9, type=int)
parser.add_argument('--exclude')
parser.add_argument('input')
args = parser.parse_args()

if args.subset is None:
    positions = None
else:
    positions = set()
    with open(args.subset, 'rb') as handle:
        for line in handle:
            contig, _, pos = line.strip().split()[:3]
            positions.add((contig, pos))

exclude = set()
if args.exclude is not None:
    for r in args.exclude.split(','):
        exclude.add(r.strip())

rdr = vcf.Reader(filename=args.input)
data = [[] for s in rdr.samples]
samples = rdr.samples
het_count = 0
for rec in rdr:
    if (positions is not None) and ((rec.CHROM, str(rec.POS)) not in positions):
        continue
    if rec.CHROM in exclude:
        continue
    if rec.POS > args.max_pos:
        continue
    if max(map(len, rec.alleles)) > 1:
        continue
    gts = []
    missed = 0
    if len(rec.samples) != len(samples):
        print rec.CHROM, rec.POS
        print "Wrong number of samples"
        raise Exception("Wrong number of samples!")
    for s in rec.samples:
        gt = s.gt_alleles
        if gt[0] is None or s.is_het:
            if s.is_het:
                het_count += 1
            g = 'N'
            missed += 1
        else:
            # Only take the first allele because we're skipping
            # heterozygotes
            g = rec.alleles[int(gt[0])]
        gts.append(str(g))
    if 1 - (missed / float(len(gts))) < args.min_gt:
        continue
    for i, g in enumerate(gts):
        data[i].append(g)

with open(args.output, 'wb') as out:
    for sample, gts in itertools.izip(samples, data):
        out.write('>' + sample + '\n')
        out.write(''.join(gts) + '\n')

print 'Found %i heterozygous calls' % het_count
