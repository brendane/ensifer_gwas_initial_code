#!/usr/bin/env python

import sys

from Bio import SeqIO

max_len = int(sys.argv[1])

for rec in SeqIO.parse(sys.stdin, 'fasta'):
    if len(rec) > max_len:
        continue
    SeqIO.write(rec, sys.stdout, 'fasta')
