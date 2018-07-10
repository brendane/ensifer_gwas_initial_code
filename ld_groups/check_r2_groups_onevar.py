#!/usr/bin/env python2.7
"""
    Use plink to check the one_variant.tsv file.

    check_r2_groups_onevar.py <plink bed file> <onevar file> <threshold>

    Optionally, check a random subset.
"""

import os
import random
import subprocess
import sys

threshold = round(float(sys.argv[3]), 3)

groups = {}
with open(sys.argv[2], 'rb') as ihndl:
    for line in ihndl:
        row = line.strip().split('\t')
        groups[row[3]] = row[-1].split(',')

if len(sys.argv) > 4:
    seeds = set(random.sample(groups, int(sys.argv[4])))
else:
    seeds = set(groups)

for seed, group in groups.iteritems():
    if seed not in seeds:
        continue
    for variant in group:
        if variant == seed:
            continue
        cmd = ['plink', '--bfile', sys.argv[1], '--allow-extra-chr',
               '--ld', seed, variant]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        if p.wait() != 0:
            raise Exception('command failed on %s and %s' % (seed, variant))
        result = [l for l in p.communicate()[0].split('\n') if 'R-sq' in l][0]
        if round(float(result.strip().split()[2]), 3) < threshold:
            raise Exception('%s and %s: r2 < %f' %
                            (seed, var, threshold))

with open('__tmp__.seeds', 'wb') as ohnd:
    for seed in seeds:
        ohnd.write(seed + '\n')

cmd = ['plink', '--bfile', sys.argv[1], '--allow-extra-chr',
       '--r2', 'inter-chr', '--ld-window-r2', '0', '--extract',
       '__tmp__.seeds', '--out', '__tmp__']
p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
if p.wait() != 0:
    raise Exception('command failed on seeds')

with open('__tmp__.ld', 'rb') as ihnd:
    ihnd.readline()
    for line in ihnd:
        fields = line.strip().split()
        r2 = float(fields[-1])
        if round(r2, 3) >= threshold:
            print '%s and %s (%f) should be grouped' % (fields[2], fields[5], r2)

os.unlink('__tmp__.ld')
os.unlink('__tmp__.nosex')
os.unlink('__tmp__.log')
os.unlink('__tmp__.seeds')

exit(0)
