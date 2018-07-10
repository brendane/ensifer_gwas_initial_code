#!/usr/bin/env python
"""
    Process a bed file from bedtools closest -D to get the nearest
    downstream feature is within a certain distance and isn't "blocked"
    by another feature.

    This version ("2") anticipates an extra column in the bed file.
"""

import csv
import itertools
import sys

max_dist = int(sys.argv[1])
multiple_genes = False
if len(sys.argv) > 2:
    if sys.argv[2] == '1':
        multiple_genes = True

def get_closest(features, max_dist=10000000, multiple_genes=False):
    features = [f for f in features if abs(int(f[14])) <= max_dist]
    features.sort(key=lambda x: abs(int(x[14])))

    if len(features) == 0:
        return None

    if int(features[0][14]) == 0:
        ret = []
        # If the variant is in the middle of a gene, then
        # we are all set
        if multiple_genes:
            for f in features:
                if int(f[14]) == 0:
                    ret.append(f)
        else:
            ret = [features[0]]
        return ret
    
    # First get the closest downstream
    closest_down = None
    for f in features:
        if int(f[14]) < 0:
            closest_down = f
            break
    # If there isn't any, then there is no match
    if closest_down is None:
        return None

    # If there is only one match, we are done
    if len(features) == 1:
        return [closest_down]

    # Check if there are any upstream matches between
    # closest_down and the variant
    up = False
    cd_start, cd_end = int(closest_down[8]), int(closest_down[9])
    a_start, a_end = int(closest_down[1]), int(closest_down[2])
    for f in features:
        if int(f[14]) < 0:
            continue
        start, end = int(f[8]), int(f[9])
        if end >= cd_end and end <= a_start:
            up = True
            break
    if up:
        return None
    else:
        return [closest_down]

rdr = csv.reader(sys.stdin, delimiter='\t')
for coords, _features in itertools.groupby(rdr, lambda x: (x[0], x[1], x[2], x[3])):
    features = list(_features)

    closest_down = get_closest(features, max_dist, multiple_genes)
    if closest_down is None:
        sys.stdout.write('\t'.join(features[0][:5] +
                                   ['.', '.','.', '.', '.', '.', '.', '.', '.'])
                         + '\n')
    else:
        for cd in closest_down:
            sys.stdout.write('\t'.join(cd[:-1]) + '\n')
