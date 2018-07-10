#!/usr/bin/env python2.7
"""
    Merge GEMMA BSLMM runs.

    UPDATE 23 May 2017: calculates beta_hat and beta_bar. Beta_hat is the
    mean beta estimate given gamma == 1, while beta_bar is just the mean
    beta across all steps (beta == 0 when gamma == 0). Should work as long
    as the all the chains are the same length.
"""

import csv
import gzip
import os
import os.path as osp
import re
import sys

iodir = sys.argv[1]

hyp_files = []
param_files = []
for fname in os.listdir(iodir):
    if fname.endswith('param.txt') or fname.endswith('param.txt.gz'):
        param_files.append(fname)
    if fname.endswith('.hyp.txt') or fname.endswith('hyp.txt.gz'):
        hyp_files.append(fname)

pdata = {}
rs = []
for i, fname in enumerate(param_files):
    run = re.sub('\\..+', '', re.sub('output_', '', fname))
    ofun = open
    if fname.endswith('.gz'):
        ofun = gzip.open
    with ofun(osp.join(iodir, fname), 'rb') as ihandle:
        rdr = csv.DictReader(ihandle, delimiter='\t')
        for row in rdr:
            if row['rs'] not in pdata:
                pdata[row['rs']] = []
            pdata[row['rs']].append(row)
            if i == 0:
                rs.append(row['rs'])

with open(osp.join(iodir, 'output.param.merged.tsv'), 'wb') as out:
    out.write('chr\trs\tps\tn_miss\talpha\tbeta_bar\tbeta_hat\tgamma\talphas\tbetas\tgammas\n')
    for r in rs:
        r_data = pdata[r]
        alphas = [float(d['alpha']) for d in r_data]
        betas = [float(d['beta']) for d in r_data]
        gammas = [float(d['gamma']) for d in r_data]
        mean_alpha = sum(alphas) / len(alphas)
        if sum(gammas) > 0:
            beta_bar = sum(b*g for b, g in zip(betas, gammas)) / len(betas)
            beta_hat = sum(b*g for b, g in zip(betas, gammas)) / sum(gammas)
        else:
            beta_bar = 0.
            beta_hat = 0.
        mean_gamma = sum(gammas) / len(gammas)
        out.write(r_data[0]['chr'] + '\t' +
                  r + '\t' +
                  r_data[0]['ps'] + '\t' +
                  r_data[0]['n_miss'] + '\t' +
                  str(mean_alpha) + '\t' +
                  str(beta_bar) + '\t' +
                  str(beta_hat) + '\t' +
                  str(mean_gamma) + '\t' +
                  ','.join(str(f) for f in alphas) + '\t' +
                  ','.join(str(f) for f in betas) + '\t' +
                  ','.join(str(f) for f in gammas) +
                  '\n')

with open(osp.join(iodir, 'output.hyp.merged.tsv'), 'wb') as out:
    for i, fname in enumerate(hyp_files):
        ofun = open
        if fname.endswith('.gz'):
            ofun = gzip.open
        with ofun(osp.join(iodir, fname), 'rb') as ihandle:
            line = ihandle.readline()
            if i == 0:
                out.write('run\t' + line.strip() + '\n')
            for line in ihandle:
                out.write(str(i) + '\t' + line.strip() + '\n')
