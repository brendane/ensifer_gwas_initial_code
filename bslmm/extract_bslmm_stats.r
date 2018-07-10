#!/usr/bin/env Rscript
#
# Extract means and HPD intervals for PGE, PVE, and n_gamma. Works
# on just one file.
#

library(coda)
library(data.table)

argv = commandArgs(trailingOnly=TRUE)

thin = as.numeric(argv[2])
infile = argv[1]
infile_log = argv[3]

make_vec = function() numeric(1)
stats = data.frame('run'=make_vec(),
                   'pve_95_l'=make_vec(),
                   'pve_95_h'=make_vec(),
                   'pve_mean'=make_vec(),
                   'pge_95_l'=make_vec(),
                   'pge_95_h'=make_vec(),
                   'pge_mean'=make_vec(),
                   'ng_95_l'=make_vec(),
                   'ng_95_h'=make_vec(),
                   'ng_mean'=make_vec(),
                   'lmm_pve'=make_vec())

if(grepl('\\.gz$', infile)) {
    x = fread(paste('zcat', infile))
} else {
    x = fread(infile)
}

mcmc_data = as.mcmc(x, start=1, end=thin*nrow(x),
                    thin=thin)
hpd_95 = HPDinterval(mcmc_data, 0.95)
means = colMeans(x)

stats[1, 'pve_95_l'] = hpd_95['pve', 'lower']
stats[1, 'pve_95_h'] = hpd_95['pve', 'upper']
stats[1, 'pge_95_l'] = hpd_95['pge', 'lower']
stats[1, 'pge_95_h'] = hpd_95['pge', 'upper']
stats[1, 'ng_95_l'] = hpd_95['n_gamma', 'lower']
stats[1, 'ng_95_h'] = hpd_95['n_gamma', 'upper']
stats[1, 'pve_mean'] = means['pve']
stats[1, 'pge_mean'] = means['pge']
stats[1, 'ng_mean'] = means['n_gamma']
stats[1, 'run'] = 1

l = scan(infile_log, what='character', sep='\n', quiet=TRUE)
stats[1, 'lmm_pve'] = as.numeric(gsub('.+= ', '', l[grepl('pve', l)][1]))


write.table(stats, file=stdout(), sep='\t', col.names=TRUE,
            row.names=FALSE, quote=FALSE)
