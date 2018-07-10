#!/usr/bin/env Rscript
#
# Make a spreadsheet with the maximum temperature, NaCl concentration,
# and PEG concentration for each strain.
#
# INPUT
#
# phenotypes: non-biolog traits in the 2016-07-15 directory
#
# This script runs very fast and uses very little memory, so
# no need for qsub.
#
# UPDATE 19 Aug 2016: now makes a spreadsheet of raw values too.
#

library(dplyr)
library(magrittr)

projdir = '/home/tiffinp/epste051/project/metab_gwas'
outdir  = file.path(projdir, 'results', 'process_phenotypes',
                    'composite_2017-01-31')
workdir = file.path(outdir, 'working')
scriptdir = file.path(outdir, 'script_copies')
infile = file.path(projdir, 'data', 'phenotype', '2016-07-15',
                   'non_biolog.tsv')

dir.create(outdir, recursive=TRUE)
dir.create(workdir, recursive=TRUE)
dir.create(scriptdir, recursive=TRUE)

## Copy this file
cargs = commandArgs(trailingOnly=FALSE)
argv0 = gsub('--file=', '', cargs[grepl('--file=', cargs)])
file.copy(argv0,
          file.path(scriptdir,
                    paste0(format(Sys.time(), '%Y-%m-%d-%H%M'),
                           '-', basename(argv0))))

## Read in the data
dat = read.csv(infile, sep='\t', as.is=TRUE, header=TRUE)

## Salt
salt_levels = sort(c(200,300,500,600,800))
max_salt = numeric(nrow(dat))
for(i in 1:nrow(dat)) {
    ms = 0
    for(sl in salt_levels) {
        if(dat[i, paste0('NaCLmM_', sl)] == 1)
            ms = sl
    }
    max_salt[i] = ms
}

## Temp
temp_levels = sort(c(20, 28, 37, 40, 43))
max_temp = numeric(nrow(dat))
for(i in 1:nrow(dat)) {
    ms = 0
    for(tl in temp_levels) {
        if(dat[i, paste0('Temp_', tl)] == 1)
            ms = tl
    }
    max_temp[i] = ms
}

## PEG
peg_levels = sort(c(0.05,0.10,0.15,0.20,0.25))
max_peg = numeric(nrow(dat))
for(i in 1:nrow(dat)) {
    ms = 0
    for(pl in peg_levels) {
        pl_ = gsub('0\\.', '.', formatC(pl, format='f', digits=2))
        if(dat[i, paste0('PEG_', pl_)] == 1)
            ms = pl
    }
    max_peg[i] = ms
}

## Write the output
write.table(data.frame('Strain_ID'=dat[, 'Strain_ID'],
                       'GWAS_ID'=dat[, 'GWAS_ID'],
                       Temp_max=max_temp, Salt_max=max_salt,
                       PEG_max=max_peg),
            file=file.path(outdir, 'pheno.tsv'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
