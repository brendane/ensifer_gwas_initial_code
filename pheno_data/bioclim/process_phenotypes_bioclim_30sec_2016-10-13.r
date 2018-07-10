#!/usr/bin/env Rscript
#
# Make a spreadsheet of PCs of Bioclim variables for all strains.
# Also make a spreadsheet of raw values.
#
# INPUT
#
# bioclim: 30 second resolution data downloaded August 2016
#
# This script runs very fast and uses very little memory, so
# no need for qsub.
#
# UPDATE 19 Aug 2016: now makes a spreadsheet of raw values too.
#
# UPDATE 06 Mar 2017: now output PC loadings.
#

library(dplyr)
library(magrittr)

projdir = '/home/tiffinp/epste051/project/metab_gwas'
outdir  = file.path(projdir, 'results', 'process_phenotypes',
                    'bioclim_30sec_2016-10-13')
workdir = file.path(outdir, 'working')
scriptdir = file.path(outdir, 'script_copies')
infile = file.path(projdir, 'data', 'strain', '2016-07-27',
                   'data_bioclim.tsv')

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

## Read in the file with the bioclim data and take only the
## the individuals that are in our current genotyping set.
## Also trim down to just the bioclim columns
dat = read.csv('/home/tiffinp/epste051/project/metab_gwas/data/strain/2016-07-27/data_bioclim.tsv',
               sep='\t', header=TRUE, as.is=TRUE) %>%
    filter(!is.na(bio_1))
pcadat = dat %>%
    select(starts_with('bio_')) %>%
    as.matrix()

## Run a PCA using the standard R functions
## As recommended in the manual, I use prcomp instead of princomp,
## and I center and scale the matrix.
pca = prcomp(pcadat, center=TRUE, scale=TRUE)

pc = pca$x
colnames(pc) = dimnames(pca$x)[[2]]

## Write the output
write.table(data.frame('Strain_ID'=dat[, 'Strain_ID'],
                       'GWAS_ID'=dat[, 'GWAS_ID'],
                       pc),
            file=file.path(outdir, 'pca.tsv'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

write.table(dat %>% select(GWAS_ID, Strain_ID, starts_with('bio_')),
            file=file.path(outdir, 'bioclim.tsv'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

write.table(data.frame('variable'=rownames(pca$rotation), pca$rotation),
            file=file.path(outdir, 'loadings.tsv'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

cat(pca$sdev^2/sum(pca$sdev^2), '\n',
    file=file.path(outdir, 'variance.txt'))
