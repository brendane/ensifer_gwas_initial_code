#!/usr/bin/env Rscript
#
# Combine all the randomly permuted phenotypes made for the 2017-07-17_random
# GEMMA LMM run into tables.
#
# UPDATE 2018-05-21: now runs bio_1 and bio_12
#

n_random = 100

cargs = commandArgs(trailingOnly=FALSE)
script_name = gsub('--file=', '', cargs[grepl('--file=', cargs)])

projdir = '/home/tiffinp/epste051/project/metab_gwas'
outdir = file.path(projdir, 'results', 'random_phenotype_tables', '2017-09-21')
indir = file.path(projdir, 'results/gwas/phenotype/gemma/2017-07-17_random')

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
file.copy(script_name, file.path(outdir,
                                 paste0(format(Sys.time(), '%Y-%m-%d-%H%M'),
                                        '-', basename(script_name))))

traits = grep('K|work|array|script|log',
              list.dirs(indir, full.names=FALSE, recursive=FALSE),
              value=TRUE, invert=TRUE)

for(trait in traits) {
    data_list = vector('list', length=n_random)
    for(i in 0:(n_random-1)) {
        data_list[[i+1]] = read.table(file.path(indir, trait,
                                                paste0('phenotype.', i, '.tsv')),
                                      header=FALSE, as.is=TRUE)[, 2:3]
        colnames(data_list[[i+1]]) = c('strain', i)
    }
    merged = Reduce(function(x, y) merge(x, y, by='strain', all=TRUE),
                    data_list[-1], data_list[[1]])
    write.table(merged, file=file.path(outdir, paste0(trait, '.tsv')),
                sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
}
