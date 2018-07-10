#!/usr/bin/env Rscript
#
# This makes the table of candidate genes for the LD groups that
# explained more variation than the null distribution. Once a variant
# is encountered that explains less than the null, this script stops
# looking for genes. (Less than the null means <= 95% of the random
# permutations.)
#
# UPDATE 12 June 2018: added variant type annotations, but may not
# use them for the ms.
#

library(data.table)
library(dplyr)

n_random = 100
n = 25

projdir = file.path(Sys.getenv('HOME'), 'project', 'metab_gwas')
indir_modsel = file.path(projdir, 'results/pve/phenotype/lm/2017-07-20_forward')
indir_gwa = file.path(projdir, 'results/gwas/phenotype/gemma/2017-07-17')
ld_file = file.path(projdir, 'results/ld/r2_groups/2017-07-04/0.95/output.tsv')
ann_file = file.path(projdir, 'results/annotate_variants/snpeff/2017-04-26/output.rs.one.tsv')

traits = c('2_Aminoethanol', 'Formic_Acid', 'Anti_Gent', 'Anti_Spec',
           'Anti_Strept', 'PEG_max', 'weight_A17', 'weight_R108', 
           'bio_1', 'bio_12')
trait_long_names = c('2_Aminoethanol'='2-Aminoethanol', 
                     'D_Trehalose'='D-Trehalose',
                     'Formic_Acid'='Formic Acid',
                     'L_Fucose'='L-Fucose',
                     'N_Acetyl_D_Glucosamine'='N-Acetyl-D-\nGlucosamine',
                     'Putrescine'='Putrescine',
                     'Anti_Gent'='Gentamicin\nResistance',
                     'Anti_Spec'='Spectinomycin\nResistance',
                     'Anti_Strept'='Streptomycin\nResistance',
                     'metal_Cd'='Cadmium\nTolerance',
                     'PEG_max'='PEG Tolerance',
                     'Temp_max'='Max. Growth Temp.',
                     'Salt_max'='NaCl Tolerance',
                     'weight_A17'='A17 Biomass',
                     'weight_R108'='R108 Biomass',
                     'nodule_A17'='A17 Nodule\nNumber',
                     'nodule_R108'='R108 Nodule\nNumber',
                     'Growth_Rate'='Growth Rate',
                     'bio_1'='Annual Mean Temp.',
                     'bio_12'='Annual Precip.')

## Read variant annotations
ann = read.csv(ann_file, sep='\t', header=FALSE, as.is=TRUE)

## Read LD grouping information
ld = read.csv(ld_file, sep='\t', header=TRUE, as.is=TRUE)

cat('trait\tchrom\tpos\trs\tld_group\tmodsel_rank\tgwa_rank\tgwa_pval\t',
    'model_r2\tmodel_resid_r2\trandom_resid_r2_u95\tseed\tID\ttype\tgenes\n',
    sep='')

for(trait in traits) {
    ## Linear model selection
    dat_modsel = read.csv(file.path(indir_modsel,
                                    trait, 'output.stats.tsv'),
                          header=TRUE, sep='\t', as.is=TRUE)
    ## GWAS top variants
    dat_gwa = read.csv(file.path(indir_gwa,
                                 trait, 'output.0.closest.genes.tsv'),
                       header=FALSE, sep='\t', as.is=TRUE, nrows=500,
                       col.names=c('chrom', 'pos0', 'pos', 'group', 'rank', 'chrom2',
                                   'source', 'type', 'pos3', 'pos4', 'blank',
                                   'strand', 'blank2', 'gene'))
    if(max(dat_gwa[, 'rank']) <= 24) {
        stop('Did not read enough lines')
    }

    dat_gwa_tests = fread(paste('zcat',
                                file.path(indir_gwa, trait, 'output',
                                          'output0.assoc.txt.gz')))

    ## Extract residual R^2 for the random permtuations for 1 - n
    ## top variants.
    randoms = matrix(ncol=n, nrow=n_random)
    for(i in 0:(n_random-1)) {
        run = paste0('random-', i)
        randoms[i+1, ] = dat_modsel[dat_modsel$run == run, 'r2_resid'][1:n]
    }

    ## Combine LD grouping, annotation, and model selection
    ## and print out table rows
    gwa_modsel_ranks = as.numeric(unlist(strsplit(dat_modsel[n, 'ranks_final'], ','))) - 1
    gwa_modsel_rs = unlist(strsplit(dat_modsel[n, 'vars_final'], ','))
    for(i in 1:n) {
        r2 = dat_modsel[dat_modsel$run == 'real', 'r2'][i]
        resid_r2 = dat_modsel[dat_modsel$run == 'real', 'r2_resid'][i]
        random_resid_r2_u95 = quantile(randoms[, i], 0.95)

        if(resid_r2 <= random_resid_r2_u95) {
            break
        }

        gwa_rows = dat_gwa[dat_gwa[, 'rank'] == gwa_modsel_ranks[i], ]
        ld_group = as.numeric(gsub('group-', '', gwa_rows[1, 'group']))

        done = character(nrow(gwa_rows))
        for(j in 1:nrow(gwa_rows)) {
            pos = gwa_rows[j, 'pos']
            chrom = gwa_rows[j, 'chrom']
            rs_s = ld[ld[, 'group'] == ld_group & ld[, 'chrom'] == chrom &
                      ld[, 'pos'] == pos, 'rs']
            if(length(rs_s) > 1) {
                cat('Warning found multiple variants at', chrom, pos,
                    'in', trait, '\n', file=stderr())
            }
            ## This loop is meant to deal with multi-allelic sites; don't know if
            ## it actually works; but the above warning should allow me to
            ## know if it was needed
            for(k in 1:length(rs_s)) {
                rs = rs_s[k]
                if(rs %in% done) {
                    next
                }
                done[which(done == '')[1]] = rs
                seed = ld[ld[, 'rs'] == rs, 'seed']
                modsel_rank = i
                gwa_rank = gwa_modsel_ranks[i]
                genes = paste(gwa_rows[gwa_rows$chrom == chrom & gwa_rows$pos == pos, 'gene'],
                              collapse='; ')
                genes_ids = gsub('ID=(.+?) :.+', '\\1', genes)
                genes = gsub('ID=.+? :', '', genes)
                gwa_pval = dat_gwa_tests[['p_lrt']][dat_gwa_tests[['rs']] == paste0('group-', ld_group)]

                var_type = ann[ann[, 5] == rs, 4]

                cat(trait, '\t', chrom, '\t', pos, '\t', rs, '\t', ld_group,
                    '\t', modsel_rank, '\t',
                    gwa_rank+1, '\t', gwa_pval, '\t', r2, '\t', resid_r2, '\t',
                    random_resid_r2_u95, '\t', seed, '\t', genes_ids, '\t', var_type, '\t',
                    genes, '\n', sep='')
            }
        }
    }
}
