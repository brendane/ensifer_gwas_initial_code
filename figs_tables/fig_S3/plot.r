#!/usr/bin/env Rscript
#
# Make the figure showing the proportion of SNPs and PAVs among the
# top 10 GWA LMM hits. Also show mean MAF.
#
# Make a > 15 column for LD group size.
#

library(data.table)
library(dplyr)

make_plot = function(trait, ld, lmm_data, lmm_data_random) {

    n_random = length(lmm_data_random[[trait]])

    ## Get the proportion and size for the top 10 LD groups
    ## real_distribution is ordered by the model selection order, random
    ## distribution is not
    real_distribution = lmm_data[[trait]] %>% left_join(rename(ld, rs=rs_seed))
    random_distribution = lapply(as.character(1:n_random), function(i)
                                 ld[ld[, 'rs_seed'] %in% lmm_data_random[[trait]][[i]][['rs']], ])

    ## And MAF
    real_maf = lmm_data[[trait]][['af']]
    random_maf = lapply(as.character(1:n_random), function(i)
                        lmm_data_random[[trait]][[i]][['af']])

    ## Make some plots
    n_snps = unname(real_distribution[, 'size'] * real_distribution[, 'prop_snp'])
    n_pavs = unname(real_distribution[, 'size'] * (1-real_distribution[, 'prop_snp']))
    mean_prop_random = sapply(random_distribution, function(x) mean(x[, 'prop_snp']))
    mean_size_random = sapply(random_distribution, function(x) mean(x[, 'size']))
    mean_maf_random = sapply(random_maf, function(x) mean(x))

    barplot(rbind(n_snps, n_pavs),
            xlab='Top LD groups from LMM GWA',
            ylab='Number of variants',
            col=c(rgb(142/255, 102/255, 62/255),  'blue'),
            border=NA,
            main=trait_long_names[trait],
            ylim=c(0, 55),
            names.arg=1:length(n_snps),
            las=1)

    hist(mean_prop_random, col='gray60', border='gray80', yaxs='i',
         xlab='Mean proportion of SNPs', ylab='Number of random runs',
         breaks=seq(0, 1, 0.1), main='', xlim=c(0, 1), ylim=c(0, 25),
         las=1)
    abline(v=mean(real_distribution[, 'prop_snp']), lwd=2)
    text(mean(real_distribution[, 'prop_snp']),
         par('usr')[4]*0.9,
         paste0(sum(mean_prop_random > mean(real_distribution[, 'prop_snp'])), '%'),
         pos=4, xpd=NA)

    ## If > 15 per LD group, lump into 16
    mean_size_random[mean_size_random > 15] = 16
    hist(mean_size_random, col='gray60', border='gray80', yaxs='i',
         xlab='Mean LD group size', ylab='Number of random runs', main='',
         breaks=seq(0, 16), xlim=c(0, 17), ylim=c(0, 80), las=1,
         xaxt='n')
    axis(side=1, at=c(0, 4.5, 9.5, 15.5),
         labels=c('0', '5', '10', '>15'))
    abline(v=mean(real_distribution[, 'size']), lwd=2)
    text(mean(real_distribution[, 'size']),
         par('usr')[4]*0.9,
         paste0(sum(mean_size_random > mean(real_distribution[, 'size'])), '%'),
         pos=4, xpd=NA)

    hist(mean_maf_random, col='gray60', border='gray80', yaxs='i',
         xlab='Mean MAF', ylab='Number of random runs', main='',
         breaks=seq(0, 0.5, 0.025),
         xlim=c(0, 0.55), ylim=c(0, 80),
         las=1)
    abline(v=mean(real_maf), lwd=2)
    text(mean(real_maf),
         par('usr')[4]*0.9,
         paste0(sum(mean_maf_random > mean(real_maf)), '%'),
         pos=4, xpd=NA)
}


n_random = 100

projdir = file.path(Sys.getenv('HOME'), 'project', 'metab_gwas')
ld_file = file.path(projdir, 'results/ld/r2_groups/2017-07-04/0.95/one_variant.tsv')
var_info_file = file.path(projdir, 'results/ld/r2_groups/2017-07-04/0.95/output.tsv')
modsel_indir = file.path(projdir, 'results/pve/phenotype/lm/2017-07-20_forward')

chosen_traits = c('2_Aminoethanol', 'Formic_Acid', 'Anti_Gent', 'Anti_Spec',
                  'Anti_Strept', 'PEG_max', 'weight_A17', 'weight_R108', 
                  'bio_1', 'bio_12')
all_traits = c('2_Aminoethanol', 'D_Trehalose', 'Formic_Acid', 'L_Fucose',
               'N_Acetyl_D_Glucosamine', 'Putrescine', 'Anti_Gent', 'Anti_Spec',
               'Anti_Strept', 'metal_Cd', 'PEG_max', 'Temp_max', 'Salt_max',
               'weight_A17', 'weight_R108', 'nodule_A17', 'nodule_R108',
               'Growth_Rate', 'bio_1', 'bio_12')
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

## Read in the LD grouping assignments and calculate group
## size and proportion that are SNPs
ld = read.csv(ld_file, sep='\t', header=FALSE, as.is=TRUE) %>%
rename(chrom=V1, pos0=V2, pos=V3, rs_seed=V4, group=V5, members=V6) %>%
mutate(.,
       size=sapply(.[, 'members'], function(x) length(unlist(strsplit(x, ',')))),
       prop_snp=sapply(.[, 'members'],
                       function(x) {
                           m = unlist(strsplit(x, ','))
                           sum(grepl('snp', m)) / length(m)
                       }))

var_info = read.csv(var_info_file, sep='\t', header=TRUE, as.is=TRUE) %>%
    filter(seed == 1) %>%
    select(rs, maf)


## Read in the association results for each trait
## Only need the first 10 LD groups
lmm_data = structure(vector('list', length=length(all_traits)),
                     names=all_traits)
lmm_data_random = structure(vector('list', length=length(all_traits)),
                            names=all_traits)
for(trait in chosen_traits) {

    modsel_data = read.csv(file.path(modsel_indir, trait, 'output.stats.tsv'),
                           sep='\t', header=TRUE, as.is=TRUE)
    real_modsel_data = modsel_data[modsel_data[, 'run'] == 'real' &
                                   modsel_data[, 'n_vars_start'] == 10, ]
    random_modsel_data = modsel_data[modsel_data[, 'run'] != 'real' &
                                     modsel_data[, 'n_vars_start'] == 10, ]

    lmm_data[[trait]] = data.frame('rs'=unlist(strsplit(real_modsel_data[, 'vars_final'], ',')),
                                   stringsAsFactors=FALSE) %>%
        left_join(var_info) %>%
        rename(af=maf)

    for(i in 1:n_random) {
        lmm_data_random[[trait]][[as.character(i)]] =
            data.frame('rs'=unlist(strsplit(random_modsel_data[i, 'vars_final'], ',')),
                       stringsAsFactors=FALSE) %>%
            left_join(var_info) %>%
            rename(af=maf)
    }

}

## This will be a multi-page figure -- I think it has to be
pdf('variant_type_proportion_figure_2018-06-12.supp.pdf', width=6, height=9)
layout(matrix(ncol=4, data=1:(2*length(chosen_traits)), byrow=TRUE))
par(mgp=par('mgp')/1.8, mar=par('mar')/1.6)
for(trait in chosen_traits) {
    make_plot(trait, ld, lmm_data, lmm_data_random)
}
dev.off()
