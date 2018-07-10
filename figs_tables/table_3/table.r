#!/usr/bin/env Rscript
#
# Make a table comparing the LMM ("BSLMM null") polygenic estimates of
# PVE to the estimates obtained from the linear model selection. The linear
# model selection estimate is from the number of variants with the greatest
# adjusted R^2 after subtracting the median of permutations. The median of
# permutations is also subtracted from the LMM estimate.
#
# This uses cumulative PVE rather than residual PVE, because we are looking
# for total variance explained.
#

projdir = file.path(Sys.getenv('HOME'), 'project', 'metab_gwas')

chosen_traits = c('2_Aminoethanol', 'Formic_Acid', 'Anti_Gent', 'Anti_Spec',
                  'Anti_Strept', 'PEG_max', 'weight_A17', 'weight_R108',
                  'bio_1', 'bio_12')
traits = c('2_Aminoethanol', 'D_Trehalose', 'Formic_Acid', 'L_Fucose',
           'N_Acetyl_D_Glucosamine', 'Putrescine', 'Anti_Gent', 'Anti_Spec',
           'Anti_Strept', 'metal_Cd', 'PEG_max', 'Temp_max', 'Salt_max',
           'weight_A17', 'weight_R108', 'nodule_A17', 'nodule_R108',
           'Growth_Rate', 'bio_1', 'bio_12')
trait_long_names = c('2_Aminoethanol'='2-Aminoethanol', 
                     'D_Trehalose'='D-Trehalose',
                     'Formic_Acid'='Formic Acid',
                     'L_Fucose'='L-Fucose',
                     'N_Acetyl_D_Glucosamine'='N-Acetyl-D-Glucosamine',
                     'Putrescine'='Putrescine',
                     'Anti_Gent'='Gentamicin Resistance',
                     'Anti_Spec'='Spectinomycin Resistance',
                     'Anti_Strept'='Streptomycin Resistance',
                     'metal_Cd'='Cadmium Tolerance',
                     'PEG_max'='PEG Tolerance',
                     'Temp_max'='Max. Growth Temp.',
                     'Salt_max'='NaCl Tolerance',
                     'weight_A17'='A17 Biomass',
                     'weight_R108'='R108 Biomass',
                     'nodule_A17'='A17 Nodule Number',
                     'nodule_R108'='R108 Nodule Number',
                     'Growth_Rate'='Growth Rate',
                     'bio_1'='Annual Mean Temp.',
                     'bio_12'='Annual Precip.')

### TOP VARIANTS -- LINEAR MODEL SELECTION

## Read in the data
lm_indir = file.path(projdir, 'results/pve/phenotype/lm/2017-07-20_forward')
lm_stats = structure(vector('list', length=length(traits)), names=traits)
for(trait in traits) {
    lm_stats[[trait]] = read.csv(file.path(lm_indir, trait, 'output.stats.tsv'),
                                 header=TRUE, as.is=TRUE, check.names=FALSE, sep='\t')
}

## Get the maximum PVE of top variants after subtracting the median
## of the permutations
max_top_pve = data.frame(trait=traits,
                         'max_pve'=numeric(length(traits)),
                         'number_of_vars_max_pve'=numeric(length(traits)),
                         stringsAsFactors=FALSE)
rownames(max_top_pve) = traits
for(trait in traits) {
    s = lm_stats[[trait]]
    s_random = s[s[, 'run'] != 'real', ]
    random_medians = aggregate(s_random[, 'r2'],
                               list(s_random[, 'n_vars_start']),
                               median)
    real_above_random = s[s[, 'run'] == 'real', 'r2'] - random_medians[['x']]
    max_top_pve[trait, 'max_pve'] = max(real_above_random)
    max_top_pve[trait, 'number_of_vars_max_pve'] = which(real_above_random == max(real_above_random))[1]
}


### LMM polygenic estimate (skip BSLMM b/c I don't have permutations for all
### of the traits)

## Read in the data
real_poly_dir = file.path(projdir, 'results/pve/phenotype/bslmm/2017-07-17')
random_poly_dir = file.path(projdir, 'results/pve/phenotype/bslmm/2017-07-17_random')
random_poly_dir2 = file.path(projdir, 'results/pve/phenotype/bslmm/2017-07-17_20_random')

bslmm_data = structure(vector('list', length=length(traits)),
                       names=traits)
for(trait in traits) {
    real_stats = read.csv(file.path(real_poly_dir, trait, 'stats.txt'),
                          sep='\t', header=TRUE, check.names=FALSE, as.is=TRUE)
    real_stats[1, 'run'] = 'real'
    random_file = file.path(random_poly_dir, trait, 'stats.txt')
    if(!file.exists(random_file)) {
        random_file = file.path(random_poly_dir2, trait, 'stats.txt')
    }
    random_stats = read.csv(random_file,
                            sep='\t', header=TRUE, check.names=FALSE, as.is=TRUE)
    random_stats[, 'run'] = paste0('random-', random_stats[, 'run'])
    bslmm_data[[trait]] = rbind(real_stats, random_stats)
}

## Get the estimates
lmm_pve = data.frame(trait=traits,
                     'lmm_pve'=numeric(length(traits)),
                     stringsAsFactors=FALSE)
rownames(lmm_pve) = traits
for(trait in traits) {
    s = bslmm_data[[trait]]
    real_pve = s[s[, 'run'] == 'real', 'pve_mean']
    random_median = median(s[s[, 'run'] != 'real', 'pve_mean'])
    lmm_pve[trait, 'lmm_pve'] = real_pve - random_median
}


### Combine
m = merge(max_top_pve, lmm_pve)
m = m[m[, 'trait'] %in% chosen_traits, ]
m[, 'lmm_minus_top'] = m[, 'lmm_pve'] - m[, 'max_pve']
m[, 'top_div_lmm'] = m[, 'max_pve'] / m[, 'lmm_pve']
m[, 'max_pve'] = round(m[, 'max_pve'], 3)
m[, 'lmm_pve'] = round(m[, 'lmm_pve'], 3)
m[, 'lmm_minus_top'] = round(m[, 'lmm_minus_top'], 3)
m[, 'top_div_lmm'] = round(m[, 'top_div_lmm'], 3)
m[, 'trait'] = trait_long_names[m[, 'trait']]

write.table(m, file='lm_lmm_pve_comparison_table_2018-05-24.tsv',
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
