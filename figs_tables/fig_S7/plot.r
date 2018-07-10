#!/usr/bin/env Rscript
#
# Make QQ plots for 20 traits.
#

library(data.table)
library(qqman)

chosen_traits = c('weight_A17', 'Anti_Spec', 'PEG_max', 'Putrescine', 'nodule_A17')
pve_traits = c('weight_A17', 'Anti_Spec', 'PEG_max',
               '2_Aminoethanol', 'Formic_Acid', 'Anti_Gent', 'Anti_Strept',
               'weight_R108', 'bio_1', 'bio_12')
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
                     'PEG_max'='Dessication\nTolerance',
                     'Temp_max'='Max. Growth Temp.',
                     'Salt_max'='NaCl Tolerance',
                     'weight_A17'='A17 Biomass',
                     'weight_R108'='R108 Biomass',
                     'nodule_A17'='A17 Nodule\nNumber',
                     'nodule_R108'='R108 Nodule\nNumber',
                     'Growth_Rate'='Growth Rate',
                     'bio_1'='Annual Mean Temp.',
                     'bio_12'='Annual Precip.')



projdir = file.path(Sys.getenv('HOME'), 'project', 'metab_gwas')
indir_lmm = file.path(projdir, 'results/gwas/phenotype/gemma/2017-07-17')

data_lmm = structure(vector('list', length=length(all_traits)),
                     names=all_traits)
for(trait in all_traits) {
    data_lmm[[trait]] = fread(paste('zcat', file.path(indir_lmm, trait, 'output',
                                                      'output0.assoc.txt.gz')))
}     


png('qq_figure_2018-06-12.supp.png', width=6*720, height=7.5*720, res=720)
par(mfrow=c(5, 4), mgp=par('mgp')/1.6, mar=par('mar')/1.3)
k = 1
for(trait in all_traits) {
    qq(data_lmm[[trait]][['p_lrt']], main='', col='black')
    title(trait_long_names[[trait]],
          col.main={if(trait %in% chosen_traits) 'red' else 'black'},
          font.main={if(trait %in% pve_traits) 2 else 1})
    k = k + 1
}
