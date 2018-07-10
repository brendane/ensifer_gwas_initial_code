#!/usr/bin/env Rscript

projdir = file.path(Sys.getenv('HOME'), 'project', 'metab_gwas')

chosen_traits = c('weight_A17', 'Anti_Spec', 'PEG_max', 'Putrescine', 'nodule_A17')
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
                     'bio_1'='Mean Annual Temp.',
                     'bio_12'='Annual Precip.')

trait_lab = c('2_Aminoethanol'='Relative Absorbance', 
              'D_Trehalose'='Relative Absorbance',
              'Formic_Acid'='Relative Absorbance',
              'L_Fucose'='Relative Absorbance',
              'N_Acetyl_D_Glucosamine'='Relative Absorbance',
              'Putrescine'='Relative Absorbance',
              'Anti_Gent'=expression('Growth w/ 20 ' * mu * 'g/ml Gent.'),
              'Anti_Spec'=expression('Growth w/ 20 ' * mu * 'g/ml Spec.'),
              'Anti_Strept'=expression('Growth w/ 20 ' * mu * 'g/ml Strept.'),
              'metal_Cd'=expression('Growth w/ 20 ' * mu * 'g/ml Cd'),
              'PEG_max'='Max PEG4000 w/ Growth (%)',
              'Temp_max'='Max Temp w/ Growth (C)',
              'Salt_max'='Max NaCl w/ Growth (mM)',
              'weight_A17'='Adjusted dry shoot mass (g)',
              'weight_R108'='Adjusted dry shoot mass (g)',
              'nodule_A17'='Adjusted nodule number',
              'nodule_R108'='Adjusted nodule number',
              'Growth_Rate'='Growth Rate (doubling / hr)',
              'bio_1'='Mean Annual Temp.',
              'bio_12'='Annual Precip.')



x = read.csv(file.path(projdir,
                       'results/tabulate_phenotypes/2017-07-12/phenotypes.mel153.tsv'),
             sep='\t', header=TRUE, as.is=TRUE, check.names=FALSE)
x[, 'PEG_max'] = x[, 'PEG_max'] * 100

pdf('trait_distribution_figure_2018-05-21.main.pdf', width=3, height=9)
par(mfrow=c(5, 1), mgp=par('mgp')/1.8, mar=par('mar')/1.6)
for(i in seq_along(chosen_traits)) {
    trait = chosen_traits[i]
    hist(x[, trait], main='', xlab=trait_lab[[trait]], ylab='',
         col='grey20',
         border='grey40',
         breaks=15, xaxs='i', yaxs='i', lwd=0.25, cex.main=1.4,
         ylim=c(0, 120))
    abline(h=0, lwd=0.5)
    title(trait_long_names[trait], xpd=NA)
}

pdf('trait_distribution_figure_2018-05-21.supp.pdf', width=6, height=7.5)
par(mfrow=c(5, 4), mgp=par('mgp')/1.8, mar=par('mar')/1.6)
for(i in seq_along(all_traits)) {
    trait = all_traits[i]
    hist(x[, trait], main='', xlab=trait_lab[[trait]], ylab='',
         col={if(trait %in% chosen_traits){ 'red' }else{ 'grey20' }},
         border={if(trait %in% chosen_traits){ 'tomato' }else{ 'grey40' }},
         breaks=15, xaxs='i', yaxs='i', lwd=0.25, cex.main=1.4,
         ylim=c(0, 120))
    abline(h=0, lwd=0.5)
    title(trait_long_names[trait], xpd=NA)
}
