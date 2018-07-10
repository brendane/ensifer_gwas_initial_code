#!/usr/bin/env Rscript
#
# This makes the PRVE figure -- the curves showing variance explained
# by the top 15 variants.
#
# 1. Forward model selection is performed on the top 25 variants from
#    the LMM association mapping. Adjusted R^2 is used as the criterion.
#    This is all done in a separate script.
# 2. For each variant, the adjusted R^2 from regressing that variant on
#    the residuals of the previously-added variants is calculated. Also
#    done in a separate script. These values are plotted (previous version
#    subtracted median of random).
#
# This script also determines which traits have evidence for top variants
# explaining more variance than chance. It looks to see if the first
# variant explains more variation than 95% of random permutations.
#

library(data.table)
library(dplyr)
library(dtplyr)

make_prve_plot = function(trait, dat, n_random, n, trait_name, r2_legend=FALSE,
                          cex_points=1, colr='black', small) {
    if(small) {
        lwd_line = 1
        lwd_seg = 2
        cex_rand = 0.75
    } else {
        lwd_line = 2
        lwd_seg = 2
        cex_rand = 1
    }

    ## Set up the plot axes and legend.
    plot(dat[dat$run == 'real', 'n_vars_start'],
         dat[dat$run == 'real', 'r2_resid'],
         xlim=c(0, n+1), ylim=c(0, 0.65), yaxs='i', xaxs='i',
         xlab='Number of variants', ylab='PRVE',
         type='n', xaxt='n', yaxt='n', bty='n')
    abline(h=0, col=colr)
    axis(side=2, at=seq(-0.2, 1, 0.2),
         labels=c('-0.2', '', '', '0.4', '', '', '1.0'),
         tick=TRUE, col=colr, las=2)
    axis(side=1, at=seq(0, n, 5), labels=seq(0, n, 5), tick=TRUE,
         col=colr)
    grid()

    ## Extract residual R^2 for the random permtuations for 1 - n
    ## top variants.
    randoms = matrix(ncol=n, nrow=n_random)
    randoms_cum = matrix(ncol=n, nrow=n_random)
    for(i in 0:(n_random-1)) {
        run = paste0('random-', i)
        randoms[i+1, ] = dat[dat$run == run, 'r2_resid'][1:n]
        randoms_cum[i+1, ] = dat[dat$run == run, 'r2'][1:n]
    }

    ## Calculate the median and lower 95% range for the random permutations
    ## at each number of top variants. Also, plot the ranges with gray
    ## lines.
    random_medians = numeric(n)
    for(i in 1:n) {
        q0 = min(randoms[, i])
        q95 = quantile(randoms[, i], 0.95)
        q50 = quantile(randoms[, i], 0.50)
        random_medians[i] = q50
        segments(i, q0, i, q95, col='grey70', lwd=lwd_seg)
        if(q95 < dat[dat$run == 'real', 'r2_resid'][i]) {
            cat(trait, i, 'greater residual r2 than random\n')
        }

        q95_cum = quantile(randoms_cum[, i], 0.95)
        if(q95_cum < dat[dat$run == 'real', 'r2'][i]) {
            cat(trait, i, 'greater cumulative r2 than random\n')
        }
    }

    ## Extract the real residual R^2 values, and plot with a line.
    d = dat[dat$run == 'real', 'r2_resid'][1:n]
    lines(dat[dat$run == 'real', 'n_vars_start'][1:n],
          d,
          pch=1, col='black', lwd=lwd_line, xpd=NA)

    ## Add title
    more_pve_than_chance = FALSE
    q95_1 = quantile(randoms_cum[, 1], 0.95)
    more_pve_than_chance = (q95_1) < (dat[dat$run == 'real', 'r2'][1])
    title(trait_name, col.main=colr, font.main={if(more_pve_than_chance) 2 else 1})
}

n_random = 100

projdir = file.path(Sys.getenv('HOME'), 'project', 'metab_gwas')
indir_modsel = file.path(projdir, 'results/pve/phenotype/lm/2017-07-20_forward')

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

data_modsel = structure(vector('list', length=length(all_traits)),
                        names=all_traits)

for(trait in all_traits) {
    ## Linear model selection
    data_modsel[[trait]] = read.csv(file.path(indir_modsel,
                                              trait, 'output.stats.tsv'),
                                    header=TRUE, sep='\t', as.is=TRUE)
}

## Make the figure with just five traits
pdf('model_selection_pve_figure_2018-06-11.main.pdf', width=3, height=9)
par(mfrow=c(5, 1))
par('mgp'=par('mgp')/1.6, 'mar'=par('mar')/1.4,
    cex.axis=1.3, cex.lab=1.5)
n = 10
for(trait in chosen_traits) {
    dat = data_modsel[[trait]]
    make_prve_plot(trait, dat, n_random, n, trait_long_names[trait],
                   trait==chosen_traits[1], 1, 'black', FALSE)
}
dev.off()

## Make the figure with all the traits
pdf('model_selection_pve_figure_2018-06-11.supp.pdf', width=6, height=7.5)
par(mfrow=c(5, 4))
par('mgp'=par('mgp')/1.6, 'mar'=par('mar')/1.4)
n=10
for(trait in all_traits) {
    dat = data_modsel[[trait]]
    make_prve_plot(trait, dat, n_random, n, trait_long_names[trait],
                   trait==all_traits[1], 0.8,
                   ifelse(trait %in% chosen_traits, 'red', 'black'),
                   TRUE)
}
dev.off()
