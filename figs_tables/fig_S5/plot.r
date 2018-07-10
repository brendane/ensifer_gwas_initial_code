#!/usr/bin/env Rscript
#
# Figures showing variance explained by the K-matrix and by BSLMM.
#

make_plot = function(trait, pve_real, pve_real_snp, pve_real_pav,
                     pve_random, pve_random_snp, pve_random_pav,
                     lgnd=FALSE, show_bslmm_random=FALSE) {

    ## Set up the plotting region and axes
    plot(1:8, rep(0, 8),
         xlim=c(0, 8), ylim=c(0, 1), yaxs='i', xaxs='i',
         main=trait_long_names[trait], xlab='', ylab='PVE',
         type='n', xaxt='n', bty='n', las=2)
    axis(side=1, at=c(2, 6), lwd=0,
         labels=c('LMM', 'BSLMM'), tick=TRUE)
    abline(h=seq(0.2, 1, 0.2), col='lightgray', lty=3, lwd=0.5)
    abline(h=0)
    abline(v=4)

    ## Plot the real LMM PVE estimates as colored points on the left
    ## side of the plot
    points(1, pve_real[[trait]][['lmm_pve']], col='black',
           pch=15, cex=2, xpd=NA)
    points(2, pve_real_snp[[trait]][['lmm_pve']], col=rgb(142/255, 102/255, 62/255),
           pch=15, cex=2, xpd=NA)
    points(3, pve_real_pav[[trait]][['lmm_pve']], col='blue',
           pch=15, cex=2, xpd=NA)

    ## Plot the permuted LMM PVE estimates as gray lines
    segments(1, min(pve_random[[trait]][, 'lmm_pve']), 1,
             quantile(pve_random[[trait]][, 'lmm_pve'], 0.95),
             lty=1, col='grey40', lwd=2)
    segments(2, min(pve_random_snp[[trait]][, 'lmm_pve']), 2,
             quantile(pve_random_snp[[trait]][, 'lmm_pve'], 0.95),
             lty=1, col='grey40', lwd=2)
    segments(3, min(pve_random_pav[[trait]][, 'lmm_pve']), 3,
             quantile(pve_random_pav[[trait]][, 'lmm_pve'], 0.95),
             lty=1, col='grey40', lwd=2)

        ## Plot the real BSLMM point estimates as colored points on the right
        ## side of the plot
        points(5, pve_real[[trait]][['pve_mean']], col='black',
               pch=15, cex=2, xpd=NA)
        points(6, pve_real_snp[[trait]][['pve_mean']], col=rgb(142/255, 102/255, 62/255),
               pch=15, cex=2, xpd=NA)
        points(7, pve_real_pav[[trait]][['pve_mean']], col='blue',
               pch=15, cex=2, xpd=NA)

    if(show_bslmm_random) {
        ## And the random BSLMM estimates
        segments(5, min(pve_random[[trait]][, 'pve_mean']), 5,
                 quantile(pve_random[[trait]][, 'pve_mean'], 0.95),
                 lty=1, col='grey40', lwd=2)
        segments(6, min(pve_random_snp[[trait]][, 'pve_mean']), 6,
                 quantile(pve_random_snp[[trait]][, 'pve_mean'], 0.95),
                 lty=1, col='grey40', lwd=2)
        segments(7, min(pve_random_pav[[trait]][, 'pve_mean']), 7,
                 quantile(pve_random_pav[[trait]][, 'pve_mean'], 0.95),
                 lty=1, col='grey40', lwd=2)
    }

    if(lgnd) {
        legend('topright', legend=c('All', 'SNP', 'PAV'),
               fill=c('black', rgb(142/255, 102/255, 62/255), 'blue'),
               border=NA)
    }
}


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

n_random = 100

projdir = file.path(Sys.getenv('HOME'), 'project', 'metab_gwas')

indir_all = file.path(projdir, 'results/pve/phenotype/bslmm/2017-07-17')
indir_snp = file.path(projdir, 'results/pve/phenotype/bslmm/2017-07-17_separate/snp')
indir_pav = file.path(projdir, 'results/pve/phenotype/bslmm/2017-07-17_separate/rdv')
indir_all_random = file.path(projdir, 'results/pve/phenotype/bslmm/2017-07-17_random')
indir_snp_random = file.path(projdir, 'results/pve/phenotype/bslmm/2017-07-17_separate_random/snp')
indir_pav_random = file.path(projdir, 'results/pve/phenotype/bslmm/2017-07-17_separate_random/rdv')
indir_all_random2 = file.path(projdir, 'results/pve/phenotype/bslmm/2017-07-17_20_random')
indir_snp_random2 = file.path(projdir, 'results/pve/phenotype/bslmm/2017-07-17_20_separate_random/snp')
indir_pav_random2 = file.path(projdir, 'results/pve/phenotype/bslmm/2017-07-17_20_separate_random/rdv')

pve_real = structure(vector('list', length=length(all_traits)),
                     names=all_traits)
pve_real_snp = structure(vector('list', length=length(all_traits)),
                         names=all_traits)
pve_real_pav = structure(vector('list', length=length(all_traits)),
                         names=all_traits)
pve_random = structure(vector('list', length=length(all_traits)),
                       names=all_traits)
pve_random_snp = structure(vector('list', length=length(all_traits)),
                           names=all_traits)
pve_random_pav = structure(vector('list', length=length(all_traits)),
                           names=all_traits)

for(trait in all_traits) {

    ## Real PVE estimates for BSLMM and the LMM null model taken from
    ## the BSLMM runs
    pve_real[[trait]] = read.csv(file.path(indir_all, trait, 'stats.txt'),
                                 header=TRUE, as.is=TRUE, sep='\t')
    pve_real_snp[[trait]] = read.csv(file.path(indir_snp, trait, 'stats.txt'),
                                     header=TRUE, as.is=TRUE, sep='\t')
    pve_real_pav[[trait]] = read.csv(file.path(indir_pav, trait, 'stats.txt'),
                                     header=TRUE, as.is=TRUE, sep='\t')
    
    ## PVE estimates from permutations. Two different sets of runs -- the
    ## BSLMM runs for the non-focal traits were really short and should
    ## not be trusted, but the LMM "null" model estimates are fine.
    rid = indir_all_random
    rids = indir_snp_random
    ridp = indir_pav_random
    if(!(trait %in% chosen_traits)) {
        rid = indir_all_random2
        rids = indir_snp_random2
        ridp = indir_pav_random2
    }
    pve_random[[trait]] = read.csv(file.path(rid, trait, 'stats.txt'),
                                 header=TRUE, as.is=TRUE, sep='\t')
    pve_random_snp[[trait]] = read.csv(file.path(rids, trait, 'stats.txt'),
                                     header=TRUE, as.is=TRUE, sep='\t')
    pve_random_pav[[trait]] = read.csv(file.path(ridp, trait, 'stats.txt'),
                                     header=TRUE, as.is=TRUE, sep='\t')
}

pdf('bslmm_lmm_pve_figure_2018-05-21.supp.pdf', width=6, height=7.5)
par(mfrow=c(5, 4))
par('mgp'=par('mgp')/1.6, 'mar'=par('mar')/1.4)
for(trait in all_traits) {
    make_plot(trait, pve_real, pve_real_snp, pve_real_pav,
              pve_random, pve_random_snp, pve_random_pav,
              (trait==all_traits[1]), (trait %in% chosen_traits))
}
dev.off()
