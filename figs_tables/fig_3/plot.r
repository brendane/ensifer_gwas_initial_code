#!/usr/bin/env Rscript
#
# This makes the PRVE example figure using A17.
#
# 1. Forward model selection is performed on the top 25 variants from
#    the LMM association mapping. Adjusted R^2 is used as the criterion.
#    This is all done in a separate script.
# 2. For each variant, the adjusted R^2 from regressing that variant on
#    the residuals of the previously-added variants is calculated. Also
#    done in a separate script.
# 3. These residual R^2 estimates are plotted after subtracting the median
#    residual R^2 estimates from the random permutations.

library(data.table)
library(dplyr)
library(dtplyr)

make_pve_plot = function(trait, dat, n_random, n, ylab, r2_legend=FALSE,
                          cex_points=1, colr='black') {

    ## Set up the plot axes and legend.
    plot(dat[dat$run == 'real', 'n_vars_start'],
         dat[dat$run == 'real', 'r2'],
         xlim=c(0, n+1), ylim=c(0, 1.0), yaxs='i', xaxs='i',
         main='', xlab='Number of variants', ylab=ylab,
         type='n', xaxt='n', yaxt='n', bty='n')
    abline(h=0, col=colr)
    axis(side=2, at=seq(0, 1, 0.2),
         labels=c('0', '', '0.4', '', '', '1.0'),
         tick=TRUE, col=colr, las=2)
    axis(side=1, at=seq(0, n, 5), labels=seq(0, n, 5), tick=TRUE,
         col=colr)
    grid()

    ## Extract R^2 for the random permtuations for 1 - n
    ## top variants and plot.
    for(i in 0:(n_random-1)) {
        run = paste0('random-', i)
        y = dat[dat$run == run, 'r2'][1:n]
        lines(1:n, y, col='gray', lwd=1)
    }

    ## Extract the real R^2 values and plot with a line and points
    ## combination. The size of the points is proportional to the amount
    ## of genome tagged by each set of variants.
    d = dat[dat$run == 'real', 'r2'][1:n]
    lines(dat[dat$run == 'real', 'n_vars_start'][1:n],
          d,
          pch=1, col='black', lwd=2, xpd=NA)
}


make_prve_plot_a = function(trait, dat, n_random, n, ylab, r2_legend=FALSE,
                            cex_points=1, colr='black') {

    ## Set up the plot axes and legend.
    plot(dat[dat$run == 'real', 'n_vars_start'],
         dat[dat$run == 'real', 'r2_resid'],
         xlim=c(0, n+1), ylim=c(-0.2, 1.0), yaxs='i', xaxs='i',
         main='', xlab='Number of variants', ylab=ylab,
         type='n', xaxt='n', yaxt='n', bty='n')
    abline(h=-0.2, col=colr)
    axis(side=2, at=seq(-0.2, 1, 0.2),
         labels=c('-0.2', '', '', '0.4', '', '', '1.0'),
         tick=TRUE, col=colr, las=2)
    axis(side=1, at=seq(0, n, 5), labels=seq(0, n, 5), tick=TRUE,
         col=colr)
    grid()

    ## Extract residual R^2 for the random permtuations for 1 - n
    ## top variants.
    randoms = matrix(ncol=n, nrow=n_random)
    for(i in 0:(n_random-1)) {
        run = paste0('random-', i)
        randoms[i+1, ] = dat[dat$run == run, 'r2_resid'][1:n]
    }

    ## Calculate the median and lower 95% range for the random permutations
    ## at each number of top variants. Also, plot the ranges with gray
    ## lines.
    for(i in 1:n) {
        q0 = min(randoms[, i])
        q95 = quantile(randoms[, i], 0.95)
        q50 = quantile(randoms[, i], 0.50)
        segments(i, q0, i, q95, col='grey70', lwd=2)
        points(i, q0, col='grey70', cex=1, pch=15)
        points(i, q95, col='grey70', cex=1, pch=15)
    }

    ## Extract the real residual R^2 values and plot with a line and
    ## points combination. The size of the points is proportional to the
    ## amount of genome tagged by each set of variants.
    d = dat[dat$run == 'real', 'r2_resid'][1:n]
    lines(dat[dat$run == 'real', 'n_vars_start'][1:n],
          d,
          pch=1, col='black', lwd=2, xpd=NA)
}



make_prve_plot = function(trait, dat, n_random, n, ylab, r2_legend=FALSE,
                          cex_points=1, colr='black') {

    ## Set up the plot axes and legend.
    plot(dat[dat$run == 'real', 'n_vars_start'],
         dat[dat$run == 'real', 'r2_resid'],
         xlim=c(0, n+1), ylim=c(-0.2, 1.0), yaxs='i', xaxs='i',
         main='', xlab='Number of variants', ylab=ylab,
         type='n', xaxt='n', yaxt='n', bty='n')
    abline(h=-0.2, col=colr)
    if(r2_legend) {
        legend('topright',
               legend=c(0.01, 0.05, 0.1, 0.25, 0.5)*100,
               pt.cex=sqrt(c(0.01, 0.05, 0.1, 0.25, 0.5))*3*cex_points,
               pch=1, title='% Genome Tagged')

    }
    axis(side=2, at=seq(-0.2, 1, 0.2),
         labels=c('-0.2', '', '', '0.4', '', '', '1.0'),
         tick=TRUE, col=colr, las=2)
    axis(side=1, at=seq(0, n, 5), labels=seq(0, n, 5), tick=TRUE,
         col=colr)
    grid()

    ## Extract residual R^2 for the random permtuations for 1 - n
    ## top variants.
    randoms = matrix(ncol=n, nrow=n_random)
    for(i in 0:(n_random-1)) {
        run = paste0('random-', i)
        randoms[i+1, ] = dat[dat$run == run, 'r2_resid'][1:n]
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
        segments(i, q0-q50, i, q95-q50, col='grey70', lwd=2)
        points(i, q0-q50, col='grey70', cex=1, pch=15)
        points(i, q95-q50, col='grey70', cex=1, pch=15)
    }

    ## Extract the real residual R^2 values, subtract the median of the
    ## random permutations, and plot with a line and points combination.
    ## The size of the points is proportional to the amount of genome
    ## tagged by each set of variants.
    d = dat[dat$run == 'real', 'r2_resid'][1:n]
    lines(dat[dat$run == 'real', 'n_vars_start'][1:n],
          d-random_medians,
          pch=1, col='black', lwd=2, xpd=NA)
}

n_random = 100

projdir = file.path(Sys.getenv('HOME'), 'project', 'metab_gwas')
indir_modsel = file.path(projdir, 'results/pve/phenotype/lm/2017-07-20_forward')

data_modsel = read.csv(file.path(indir_modsel, 'weight_A17/output.stats.tsv'),
                       header=TRUE, sep='\t', as.is=TRUE)

pdf('model_selection_pve_figure_2018-02-23.example.pdf',
    width=4, height=9, useDingbats=FALSE)
par(mfrow=c(3, 1))
par(cex.axis=1.3, cex.lab=2)
n = 10

## PVE
make_pve_plot('weight_A17', data_modsel, n_random, n,
              'PVE', FALSE, 1, 'black')

## Residual R^2
make_prve_plot_a('weight_A17', data_modsel, n_random, n,
               'Proportion residual\nvariance explained (PRVE)',
               FALSE, 1, 'black')

## Final plot
make_prve_plot('weight_A17', data_modsel, n_random, n,
               'PRVE - median(permutations)', FALSE, 1, 'black')

dev.off()
