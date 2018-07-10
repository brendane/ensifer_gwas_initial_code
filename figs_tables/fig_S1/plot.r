#!/usr/bin/env Rscript
#
# Plot LD group size and calculate statistics. Also plot Venn diagrams
# of the type of variant in each LD group.
#
# This version separates single-variant groups from the others. Also, an
# update to the calculation of distance between variants was made.
#

options('nwarnings'=10000)

library(data.table)
library(gplots)

plot_colors = c('SNP'=rgb(142/255, 102/255, 62/255),
                'PAV'='blue',
                'Chr'=rgb(0.75, 0.5, 0.5),
                'pSymB'=rgb(0.5, 0.75, 0.5),
                'pSymA'=rgb(0.5, 0.5, 0.75))

replicon_lengths = c('Chr'=3671869,
                     'pSymB'=1680729,
                     'pSymA'=1363628,
                     'SNP'=10000000,
                     'PAV'=10000000)

circ_dist = function(x, y, l) {
    d = abs(y - x)
    d = ifelse(d > l / 2, l - d, d)
    d
}

calc_stats = function(memb) {
    n = length(unique(memb[['group']]))
    n_nosingle = length(unique(memb[['group']][memb[['size']] > 1]))
    sizes = memb[['size']]
    sizes_nosingle = sizes[memb[['size']] != 1]
    total = sum(sizes)
    x = 0
    n50 = 0
    for(s in sort(sizes_nosingle, decreasing=TRUE)) {
        x = x + s
        if(x >= total / 2) {
            n50 = s
            break
        }
    }
    list('n'=n, 
         'n_nosingle'=n_nosingle,
         'single'=sum(sizes == 1),
         'max'=max(sizes_nosingle), 'min'=min(sizes_nosingle),
         'mean'=mean(sizes_nosingle), 'median'=median(sizes_nosingle),
         'n50'=n50)
}

## Assumes all variants are on the same replicon
calc_dist_stats = function(ldg, groups, len) {
    ldg = ldg[ldg[['group']] %in% groups, ]
    aggregate(ldg[['pos']], list(ldg[['group']]),
              function(x) {
                  if(length(x) > 1) {
                      d = sapply(x, function(y) circ_dist(x, y, len))
                      d = d[upper.tri(d)]
                      c('min'=min(d), 'max'=max(d), 'mean'=mean(d),
                        'median'=median(d))
                  } else {
                      c('min'=NaN, 'max'=NaN, 'mean'=NaN, 'median'=NaN)
                  }
              })
}

projdir = file.path(Sys.getenv('HOME'), 'project/metab_gwas')
ld_groups_95 = fread(file.path(projdir, 'results/ld/r2_groups/2017-07-04/0.95/output.tsv'))
ld_groups_80 = fread(file.path(projdir, 'results/ld/r2_groups/2017-07-04/0.80/output.tsv'))

membs_95 = aggregate(ld_groups_95[['rs']],
                     list(ld_groups_95[['group']]),
                     function(r) {
                         c('snp'=any(grepl('snp', r)),
                           'pav'=any(grepl('rdv', r)),
                           'chr'=any(grepl('Uni0', r)),
                           'psymb'=any(grepl('Uni1', r)),
                           'psyma'=any(grepl('Uni2', r)),
                           'size'=length(r))
                     })
membs_95 = data.frame('group'=membs_95[[1]],
                      'SNP'=as.logical(membs_95[[2]][, 'snp']),
                      'PAV'=as.logical(membs_95[[2]][, 'pav']),
                      'Chr'=as.logical(membs_95[[2]][, 'chr']),
                      'pSymB'=as.logical(membs_95[[2]][, 'psymb']),
                      'pSymA'=as.logical(membs_95[[2]][, 'psyma']),
                      'size'=membs_95[[2]][, 'size'])

membs_80 = aggregate(ld_groups_80[['rs']],
                     list(ld_groups_80[['group']]),
                     function(r) {
                         c('snp'=any(grepl('snp', r)),
                           'pav'=any(grepl('rdv', r)),
                           'chr'=any(grepl('Uni0', r)),
                           'psymb'=any(grepl('Uni1', r)),
                           'psyma'=any(grepl('Uni2', r)),
                           'size'=length(r))
                     })
membs_80 = data.frame('group'=membs_80[[1]],
                      'SNP'=as.logical(membs_80[[2]][, 'snp']),
                      'PAV'=as.logical(membs_80[[2]][, 'pav']),
                      'Chr'=as.logical(membs_80[[2]][, 'chr']),
                      'pSymB'=as.logical(membs_80[[2]][, 'psymb']),
                      'pSymA'=as.logical(membs_80[[2]][, 'psyma']),
                      'size'=membs_80[[2]][, 'size'])

stats_95 = calc_stats(membs_95)
stats_80 = calc_stats(membs_80)


## Venn diagrams
pdf('ld_variant_type_figure_2018-06-26.pdf', width=5, height=6)
par(mfrow=c(2, 2))
venn(membs_95[membs_95[, 'size'] > 1, c('SNP', 'PAV')])
venn(membs_95[membs_95[, 'size'] > 1, c('Chr', 'pSymB', 'pSymA')])
venn(membs_80[membs_80[, 'size'] > 1, c('SNP', 'PAV')])
venn(membs_80[membs_80[, 'size'] > 1, c('Chr', 'pSymB', 'pSymA')])
dev.off()


## Size distribution of all LD groups, comparing two different
## R^2 thresholds.
pdf('ld_group_size_figure_2018-06-26.main.pdf', width=6, height=2.5)
par(mfrow=c(1, 2))
par('mgp'=par('mgp')/1.6, 'mar'=par('mar')/1.4)
h = hist(membs_95[membs_95[, 'size'] > 1, 'size'], breaks=seq(0, 10000, 1),
         plot=FALSE)
log_counts = log10(h[['counts']])
log_counts[is.infinite(log_counts)] = NaN
log_counts_plot = log_counts[1:15]
log_counts_plot[15] = log10(sum(h[['counts']][15:length(h[['counts']])]))
barplot(log_counts_plot, xlim=c(1, 15.5),
        space=0, col='gray20', border='grey50',
        names.arg=c(1:14, '≥15'),
        ylab='Number of LD groups',
        xlab='Size of LD group (# variants)',
        xaxs='i',  yaxt='n', yaxs='i', cex.axis=0.9,
        main=expression(R^2 * ' = 0.95'), ylim=c(0, 5))
axis(side=2, at=1:4, labels=c('10', '100', '1000', '10000'),
     cex=0.9)

stat_string = paste0('N. groups total = ', stats_95$n, '\n',
                     'N. groups = ', stats_95$n_nosingle, '\n',
                     'Median size = ', stats_95$median, '\n',
                     'Mean size = ', round(stats_95$mean, 1), '\n',
                     'Max. size = ', stats_95$max, '\n',
                     'N50 = ', stats_95$n50)
text(par('usr')[2] * 0.4, par('usr')[4] * 0.95, stat_string,
     adj=c(0, 1), cex=0.8, 'xpd'=NA)

h = hist(membs_80[membs_80[, 'size'] > 1, 'size'], breaks=seq(0, 10000, 1),
         plot=FALSE)
log_counts = log10(h[['counts']])
log_counts[is.infinite(log_counts)] = NaN
log_counts_plot = log_counts[1:15]
log_counts_plot[15] = log10(sum(h[['counts']][15:length(h[['counts']])]))
barplot(log_counts_plot, xlim=c(1, 15.5),
        space=0, col='gray20', border='grey50',
        names.arg=c(1:14, '≥15'),
        ylab='Number of LD groups',
        xlab='Size of LD group (# variants)',
        xaxs='i',  yaxt='n', yaxs='i', cex.axis=0.9,
        main=expression(R^2 * ' = 0.80'), ylim=c(0, 4.5))
axis(side=2, at=1:4, labels=c('10', '100', '1000', '10000'),
     cex=0.9)

stat_string = paste0('N. groups total = ', stats_80$n, '\n',
                     'N. groups = ', stats_80$n_nosingle, '\n',
                     'Median size = ', stats_80$median, '\n',
                     'Mean size = ', round(stats_80$mean, 1), '\n',
                     'Max. size = ', stats_80$max, '\n',
                     'N50 = ', stats_80$n50)
text(par('usr')[2] * 0.4, par('usr')[4] * 0.95, stat_string,
     adj=c(0, 1), cex=0.8, 'xpd'=NA)
dev.off()


## Supplemental figure with distance stats and divided into variant
## types. Just for R^2 = 0.95 to keep it simple.
## Note that only groups with size > 1 are used.
pdf('ld_group_size_figure_2018-06-26.supp.pdf', width=6, height=7.5)
par(mfrow=c(3, 2))
par('mgp'=par('mgp')/1.6, 'mar'=par('mar')/1.4)

cat('Number of variants by themselves at R^2 = 0.95\n')
cat('Total:', sum(membs_95[, 'size'] == 1), '\n')
cat('SNPs:', sum(membs_95[, 'size'] == 1 & membs_95[, 'SNP'] & !(membs_95[, 'PAV'])), '\n')
cat('PAVs:', sum(membs_95[, 'size'] == 1 & !(membs_95[, 'SNP']) & membs_95[, 'PAV']), '\n')

cat('\nNumber of variants by themselves at R^2 = 0.80\n')
cat('Total:', sum(membs_80[, 'size'] == 1), '\n')
cat('SNPs:', sum(membs_80[, 'size'] == 1 & membs_80[, 'SNP'] & !(membs_80[, 'PAV'])), '\n')
cat('PAVs:', sum(membs_80[, 'size'] == 1 & !(membs_80[, 'SNP']) & membs_80[, 'PAV']), '\n')

cat('\nNumber of SNPs by themselves at R^2 = 0.95\n')
cat('Total:', sum(membs_95[, 'size'] == 1), '\n')
cat('Chr SNPs:', sum(membs_95[, 'size'] == 1 & membs_95[, 'SNP'] & !(membs_95[, 'PAV']) &
                     (membs_95[, 'Chr']) & !(membs_95[, 'pSymB']) & !(membs_95[, 'pSymA'])), '\n')
cat('pSymB SNPs:', sum(membs_95[, 'size'] == 1 & membs_95[, 'SNP'] & !(membs_95[, 'PAV']) &
                       !(membs_95[, 'Chr']) & (membs_95[, 'pSymB']) & !(membs_95[, 'pSymA'])), '\n')
cat('pSymA SNPs:', sum(membs_95[, 'size'] == 1 & membs_95[, 'SNP'] & !(membs_95[, 'PAV']) &
                       !(membs_95[, 'Chr']) & !(membs_95[, 'pSymB']) & (membs_95[, 'pSymA'])), '\n')

cat('\nStats for groups with both SNPs and PAVs\n')
m = membs_95[(membs_95[, 'SNP']) & (membs_95[, 'PAV']) & (membs_95[, 'size'] > 1), ]
stats = calc_stats(m)
cat(paste0('N. groups = ', stats$n_nosingle, '\n',
           'Median size = ', stats$median, '\n',
           'Mean size = ', round(stats$mean, 1), '\n',
           'Max. size = ', stats$max, '\n',
           'N50 = ', stats$n50), '\n')


groups = vector('list', length=5)
names(groups) = c('SNP', 'PAV', 'Chr', 'pSymB', 'pSymA')
groups[['SNP']] = membs_95[membs_95[, 'SNP'] & !(membs_95[, 'PAV']) & (membs_95[, 'size'] > 1), 'group']
groups[['PAV']] = membs_95[membs_95[, 'PAV'] & !(membs_95[, 'SNP']) & (membs_95[, 'size'] > 1), 'group']
groups[['Chr']] = membs_95[membs_95[, 'Chr'] & !(membs_95[, 'PAV']) &
                           !(membs_95['pSymB']) & !(membs_95[, 'pSymA']) & (membs_95[, 'size'] > 1), 'group']
groups[['pSymB']] = membs_95[membs_95[, 'pSymB'] & !(membs_95[, 'PAV']) &
                             !(membs_95['Chr']) & !(membs_95[, 'pSymA']) & (membs_95[, 'size'] > 1), 'group']
groups[['pSymA']] = membs_95[membs_95[, 'pSymA'] & !(membs_95[, 'PAV']) &
                             !(membs_95['pSymB']) & !(membs_95[, 'Chr']) & (membs_95[, 'size'] > 1), 'group']

for(vt in names(groups)) {
    g = groups[[vt]]
    m = membs_95[membs_95[, 'group'] %in% g, ]
    rep_len = replicon_lengths[vt]

    ## Plot the distribution of LD group sizes for groups containing just
    ## the target variant type
    h = hist(m[, 'size'], breaks=seq(0, 10000, 1),
             plot=FALSE)
    log_counts = log10(h[['counts']])
    log_counts[is.infinite(log_counts)] = NaN
    log_counts_plot = log_counts[1:15]
    log_counts_plot[15] = log10(sum(h[['counts']][15:length(h[['counts']])]))
    barplot(log_counts_plot, xlim=c(1, 15.5),
            space=0, col=plot_colors[vt],
            border='grey50',
            names.arg=c(1:14, '≥15'),
            ylab='Number of LD groups',
            xlab='Size of LD group (# variants)',
            xaxs='i',  yaxt='n', yaxs='i', cex.axis=0.9,
            main=vt, ylim=c(0, 5))
    axis(side=2, at=1:4, labels=c('10', '100', '1000', '10000'),
         cex=0.9)

    ## Then calculate size statistics
    stats = calc_stats(m)

    ## And distance stats
    d_stats = calc_dist_stats(ld_groups_95, g, rep_len)

    ## Put the stats on the plot
    ## Note that the number of groups for this plot does not include
    ## groups with one variant.
    stat_string = paste0('N. groups = ', stats$n_nosingle, '\n',
                         'Median size = ', stats$median, '\n',
                         'Mean size = ', round(stats$mean, 1), '\n',
                         'Max. size = ', stats$max, '\n',
                         'N50 = ', stats$n50, '\n',
                         'Median(median dist.) = ', median(d_stats[[2]][, 'median'], na.rm=TRUE), '\n',
                         'Median(max dist.) = ', median(d_stats[[2]][, 'max'], na.rm=TRUE))
    text(par('usr')[2] * 0.4, par('usr')[4] * 0.95, stat_string,
         adj=c(0, 1), cex=0.8)
}

dev.off()

warnings()
