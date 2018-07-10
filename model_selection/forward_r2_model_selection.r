#!/usr/bin/env Rscript
###
### Perform simple forward model selection, always choosing to add the
### variant that increases adjusted R2 by the greatest amount (or
### decreases adj. R2 by the least).
###
### Imputes missing genotype data by using the mean; ignores strains with
### missing phenotype data.
###

library(data.table)
library(MASS)

argv = commandArgs(trailingOnly=TRUE)
genotype_file = argv[1]
phenotype_file = argv[2]
variant_file = argv[3]
n_vars = as.numeric(argv[4])
min_maf = as.numeric(argv[5])
output_prefix = argv[6]

oh = file(paste0(output_prefix, '.stats.tsv'), 'wb')
ohf = file(paste0(output_prefix, '.fitted.tsv'), 'wb')

n_random_vars = 10000

genotypes = fread(genotype_file, header=TRUE)
setkey(genotypes, 'rs')
phenotypes = fread(phenotype_file)
variants = scan(variant_file, what='character')

## Filter by MAF
maf = apply(genotypes[, -(1:3)], 1,
            function(x) min(sum(x == 0, na.rm=TRUE), sum(x == 1, na.rm=TRUE))/sum(!is.na(x)))
genotypes = genotypes[maf >= min_maf, ]
cat(nrow(genotypes), 'variants left after MAF filtering\n', file=stderr())

if(length(variants) < n_vars) {
    stop("Can't choose that many variants")
}

## Sort genotype data by the order in variant_file and 
## impute missing genotype data.
X = apply(t(as.matrix(genotypes[variants, -(1:3)])), 2,
          function(x) { x[is.na(x)] = mean(x, na.rm=TRUE); x})
p = phenotypes[match(rownames(X), as.character(phenotypes[[ 2]])), ]
y = p[[3]]

## Choose random variants for calculating LD
X_rand = t(as.matrix(genotypes[sample(1:nrow(genotypes), n_random_vars, FALSE), -(1:3)]))

## Print output file headers
cat('n_vars_start\tr2\tr2_unadj\tr2_resid\tgenomic_r2\tgenomic_r2_unadj\tn_vars_final\t',
    'ranks_final\tvars_final\tcoefficients\n',
    file=oh, sep='')
cat('n_vars_start\tphenotype\tfitted\n', file=ohf)

### Run the model selection
vars_included = character(n_vars)
XX = as.matrix(X[, 1:n_vars]) # Genotype data with just the top n_vars variants
colnames(XX) = paste0('XX', 1:n_vars)
XX = as.matrix(XX[, !is.na(apply(XX, 2, mean))])
for(nv in 1:n_vars) {

    ## Figure out which variant gives the greatest R2
    max_r2a = -Inf
    current_var = ''
    for(i in 1:n_vars) {
        test_var_name = paste0('XX', i)
        if(test_var_name %in% vars_included) {
            next
        }
        vars_in_this_model = c(vars_included[vars_included != ''],
                               test_var_name)
        XXX = XX[, vars_in_this_model, drop=FALSE]
        if(ncol(XXX) == 0) {
            next
        }
        model = lm(y ~ XXX)
        r2a_this_variant = summary(model)$adj.r.squared
        if(r2a_this_variant > max_r2a) {
            current_var = test_var_name
            max_r2a = r2a_this_variant
        }
    }

    ## If nothing found
    if(current_var == '') {
        break
    }

    ## Get the stats for the best model
    vars_included[nv] = current_var
    vars_in_this_model = vars_included[vars_included != '']
    XXX = XX[, vars_in_this_model]
    model = lm(y ~ XXX)
    s = summary(model)
    r2 = s$r.squared
    r2a = s$adj.r.squared
    rhs_final = as.character(model$call$formula)[3]
    ranks_final = as.numeric(gsub('XXXXX', '',
                                  names(model$coefficients)[-1]))

    ## Measure the R2 on the residuals from the previous model
    if(nv > 1) {
        previous_vars_model = vars_in_this_model[1:(nv-1)]
        pv = XX[, previous_vars_model]
        rs = residuals(lm(y ~ pv))
        cv = XX[!is.na(y), current_var]
        model2 = lm(rs ~ cv)
        r2_resid = summary(model2)$adj.r.squared
    } else {
        r2_resid = r2a
    }

    coefs_final = model$coefficients
    variants_final = variants[ranks_final]
    ft = model$fitted.values
    gr2 = 0
    gr2c = 0
    gr2a = 0
    gr2ac = 0
    for(i in 1:n_random_vars) {
        ss = summary(lm(paste('X_rand[, i] ~', rhs_final)))
        if(!is.na(ss$r.squared)) {
            gr2 = gr2 + ss$r.squared
            gr2c = gr2c + 1
        }
        if(!is.na(ss$adj.r.squared)) {
            gr2a = gr2a + ss$adj.r.squared
            gr2ac = gr2ac + 1
        }
    }

    cat(nv, '\t', round(r2a, 3), '\t', round(r2, 3), 
        '\t', round(r2_resid, 3),
        '\t', round(gr2a/gr2ac, 3), '\t', round(gr2/gr2c, 3),
        '\t', length(variants_final),
        '\t', paste(ranks_final, collapse=','),
        '\t', paste(variants_final, collapse=','),
        '\t', paste(coefs_final, collapse=','),
        '\n', sep='', file=oh)
    write.table(cbind(rep(nv, sum(!is.na(y))), y[!is.na(y)], ft),
                col.names=FALSE, row.names=FALSE, quote=FALSE,
                sep='\t', file=ohf)
    flush(oh)
    flush(ohf)
}

close(oh)
close(ohf)
