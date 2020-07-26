#!/usr/bin/env Rscript
#
# Compare permutations to actual values.
#

N_PERMUTATIONS = 1000

projdir = '/home/tiffinp/epste051/project/symbiosis_gene_evolution' 

indir = file.path(projdir, 'results/gene_comparison/73strains_alpha_complete/2020-05-11/results')

sym_sets = c('symbiosis', 'fixation', 'signaling', 'gwas', 'biomass',
             'rhizobial_fitness')
sym_sets_names = c('fixation'='Fixation', 'signaling'='Signaling',
                   'biomass'='Host Benefit',
                   'rhizobial_fitness'='Symbiont Fitness')
stats = c('n_genes', 'strain_count', 'copy_number', 'dup_rate', 'loss_rate',
          'transfer_rate', 'transfer_minus_dup_rate', 'log2_median_kaks', 'log2_median_kaks_3',
          'r2_scc', 'delta', 'delta_pav', 'delta_pav_3')

breaks = list('transfer_rate'=c(0, 0.6),
              'strain_count'=c(0, 27),
              'r2_scc'=c(0, 1),
              'log2_median_kaks'=c(-4.5, -1.5),
              'log2_median_kaks_3'=c(-4.5, -1.5),
              'dup_rate'=c(0, 0.4),
              'loss_rate'=c(0, 1),
              'transfer_minus_dup_rate'=c(-0.6, 0.6),
              'delta'=c(0.5, 8.1),
              'delta_pav'=c(0.5, 8.1),
              'delta_pav_3'=c(0.15, 8.15),
              'copy_number'=c(1, 3),
              'n_genes'=c(1, 300))


output = structure(vector('list', length=2), names=c('all', 'matched'))
for(comp in names(output)) {
    output[[comp]] = matrix(nrow=length(stats), ncol=length(sym_sets), data=NaN,
                            dimnames=list(stats, sym_sets))
    for(st in stats) {
        pdf(paste0(st, '.', comp, '.pdf'), height=4, width=4)
        par(mfrow=c(2,2), mgp=c(1.5, 0.5, 0), mar=c(2.1, 3.1, 3.1, 2.1))
        for(ss in sym_sets) {
            func = median
            if(st == 'copy_number') func = mean
            fname = file.path(indir, ss, paste0('table.', comp, '.', st, '.tsv'))
            if(file.exists(fname)) {
                handle = file(fname, 'r')
                null = scan(handle, nlines=1, quiet=TRUE, sep='\t', what='character')
                real_value = func(scan(handle, nlines=1, quiet=TRUE, sep='\t'), na.rm=TRUE)
                permuted_values = numeric(N_PERMUTATIONS) * NaN
                for(i in seq_along(permuted_values)) {
                    null = scan(handle, nlines=1, quiet=TRUE, sep='\t', what='character')
                    permuted_values[i] = func(scan(handle, nlines=1, quiet=TRUE, sep='\t'), na.rm=TRUE)
                }
                close(handle)
                permuted_values = permuted_values[!is.na(permuted_values)]
                output[[comp]][st, ss] = sum(permuted_values > real_value) / length(permuted_values)
                if(ss %in% c('fixation', 'biomass', 'rhizobial_fitness', 'signaling')) {
                    pv = permuted_values
                    if(max(pv) > max(breaks[[st]])) pv[pv > max(breaks[[st]])] = max(breaks[[st]]) # Only for a few delta values
                    hist(pv, col='gray', border='white', main=sym_sets_names[ss],
                         cex.axis=1.25, cex.main=1.6, cex.lab=1.5, xaxs='i',
                         yaxs='i', ylab='', xlab='',
                         breaks=seq(breaks[[st]][1], breaks[[st]][2], length.out=20))
                    abline(v=real_value, lwd=2, col='gray10', lty=2)
                }
            } else {
                print(paste('Skipping', ss, st))
            }
        }
    }
    dev.off()
}


## Look at correlations between median Ka/Ks and some indications of HGT
## using the all and matched samples
kaks_cor = structure(vector('list', length=2), names=c('all', 'matched'))
for(comp in names(kaks_cor)) {
    kaks_cor[[comp]] = matrix(nrow=6, ncol=length(sym_sets), data=NaN,
                            dimnames=list(c('dup_rate', 'transfer_rate', 'r2_scc', 'delta', 'delta_pav', 'delta_pav_3'), sym_sets))
    for(st in rownames(kaks_cor[[comp]])) {
        for(ss in sym_sets) {
            fname = file.path(indir, ss, paste0('table.', comp, '.', st, '.tsv'))
            fname2 = file.path(indir, ss, paste0('table.', comp, '.log2_median_kaks_3.tsv'))
            if(file.exists(fname)) {
                handle = file(fname, 'r')
                handle2 = file(fname2, 'r')
                agenes = scan(handle, nlines=1, quiet=TRUE, sep='\t', what='character')
                bgenes = scan(handle2, nlines=1, quiet=TRUE, sep='\t', what='character')
                a = scan(handle, nlines=1, quiet=TRUE, sep='\t')
                b = scan(handle2, nlines=1, quiet=TRUE, sep='\t')
                a = a[match(bgenes, agenes)]
                real_value = cor(a, b, use='complete.obs')
                permuted_values = numeric(N_PERMUTATIONS) * NaN
                for(i in seq_along(permuted_values)) {
                    try({
                        agenes = scan(handle, nlines=1, quiet=TRUE, sep='\t', what='character')
                        bgenes = scan(handle2, nlines=1, quiet=TRUE, sep='\t', what='character')
                        a = scan(handle, nlines=1, quiet=TRUE, sep='\t')
                        b = scan(handle2, nlines=1, quiet=TRUE, sep='\t')
                        a = a[match(bgenes, agenes)]
                        permuted_values[i] = cor(a, b, use='complete.obs')
                    })
                }
                close(handle)
                permuted_values = permuted_values[!is.na(permuted_values)]
                kaks_cor[[comp]][st, ss] = sum(permuted_values > real_value, na.rm=TRUE) / length(permuted_values)
            } else {
                print(paste('Skipping', ss, st))
            }
        }
    }
    dev.off()
}


write_readable_table = function(object, fname=NULL, handle=NULL) {
    if(is.null(handle)) {
        handle_ = file(fname, 'w')
    } else {
        handle_ = handle
    }
    for(i in 1:ncol(object)) {
        cat('\n\t', file=handle_)
        write.table(object[, i, drop=FALSE], handle_, sep='\t',
                    col.names=TRUE, row.names=TRUE, quote=FALSE)
    }
    if(is.null(handle)) close(handle_)
}

handle = file('output.tsv', 'w')
cat('All:\n', file=handle)
write_readable_table(output[['all']], handle=handle)
cat('\nMatched:\n', file=handle)
write_readable_table(output[['matched']], handle=handle)
close(handle)

handle = file('output_kaks_cor.tsv', 'w')
cat('All:\n', file=handle)
write_readable_table(kaks_cor[['all']], handle=handle)
cat('\nMatched:\n', file=handle)
write_readable_table(kaks_cor[['matched']], handle=handle)
close(handle)
