#!/usr/bin/env Rscript
#
# Make table of medians (or means) and t-test p-values for several key
# comparisons. Has a lot of overlap with the gene_count_table.
#
# Also reports some correlations.
#

projdir = '/home/tiffinp/epste051/project/symbiosis_gene_evolution' 

indir = file.path(projdir, 'notes/table')

genes_files = list('Classic Signaling'='symbiosis_signaling_genes.edited.2020-05-07.txt',
                   'Classic Fixation'='symbiosis_fixation_genes.edited.2020-05-07.txt',
                   'Host Benefit'=c('gwas_candidates/2020-05-11/plant_biomass.A17.top10.txt',
                                    'gwas_candidates/2020-05-11/plant_biomass.R108.top10.txt'),
                   'Symbiont Fitness'=c('gwas_candidates/2020-05-11/rhizobial_fitness.A17.top10.txt',
                                        'gwas_candidates/2020-05-11/rhizobial_fitness.R108.top10.txt'))
sym_sets = names(genes_files)

genes = structure(vector('list', length=length(genes_files)),
                  names=names(genes_files))
for(geneset in names(genes)) {
    for(gf in genes_files[[geneset]]) {
        file.copy(file.path(indir, gf), '.', overwrite=TRUE)
    }
    genes[[geneset]] = unique(as.character(unlist(sapply(genes_files[[geneset]],
                                                         function(fn) {
                                                             scan(basename(fn), what='character', sep='\n')
                                                         }, USE.NAMES=FALSE))))
}

file.copy(file.path(projdir, 'results/gene_comparison/73strains_alpha_complete/2020-05-11/data.tsv'), '.',
          overwrite=TRUE)
gene_data = read.csv('data.tsv', sep='\t', comment.char='#', header=TRUE, as.is=TRUE, check.names=FALSE)
gene_data[, 'ka_ks.n.all'] = ifelse(is.na(gene_data[, 'ka_ks.n.all']), 0, gene_data[, 'ka_ks.n.all'])

genes[['Genome']] = unique(gene_data[, 'gene'])
genes[['Non-Symbiosis']] = genes[['Genome']][!(genes[['Genome']] %in% unlist(genes[sym_sets]))]

output = matrix(ncol=14, nrow=7, data=NaN,
                dimnames=list(c('Non-Symbiosis', 'Classic Signaling', 'Classic Fixation', 't_test_sig_fix',
                                'Host Benefit', 'Symbiont Fitness', 't_test_ben_fit'),
                              c('genomes', 'copies_per_genome', 'duplication', 'transfer',
                                'median_ka_ks', 'mean_ka_ks', 'median_ka_ks_3', 'mean_ka_ks_3',
                                'r2', 'delta', 'delta_pav', 'delta_pav_3', 'loss', 'n_gene')))

output_n = matrix(ncol=2, nrow=7, data=NaN,
                dimnames=list(c('Non-Symbiosis', 'Classic Signaling', 'Classic Fixation', 't_test_sig_fix',
                                'Host Benefit', 'Symbiont Fitness', 't_test_ben_fit'),
                              c( 'median_ka_ks_3', 'delta_pav_3')))

output_mean = matrix(ncol=14, nrow=5, data=NaN,
                     dimnames=list(c('Non-Symbiosis', 'Classic Signaling', 'Classic Fixation',
                                     'Host Benefit', 'Symbiont Fitness'),
                                   c('genomes', 'copies_per_genome', 'duplication', 'transfer',
                                     'median_ka_ks', 'mean_ka_ks', 'median_ka_ks_3', 'mean_ka_ks_3',
                                     'r2', 'delta', 'delta_pav', 'delta_pav_3', 'loss', 'n_gene')))

# P-values for test for difference with Non-Symbiosis
output_p = matrix(ncol=15, nrow=10, data=NaN,
                     dimnames=list(c('Classic Signaling', 'Classic Fixation',
                                     'Host Benefit', 'Symbiont Fitness',
                                     'fix_vs_sig', 'ben_vs_fit', 
                                     'sig_vs_ben', 'sig_vs_fit', 'fix_vs_ben', 'fix_vs_fit'),
                                   c('genomes', 'copies_per_genome', 'duplication', 'transfer',
                                     'transfer_minus_dup',
                                     'median_ka_ks', 'mean_ka_ks', 'median_ka_ks_3', 'mean_ka_ks_3',
                                     'r2', 'delta', 'delta_pav', 'delta_pav_3', 'loss', 'n_gene')))

# Correlations with median pairwise Ka/Ks
output_kaks_cor = matrix(ncol=12, nrow=5, data=NaN,
                         dimnames=list(c('Non-Symbiosis', 'Classic Signaling', 'Classic Fixation',
                                         'Host Benefit', 'Symbiont Fitness'),
                                       c('transfer', 'duplication', 'delta', 'delta_pav', 'delta_pav_3', 'r2',
                                         'transfer_df', 'duplication_df', 'delta_df', 'delta_pav_df', 'delta_pav_3_df', 'r2_df')))

## Median and mean stats for each measurement
for(ss in c('Non-Symbiosis', sym_sets)) {
    d = gene_data[gene_data[, 'gene'] %in% genes[[ss]], ]
    output[ss, 'genomes'] = median(d[, 'n_strains'])
    output[ss, 'copies_per_genome'] = median(d[, 'n_genes'] / d[, 'n_strains'])
    output[ss, 'duplication'] = median(d[, 'duplication'], na.rm=TRUE)
    output[ss, 'transfer'] = median(d[, 'transfer'], na.rm=TRUE)
    output[ss, 'median_ka_ks'] = median(d[, 'ka_ks.median.all'], na.rm=TRUE)
    output[ss, 'mean_ka_ks'] = median(d[, 'ka_ks.mean.all'], na.rm=TRUE)
    output[ss, 'median_ka_ks_3'] = median(d[d[, 'ka_ks.n.all'] > 2, 'ka_ks.median.all'], na.rm=TRUE)
    output[ss, 'mean_ka_ks_3'] = median(d[d[, 'ka_ks.n.all'] > 2, 'ka_ks.mean.all'], na.rm=TRUE)
    output[ss, 'r2'] = median(d[, 'r2'], na.rm=TRUE)
    output[ss, 'delta'] = median(d[, 'delta'], na.rm=TRUE)
    output[ss, 'delta_pav'] = median(d[, 'delta_pav'], na.rm=TRUE)
    output[ss, 'delta_pav_3'] = median(d[d[, 'n_genes'] > 2, 'delta_pav'], na.rm=TRUE)
    output[ss, 'loss'] = median(d[, 'loss'], na.rm=TRUE)
    output[ss, 'n_gene'] = median(d[, 'n_genes'], na.rm=TRUE)

    output_n[ss, 'median_ka_ks_3'] = sum(!is.na(d[d[, 'ka_ks.n.all'] > 2, 'ka_ks.median.all']))
    output_n[ss, 'delta_pav_3'] = sum(!is.na(d[d[, 'n_genes'] > 2, 'delta_pav']))
}

for(ss in c('Non-Symbiosis', sym_sets)) {
    d = gene_data[gene_data[, 'gene'] %in% genes[[ss]], ]
    output_mean[ss, 'genomes'] = mean(d[, 'n_strains'])
    output_mean[ss, 'copies_per_genome'] = mean(d[, 'n_genes'] / d[, 'n_strains'])
    output_mean[ss, 'duplication'] = mean(d[, 'duplication'], na.rm=TRUE)
    output_mean[ss, 'transfer'] = mean(d[, 'transfer'], na.rm=TRUE)
    output_mean[ss, 'median_ka_ks'] = mean(d[, 'ka_ks.median.all'], na.rm=TRUE)
    output_mean[ss, 'mean_ka_ks'] = mean(d[, 'ka_ks.mean.all'], na.rm=TRUE)
    output_mean[ss, 'median_ka_ks_3'] = mean(d[d[, 'ka_ks.n.all'] > 2, 'ka_ks.median.all'], na.rm=TRUE)
    output_mean[ss, 'mean_ka_ks_3'] = mean(d[d[, 'ka_ks.n.all'] > 2, 'ka_ks.mean.all'], na.rm=TRUE)
    output_mean[ss, 'r2'] = mean(d[, 'r2'], na.rm=TRUE)
    output_mean[ss, 'delta'] = mean(d[, 'delta'], na.rm=TRUE)
    output_mean[ss, 'delta_pav'] = mean(d[, 'delta_pav'], na.rm=TRUE)
    output_mean[ss, 'delta_pav_3'] = mean(d[d[, 'n_genes'] > 2, 'delta_pav'], na.rm=TRUE)
    output_mean[ss, 'loss'] = mean(d[, 'loss'], na.rm=TRUE)
    output_mean[ss, 'n_gene'] = mean(d[, 'n_genes'], na.rm=TRUE)
}

for(ss in rownames(output_kaks_cor)) {
    d = gene_data[gene_data[, 'gene'] %in% genes[[ss]], ]
    d = d[d[, 'ka_ks.n.all'] > 2, ]
    for(st in grep('_df', colnames(output_kaks_cor), value=TRUE, invert=TRUE)) {
        stc = st
        if(st == 'delta_pav_3') stc = 'delta_pav'
        if(sum(!is.na(d[, stc]) & !is.na(d[, 'ka_ks.median.all'])) < 3) {
            cat('Skipping', st, 'correlation because not enough complete observations\n')
        } else {
            ct = cor.test(d[, stc], d[, 'ka_ks.median.all'])
            output_kaks_cor[ss, st] = ct[['estimate']]
            output_kaks_cor[ss, paste0(st, '_df')] = ct[['parameter']]
        }
    }
}


## T-test (unequal variances)
tt = function(x1, x2, ...) {
    t.test(x1, x2, var.equal=FALSE, paired=FALSE, alternative='two.sided', ...)[['p.value']]
}

dns = gene_data[gene_data[, 'gene'] %in% genes[['Non-Symbiosis']], ]
for(ss in sym_sets) {
    d = gene_data[gene_data[, 'gene'] %in% genes[[ss]], ]
    output_p[ss, 'genomes'] = tt(d[, 'n_strains'], dns[, 'n_strains'])
    output_p[ss, 'copies_per_genome'] = tt(d[, 'n_genes'] / d[, 'n_strains'],
                                           dns[, 'n_genes'] / dns[, 'n_strains'])
    output_p[ss, 'duplication'] = tt(d[, 'duplication'], dns[, 'duplication'])
    output_p[ss, 'transfer'] = tt(d[, 'transfer'], dns[, 'transfer'])
    output_p[ss, 'median_ka_ks'] = tt(d[, 'ka_ks.median.all'], dns[, 'ka_ks.median.all'])
    output_p[ss, 'median_ka_ks_3'] = tt(d[d[, 'ka_ks.n.all'] > 2, 'ka_ks.median.all'],
                                        dns[dns[, 'ka_ks.n.all'] > 2, 'ka_ks.median.all'])
    output_p[ss, 'r2'] = tt(d[, 'r2'], dns[, 'r2'])
    output_p[ss, 'delta'] = tt(d[, 'delta'], dns[, 'delta'])
    output_p[ss, 'delta_pav'] = tt(d[, 'delta_pav'], dns[, 'delta_pav'])
    output_p[ss, 'delta_pav_3'] = tt(d[d[, 'n_genes'] > 2, 'delta_pav'], dns[dns[, 'n_genes'] > 2, 'delta_pav'])
    output_p[ss, 'transfer_minus_dup'] = tt(d[, 'transfer_minus_dup'], dns[, 'transfer_minus_dup'])
    output_p[ss, 'loss'] = tt(d[, 'loss'], dns[, 'loss'])
    output_p[ss, 'n_gene'] = tt(d[, 'n_genes'], dns[, 'n_genes'])
}

d1 = gene_data[gene_data[, 'gene'] %in% genes[['Classic Fixation']], ]
d2 = gene_data[gene_data[, 'gene'] %in% genes[['Classic Signaling']], ]
output['t_test_sig_fix', 'genomes'] = tt(d1[, 'n_strains'], d2[, 'n_strains'])
output['t_test_sig_fix', 'copies_per_genome'] = tt(d1[, 'n_genes']/d1[, 'n_strains'],
                                                   d2[, 'n_genes']/d2[, 'n_strains'])
output['t_test_sig_fix', 'duplication'] = tt(d1[, 'duplication'], d2[, 'duplication'])
output['t_test_sig_fix', 'transfer'] = tt(d1[, 'transfer'], d2[, 'transfer'])
output['t_test_sig_fix', 'median_ka_ks'] = tt(d1[, 'ka_ks.median.all'], d2[, 'ka_ks.median.all'])
output['t_test_sig_fix', 'median_ka_ks_3'] = tt(d1[d1[, 'ka_ks.n.all'] > 2, 'ka_ks.median.all'],
                                                d2[d2[, 'ka_ks.n.all'] > 2, 'ka_ks.median.all'])
output['t_test_sig_fix', 'r2'] = tt(d1[, 'r2'], d2[, 'r2'])
output['t_test_sig_fix', 'delta'] = tt(d1[, 'delta'], d2[, 'delta'])
output['t_test_sig_fix', 'delta_pav'] = tt(d1[, 'delta_pav'], d2[, 'delta_pav'])
output['t_test_sig_fix', 'delta_pav_3'] = tt(d1[d[, 'n_genes'] > 2, 'delta_pav'],
                                             d2[d2[, 'n_genes'] > 2, 'delta_pav'])


for(comp in list(c('Classic Fixation', 'Classic Signaling', 'fix_vs_sig'),
                 c('Host Benefit', 'Symbiont Fitness', 'ben_vs_fit'),
                 c('Classic Fixation', 'Host Benefit', 'fix_vs_ben'),
                 c('Classic Fixation', 'Symbiont Fitness', 'fix_vs_fit'),
                 c('Classic Signaling', 'Host Benefit', 'sig_vs_ben'),
                 c('Classic Signaling', 'Symbiont Fitness', 'sig_vs_fit'))) {
    ss = comp[3]
    d1 = gene_data[gene_data[, 'gene'] %in% genes[[comp[1]]], ]
    d2 = gene_data[gene_data[, 'gene'] %in% genes[[comp[2]]], ]
    output_p[ss, 'genomes'] = tt(d1[, 'n_strains'], d2[, 'n_strains'])
    output_p[ss, 'copies_per_genome'] = tt(d1[, 'n_genes'] / d1[, 'n_strains'],
                                           d2[, 'n_genes'] / d2[, 'n_strains'])
    output_p[ss, 'duplication'] = tt(d1[, 'duplication'], d2[, 'duplication'])
    output_p[ss, 'transfer'] = tt(d1[, 'transfer'], d2[, 'transfer'])
    output_p[ss, 'median_ka_ks'] = tt(d1[, 'ka_ks.median.all'], d2[, 'ka_ks.median.all'])
    output_p[ss, 'median_ka_ks_3'] = tt(d1[d1[, 'ka_ks.n.all'] > 2, 'ka_ks.median.all'],
                                        d2[d2[, 'ka_ks.n.all'] > 2, 'ka_ks.median.all'])
    output_p[ss, 'r2'] = tt(d1[, 'r2'], d2[, 'r2'])
    output_p[ss, 'delta'] = tt(d1[, 'delta'], d2[, 'delta'])
    output_p[ss, 'delta_pav'] = tt(d1[, 'delta_pav'], d2[, 'delta_pav'])
    output_p[ss, 'delta_pav_3'] = tt(d1[d1[, 'n_genes'] > 2, 'delta_pav'],
                                   d2[d2[, 'n_genes'] > 2, 'delta_pav'])
    output_p[ss, 'transfer_minus_dup'] = tt(d1[, 'transfer_minus_dup'], d2[, 'transfer_minus_dup'])
    output_p[ss, 'loss'] = tt(d1[, 'loss'], d2[, 'loss'])
    output_p[ss, 'n_gene'] = tt(d1[, 'n_genes'], d2[, 'n_genes'])
}

write_readable_table = function(object, fname) {
    handle = file(fname, 'w')
    for(i in 1:ncol(object)) {
        cat('\n\t', file=handle)
        write.table(object[, i, drop=FALSE], handle, sep='\t',
                    col.names=TRUE, row.names=TRUE, quote=FALSE)
    }
    close(handle)
}

write_readable_table(output, 'medians.tsv')
write_readable_table(output_n, 'n.tsv')
write_readable_table(output_mean, 'means.tsv')
write_readable_table(output_p, 'pvalues_vs_all.tsv')
write_readable_table(output_kaks_cor, 'correlations_with_median_kaks.tsv')
