#!/usr/bin/env Rscript
#
# Report the squared correlation coefficient between the median single-copy
# core distance among pairs of strains and the pairwise distances among
# sequences in a gene.
#

library(data.table)

argv = commandArgs(trailingOnly=TRUE)

## Read pairwise divergence data and gene summary data
pairwise_div = fread(argv[1])
strains = unique(c(pairwise_div[['strain1']], pairwise_div[['strain2']]))
gene_data = read.csv(argv[2], sep='\t', comment.char='#', header=TRUE, as.is=TRUE,
                     check.names=FALSE)
scc_genes = gene_data[gene_data[, 'n_strains'] == length(strains) &
                      gene_data[, 'n_genes'] == length(strains), 'subset']

## Create a matrix of median pairwise distances for the single copy core genes
scc_median_dist = matrix(nrow=length(strains), ncol=length(strains),
                         dimnames=list(strains, strains), data=NaN)
scc_pairwise = pairwise_div[orthoset %in% scc_genes]
for(i in seq_along(strains)) {
    s1 = strains[i]
    for(j in seq_along(strains)) {
        s2 = strains[j]
        if(j > i) next
        scc_median_dist[i, j] = median(c(scc_pairwise[strain1 == s1 & strain2 == s2][['pairwise_prot_dist']],
                                         scc_pairwise[strain1 == s2 & strain2 == s1][['pairwise_prot_dist']]))
        scc_median_dist[j, i] = scc_median_dist[i, j]
    }
}
diag(scc_median_dist) = 0


cat('gene\tr2\tr2_no_paralogs\n', file=stdout())
for(i in 1:nrow(gene_data)) {
    if(i %% 100 == 0) cat(i, 'genes done\n', file=stderr())
    gene = gene_data[i, 'subset']
    pd = pairwise_div[orthoset == gene]
    gene_dist = pd[['pairwise_prot_dist']]
    para = pd[['strain1']] == pd[['strain2']]
    species_dist = numeric(length(gene_dist)) * NaN
    for(j in seq_along(species_dist)) {
        species_dist[j] = scc_median_dist[pd[j][['strain1']],
                                          pd[j][['strain2']]]
    }
    r2 = cor(species_dist, gene_dist)^2
    r2_nopara = NaN
    if(sum(!para) > 0) {
        r2_nopara = cor(species_dist[!para], gene_dist[!para])^2
    }
    cat(gene, '\t', r2, '\t', r2_nopara, '\n', sep='', file=stdout())
}
