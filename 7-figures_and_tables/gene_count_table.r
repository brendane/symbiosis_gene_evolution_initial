#!/usr/bin/env Rscript
#
# Make of gene statistics for symbiosis genes and for all genes
#

### Symbiosis gene lists and overlap between them

projdir = '/home/tiffinp/epste051/project/symbiosis_gene_evolution' 

indir = file.path(projdir, 'notes/table')

genes_files = list('Classic Signaling'='symbiosis_signaling_genes.edited.2020-05-07.txt',
                   'Classic Fixation'='symbiosis_fixation_genes.edited.2020-05-07.txt',
                   'GWAS Plant Biomass'=c('gwas_candidates/2020-05-11/plant_biomass.A17.top10.txt',
                                          'gwas_candidates/2020-05-11/plant_biomass.R108.top10.txt'),
                   'GWAS Rhizobial Fitness'=c('gwas_candidates/2020-05-11/rhizobial_fitness.A17.top10.txt',
                                              'gwas_candidates/2020-05-11/rhizobial_fitness.R108.top10.txt'))

orthosets = structure(vector('list', length=length(genes_files)),
                      names=names(genes_files))
for(geneset in names(orthosets)) {
    orthosets[[geneset]] = unlist(sapply(genes_files[[geneset]],
                                         function(fn) {
                                             scan(file.path(indir, fn), what='character', sep='\n')
                                         }, USE.NAMES=FALSE))
}

overlap = matrix(nrow=length(orthosets)+1, ncol=length(orthosets)+1,
                 dimnames=list(c(names(orthosets), 'Total'), c(names(orthosets), 'Total')))
for(geneset in names(orthosets)) {
    overlap[geneset, 'Total'] = length(orthosets[[geneset]])
    overlap['Total', geneset] = length(orthosets[[geneset]])
    for(geneset2 in names(orthosets)) {
        overlap[geneset, geneset2] = sum(orthosets[[geneset]] %in% orthosets[[geneset2]])
    }
}

overlap['Total', 'Total'] = length(unique(unlist(orthosets)))
for(i in 1:nrow(overlap)) {
    for(j in 1:ncol(overlap)) {
        if(nchar(overlap[i, j]) < 4) {
            overlap[i, j] = paste0(paste(rep(' ', 4-nchar(overlap[i, j])),
                                         collapse=''), overlap[i, j])
        }
    }
}

write.table(overlap, 'overlap.tsv', row.names=TRUE, col.names=TRUE,
            quote=TRUE, sep='\t')


### Characteristics of symbiosis genes

gene_data = read.csv(file.path(projdir,
                               'results/gene_comparison/73strains_alpha_complete/2020-05-11/data.tsv'),
                     sep='\t', comment.char='#', header=TRUE, as.is=TRUE, check.names=FALSE)

orthosets[['Genome']] = unique(gene_data[, 'gene'])

gene_table = as.data.frame(matrix(nrow=9, ncol=length(names(orthosets)),
                                  dimnames=list(c('Number of genes',
                                                  'Single copy genes',
                                                  'Mean family size',
                                                  'Median family size',
                                                  'Number found in all genomes',
                                                  'Number found in one genome',
                                                  'Mean number of genomes',
                                                  'Median number of genomes',
                                                  'Mean copy number'),
                                                names(orthosets))),
                           stringsAsFactors=FALSE, check.names=FALSE)
gene_table['Number of genes', ] = round(sapply(orthosets, length), 0)
gene_table['Mean family size', ] = sapply(orthosets,
                                                 function(o) {
                                                     round(mean(gene_data[gene_data[, 'gene'] %in% o, 'n_genes']), 1)
                                                 })
gene_table['Median family size', ] = sapply(orthosets,
                                                 function(o) {
                                                     round(median(gene_data[gene_data[, 'gene'] %in% o, 'n_genes']), 1)
                                                 })
gene_table['Mean number of genomes', ] = sapply(orthosets,
                                                function(o) {
                                                    round(mean(gene_data[gene_data[, 'gene'] %in% o, 'n_strains']), 1)
                                                })
gene_table['Median number of genomes', ] = sapply(orthosets,
                                                function(o) {
                                                    round(median(gene_data[gene_data[, 'gene'] %in% o, 'n_strains']), 1)
                                                })
gene_table['Mean copy number', ] = sapply(orthosets,
                                          function(o) {
                                              round(mean(gene_data[gene_data[, 'gene'] %in% o, 'n_genes'] /
                                                         gene_data[gene_data[, 'gene'] %in% o, 'n_strains']), 2)
                                          })

gene_table['Single copy genes', ] = sapply(orthosets,
                                           function(o) {
                                                     sum(gene_data[gene_data[, 'gene'] %in% o, 'n_genes'] ==
                                                         gene_data[gene_data[, 'gene'] %in% o, 'n_strains'])
                                           })

gene_table['Number found in all genomes', ] = sapply(orthosets,
                                                      function(o) {
                                                          sum(gene_data[gene_data[, 'gene'] %in% o, 'n_strains'] == 27)
                                                      })

gene_table['Number found in one genome', ] = sapply(orthosets,
                                                      function(o) {
                                                          sum(gene_data[gene_data[, 'gene'] %in% o, 'n_strains'] == 1)
                                                      })

write.table(gene_table, 'counts.tsv',
            sep='\t', col.names=TRUE, row.names=TRUE, quote=FALSE)
