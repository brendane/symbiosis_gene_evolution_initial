#!/usr/bin/env Rscript
#
# Combine a bunch of pairwise distance information into a couple of tables
# that have a summary for each gene (orthoset).
#
# UPDATE 2020-01-18 added na.rm=TRUE to all mean, median, min, and max calls
#
# UPDATE 2020-05-05: fixed major bug in annotation
#

library(dplyr)
library(optparse)

rcv = function(f, ...) {
    read.csv(f, sep='\t', header=TRUE, comment.char='#', as.is=TRUE, check.names=FALSE, ...)
}


optlist = list(
               make_option('--output'),
               make_option('--exclude',
                           help='text file with one line per sequence to exclude'),
               make_option('--strains',
                           help='text file with one strain per line'),
               make_option('--single-copy-core-distances',
                           help='phylip square matrix with distances between strains to use as a reference'),
               make_option('--kaks',
                           help='tsv file with kaks for all pairs of sequences within each orthoset'),
               make_option('--annotation',
                           help='tsv file with kaks for all pairs of sequences within each orthoset'),
               make_option('--protein-divergence-directory',
                           help='directory with a pairwise protein divergence file for each orthogroup'),
               make_option('--species',
                           help='tsv file with species names; columns named "species" and "strain"')
               )
parser = OptionParser(option_list=optlist)
opts = parse_args(parser, positional_arguments=TRUE)


outpre = opts[['options']][['output']]
ortho_file = opts[['args']][1] # '../results/ortho/73strains/orthofinder/2019-05-02/Results_Directory/reconciliation_May06/orthosets.73strains_alpha_complete.tsv'
osd_file = opts[['args']][2] # '../results/ortho/73strains/orthofinder/2019-05-02/Results_Directory/reconciliation_May06/orthosets.73strains_alpha_complete.dist_0.5.tsv'

strain_file = opts[['options']][['strains']]
spp_file = opts[['options']][['species']] # ('table/73strains_alpha_complete.species.tsv')
exclude_file = opts[['options']][['exclude']] # '../results/gene_aligns/73strains_alpha_complete/muscle/2019-06-09_dist_subsets/pseudogenes.txt'
sc_dist_file = opts[['options']][['single-copy-core-distances']] # '../results/phylo/73strains_alpha_complete/iqtree/2019-05-22/single_copy_core_tree/sc.mldist'
kaks_file = opts[['options']][['kaks']] # '../results/phylo/73strains_alpha_complete/gestimator_dists/2019-07-09_orthosets/dists.tsv'
annot_file = opts[['options']][['annotation']] # '../results/ortho/73strains/orthofinder/2019-05-02/Results_Directory/reconciliation_May06/annotation.73strains_alpha_complete.tsv'
prot_div_dir = opts[['options']][['protein-divergence-directory']] # '../results/phylo/73strains_alpha_complete/fastme/2019-05-31_orthogroups'


## Read in data
ortho = rcv(ortho_file)
ortho = ortho[order(ortho[, 'subset']), ]
osd = rcv(osd_file)

core_strains = scan(strain_file, what='character', sep='\n')

s = rcv(spp_file)
species = structure(s[, 'species'], names=s[, 'strain'])

if(is.null(exclude_file)) {
    xgenes = character()
} else {
    xgenes = scan(exclude_file, what='character', sep='\n')
    xgenes = xgenes[gsub('_.+', '', xgenes) %in% core_strains]
}

sc_trimmed_mldist = read.table(sc_dist_file, skip=1, as.is=TRUE, row.names=1) %>%
    as.matrix()
    colnames(sc_trimmed_mldist) = rownames(sc_trimmed_mldist)

kaks = read.csv(kaks_file, sep='\t', header=TRUE, as.is=TRUE)
kaks = kaks[!(kaks[, 'gene1'] %in% xgenes | kaks[, 'gene2'] %in% xgenes), ]

annot = rcv(annot_file)



## Combine annotations (the gene names and descriptions) across all
## orthosets within an orthogroup, because they probably share functions.
ogs = unique(annot[, 'orthogroup'])
annot_by_og = data.frame('orthogroup'=ogs,
                         'gene'=character(length(ogs)),
                         'descriptions'=character(length(ogs)),
                         stringsAsFactors=FALSE)
gene_annot = aggregate(annot[, 'gene'], list(annot[, 'orthogroup']),
                       function(x) paste(unique(unlist(strsplit(as.character(x), ','))), collapse=','))
descr_annot = aggregate(annot[, 'descriptions'], list(annot[, 'orthogroup']),
                        function(x) paste(unique(unlist(strsplit(as.character(x), ','))), collapse=','))
annot_by_og[, 'gene'] = gene_annot[match(annot_by_og[, 'orthogroup'], gene_annot[, 1]), 2]
annot_by_og[, 'descriptions'] = descr_annot[match(annot_by_og[, 'orthogroup'], descr_annot[, 1]), 2]


## Identify symbiosis genes
ogs = sort(unique(ortho[, 'orthogroup']))
name_pattern = 'nod|noe|nol|nif|fix|nfe|nop'
description_pattern = 'nodul|symbio|fixation|nitrogenase'
description_pattern2 = 'Nod|Noe|Nol|Nif|Fix|Nfe|Nop'
symbiosis_ogs = annot_by_og[grepl(name_pattern, annot_by_og[, 'gene'], ignore.case=TRUE) |
                            grepl(description_pattern2, annot_by_og[, 'descriptions'], ignore.case=FALSE) |
                            grepl(description_pattern, annot_by_og[, 'descriptions'], ignore.case=TRUE), 'orthogroup']
sym_oss = ortho[ortho[, 'orthogroup'] %in% symbiosis_ogs, ]




## Data frames to keep track of orthogroup distance information.
## _scc = distances divided by the concatenated single copy core distance
##        for each pair of strains
## _dups = distances for paralogs (within strain)
## expected_ = distances measured in the concatenated single copy core
##             genes for the strains present in the orthogroup or
##             orthoset; counts each strain once

og_info = data.frame('orthogroup'=ogs,
                     'mean'=numeric(length(ogs)) * NaN,
                     'median'=numeric(length(ogs)) * NaN,
                     'max'=numeric(length(ogs)) * NaN,
                     'min'=numeric(length(ogs)) * NaN,
                     'sd'=numeric(length(ogs)) * NaN,
                     'mean_scc'=numeric(length(ogs)) * NaN,
                     'median_scc'=numeric(length(ogs)) * NaN,
                     'max_scc'=numeric(length(ogs)) * NaN,
                     'min_scc'=numeric(length(ogs)) * NaN,
                     'sd_scc'=numeric(length(ogs)) * NaN,
                     'mean_dups'=numeric(length(ogs)) * NaN,
                     'median_dups'=numeric(length(ogs)) * NaN,
                     'max_dups'=numeric(length(ogs)) * NaN,
                     'min_dups'=numeric(length(ogs)) * NaN,
                     'sd_dups'=numeric(length(ogs)) * NaN,
                     'symbiotic'=numeric(length(ogs)) * NaN,
                     'core'=numeric(length(ogs)) * NaN,
                     'n_strains'=numeric(length(ogs)) * 0,
                     'n_genes'=numeric(length(ogs)) * 0,
                     'expected_mean'=numeric(length(ogs)) * NaN,
                     'expected_median'=numeric(length(ogs)) * NaN,
                     'expected_max'=numeric(length(ogs)) * NaN,
                     'expected_min'=numeric(length(ogs)) * NaN,
                     'expected_mean_copynumber'=numeric(length(ogs)) * NaN,
                     'expected_median_copynumber'=numeric(length(ogs)) * NaN,
                     'expected_max_copynumber'=numeric(length(ogs)) * NaN,
                     'expected_min_copynumber'=numeric(length(ogs)) * NaN,
                     stringsAsFactors=FALSE)

os_info = data.frame('orthogroup'=ortho[, 'orthogroup'],
                     'subset'=ortho[, 'subset'],
                     'mean'=numeric(nrow(ortho)) * NaN,
                     'median'=numeric(nrow(ortho)) * NaN,
                     'max'=numeric(nrow(ortho)) * NaN,
                     'min'=numeric(nrow(ortho)) * NaN,
                     'sd'=numeric(nrow(ortho)) * NaN,
                     'mean_scc'=numeric(nrow(ortho)) * NaN,
                     'median_scc'=numeric(nrow(ortho)) * NaN,
                     'max_scc'=numeric(nrow(ortho)) * NaN,
                     'min_scc'=numeric(nrow(ortho)) * NaN,
                     'sd_scc'=numeric(nrow(ortho)) * NaN,
                     'mean_dups'=numeric(nrow(ortho)) * NaN,
                     'median_dups'=numeric(nrow(ortho)) * NaN,
                     'max_dups'=numeric(nrow(ortho)) * NaN,
                     'min_dups'=numeric(nrow(ortho)) * NaN,
                     'sd_dups'=numeric(nrow(ortho)) * NaN,
                     'n_kaks'=numeric(nrow(ortho)) * 0,
                     'mean_kaks'=numeric(nrow(ortho)) * NaN,
                     'median_kaks'=numeric(nrow(ortho)) * NaN,
                     'max_kaks'=numeric(nrow(ortho)) * NaN,
                     'min_kaks'=numeric(nrow(ortho)) * NaN,
                     'n_dups_kaks'=numeric(nrow(ortho)) * 0,
                     'mean_dups_kaks'=numeric(nrow(ortho)) * NaN,
                     'median_dups_kaks'=numeric(nrow(ortho)) * NaN,
                     'max_dups_kaks'=numeric(nrow(ortho)) * NaN,
                     'min_dups_kaks'=numeric(nrow(ortho)) * NaN,
                     'symbiotic'=numeric(nrow(ortho)) * NaN,
                     'core'=numeric(nrow(ortho)) * NaN,
                     'n_strains'=numeric(nrow(ortho)) * 0,
                     'n_genes'=numeric(nrow(ortho)) * 0,
                     'expected_mean'=numeric(nrow(ortho)) * NaN,
                     'expected_median'=numeric(nrow(ortho)) * NaN,
                     'expected_max'=numeric(nrow(ortho)) * NaN,
                     'expected_min'=numeric(nrow(ortho)) * NaN,
                     'expected_mean_copynumber'=numeric(nrow(ortho)) * NaN,
                     'expected_median_copynumber'=numeric(nrow(ortho)) * NaN,
                     'expected_max_copynumber'=numeric(nrow(ortho)) * NaN,
                     'expected_min_copynumber'=numeric(nrow(ortho)) * NaN,
                     stringsAsFactors=FALSE)

osd_info = data.frame('orthogroup'=osd[, 'orthogroup'],
                      'subset'=osd[, 'subset'],
                      'subsubset'=osd[, 'subsubset'],
                      'mean'=numeric(nrow(osd)) * NaN,
                      'median'=numeric(nrow(osd)) * NaN,
                      'max'=numeric(nrow(osd)) * NaN,
                      'min'=numeric(nrow(osd)) * NaN,
                      'sd'=numeric(nrow(osd)) * NaN,
                      'mean_scc'=numeric(nrow(osd)) * NaN,
                      'median_scc'=numeric(nrow(osd)) * NaN,
                      'max_scc'=numeric(nrow(osd)) * NaN,
                      'min_scc'=numeric(nrow(osd)) * NaN,
                      'sd_scc'=numeric(nrow(osd)) * NaN,
                      'mean_dups'=numeric(nrow(osd)) * NaN,
                      'median_dups'=numeric(nrow(osd)) * NaN,
                      'max_dups'=numeric(nrow(osd)) * NaN,
                      'min_dups'=numeric(nrow(osd)) * NaN,
                      'sd_dups'=numeric(nrow(osd)) * NaN,
                      'n_kaks'=numeric(nrow(osd)) * 0,
                      'mean_kaks'=numeric(nrow(osd)) * NaN,
                      'median_kaks'=numeric(nrow(osd)) * NaN,
                      'max_kaks'=numeric(nrow(osd)) * NaN,
                      'min_kaks'=numeric(nrow(osd)) * NaN,
                      'n_dups_kaks'=numeric(nrow(osd)) * 0,
                      'mean_dups_kaks'=numeric(nrow(osd)) * NaN,
                      'median_dups_kaks'=numeric(nrow(osd)) * NaN,
                      'max_dups_kaks'=numeric(nrow(osd)) * NaN,
                      'min_dups_kaks'=numeric(nrow(osd)) * NaN,
                      'symbiotic'=numeric(nrow(osd)) * NaN,
                      'core'=numeric(nrow(osd)) * NaN,
                      'n_strains'=numeric(nrow(osd)) * 0,
                      'n_genes'=numeric(nrow(osd)) * 0,
                      'expected_mean'=numeric(nrow(osd)) * NaN,
                      'expected_median'=numeric(nrow(osd)) * NaN,
                      'expected_max'=numeric(nrow(osd)) * NaN,
                      'expected_min'=numeric(nrow(osd)) * NaN,
                      'expected_mean_copynumber'=numeric(nrow(osd)) * NaN,
                      'expected_median_copynumber'=numeric(nrow(osd)) * NaN,
                      'expected_max_copynumber'=numeric(nrow(osd)) * NaN,
                      'expected_min_copynumber'=numeric(nrow(osd)) * NaN,
                      stringsAsFactors=FALSE)


os_pairwise_data = data.frame('orthogroup'=gsub('\\..+', '', kaks[, 'orthoset']),
                              'strain1'=gsub('_.+', '', kaks[, 'gene1']),
                              'strain2'=gsub('_.+', '', kaks[, 'gene2']),
                              kaks,
                              'subsubset'=character(nrow(kaks)),
                              'pairwise_prot_dist'=numeric(nrow(kaks)) * NaN,
                              'scc_strain_prot_dist'=numeric(nrow(kaks)) * NaN,
                              'symbiotic'=numeric(nrow(kaks)) * NaN,
                              'core'=numeric(nrow(kaks)) * NaN,
                              'dup'=numeric(nrow(kaks)) * NaN,
                              stringsAsFactors=FALSE)

for(i in seq_along(ogs)) {
    og = ogs[i]
    dists = read.table(paste0(prot_div_dir, '/', og, '.dists'),
                       skip=1, as.is=TRUE, row.names=1) %>%
    as.matrix()
    colnames(dists) = rownames(dists)
    og_genes = unlist(strsplit(ortho[ortho[, 'orthogroup'] == og, 'genes'], ','))
    og_genes = og_genes[!(og_genes %in% xgenes)]
    if(length(og_genes) == 0) next
    dists = dists[og_genes, og_genes, drop=FALSE]
    dist_strains = sapply(colnames(dists), function(x) unlist(strsplit(x, '_'))[1])
    ns = length(unique(dist_strains))
    dists_over_scc = dists / sc_trimmed_mldist[dist_strains, dist_strains]
    dists_over_scc[is.infinite(dists_over_scc)] = NaN
    ut = upper.tri(dists)
    expected_dists = sc_trimmed_mldist[unique(dist_strains), unique(dist_strains)] %>%
    .[upper.tri(.)]
    expected_dists_cn = sc_trimmed_mldist[dist_strains, dist_strains] %>%
    .[upper.tri(.)]
    same_strain = matrix(apply(expand.grid(dist_strains, dist_strains),
                               1, function(x) x[1] == x[2]),
                         nrow=nrow(dists), ncol=nrow(dists),
                         dimnames=list(og_genes, og_genes))
    og_info[i, 'symbiotic'] = as.numeric(og %in% sym_oss[, 'orthogroup'])
    og_info[i, 'n_strains'] = ns
    og_info[i, 'n_genes'] = length(og_genes)
    og_info[i, 'core'] = as.numeric(ns == 27)
    if(length(og_genes) > 1) {
        og_info[i, 'mean'] = mean(dists[ut], na.rm=TRUE)
        og_info[i, 'max'] = max(dists[ut], na.rm=TRUE)
        og_info[i, 'median'] = median(dists[ut], na.rm=TRUE)
        og_info[i, 'min'] = min(dists[ut], na.rm=TRUE)
        og_info[i, 'sd'] = sd(dists[ut], na.rm=TRUE)
        if(ns > 1) {
            og_info[i, 'mean_scc'] = mean(dists_over_scc[ut], na.rm=TRUE)
            og_info[i, 'max_scc'] = max(dists_over_scc[ut], na.rm=TRUE)
            og_info[i, 'min_scc'] = min(dists_over_scc[ut], na.rm=TRUE)
            og_info[i, 'median_scc'] = median(dists_over_scc[ut], na.rm=TRUE)
            og_info[i, 'sd_scc'] = sd(dists_over_scc[ut], na.rm=TRUE)
            og_info[i, 'expected_mean'] = mean(expected_dists, na.rm=TRUE)
            og_info[i, 'expected_median'] = median(expected_dists, na.rm=TRUE)
            og_info[i, 'expected_max'] = max(expected_dists, na.rm=TRUE)
            og_info[i, 'expected_min'] = min(expected_dists, na.rm=TRUE)
        }
        if(sum(same_strain & ut) > 0) {
            og_info[i, 'mean_dups'] = mean(dists[ut & same_strain], na.rm=TRUE)
            og_info[i, 'max_dups'] = max(dists[ut & same_strain], na.rm=TRUE)
            og_info[i, 'min_dups'] = min(dists[ut & same_strain], na.rm=TRUE)
            og_info[i, 'median_dups'] = median(dists[ut & same_strain], na.rm=TRUE)
            og_info[i, 'sd_dups'] = sd(dists[ut & same_strain], na.rm=TRUE)
        }
        og_info[i, 'expected_mean_copynumber'] = mean(expected_dists_cn, na.rm=TRUE)
        og_info[i, 'expected_median_copynumber'] = median(expected_dists_cn, na.rm=TRUE)
        og_info[i, 'expected_max_copynumber'] = max(expected_dists_cn, na.rm=TRUE)
        og_info[i, 'expected_min_copynumber'] = min(expected_dists_cn, na.rm=TRUE)
    }


    oss = ortho[ortho[, 'orthogroup'] == og, 'subset']
    for(os in oss) {
        kki = which(kaks[, 'orthoset'] == os)
        kk = kaks[kki, ]
        s1_kk = gsub('_.+', '', kaks[kki, 'gene1'])
        s2_kk = gsub('_.+', '', kaks[kki, 'gene2'])
        stopifnot(all(os_pairwise_data[kki, 'strain1'] == s1_kk))
        stopifnot(all(os_pairwise_data[kki, 'strain2'] == s2_kk))
        j = which(os_info[, 'subset'] == os)
        os_genes = unlist(strsplit(ortho[j, 'genes'], ','))
        os_genes = os_genes[!(os_genes %in% xgenes)]
        if(length(os_genes) == 0) next
        os_strains = sapply(os_genes, function(x) unlist(strsplit(x, '_'))[1])
        d = dists[os_genes, os_genes, drop=FALSE]
        d_scc = dists_over_scc[os_genes, os_genes, drop=FALSE]
        os_expected_d = sc_trimmed_mldist[unique(os_strains), unique(os_strains)] %>%
        .[upper.tri(.)]
        os_expected_d_cn = sc_trimmed_mldist[os_strains, os_strains] %>%
        .[upper.tri(.)]
        same_strain_os = same_strain[os_genes, os_genes]
        ut = upper.tri(d)
        ns = length(unique(sapply(os_genes, function(x) unlist(strsplit(x, '_'))[1])))
        os_info[j, 'symbiotic'] = as.numeric(os %in% sym_oss[, 'subset'])
        os_info[j, 'n_strains'] = ns
        os_info[j, 'n_genes'] = length(os_genes)
        os_info[j, 'core'] = as.numeric(ns == 27)
        if(length(os_genes) > 1) {
            os_info[j, 'mean'] = mean(d[ut], na.rm=TRUE)
            os_info[j, 'max'] = max(d[ut], na.rm=TRUE)
            os_info[j, 'median'] = median(d[ut], na.rm=TRUE)
            os_info[j, 'min'] = min(d[ut], na.rm=TRUE)
            os_info[j, 'sd'] = sd(d[ut], na.rm=TRUE)
            if(ns > 1) {
                os_info[j, 'mean_scc'] = mean(d_scc[ut], na.rm=TRUE)
                os_info[j, 'max_scc'] = max(d_scc[ut], na.rm=TRUE)
                os_info[j, 'min_scc'] = min(d_scc[ut], na.rm=TRUE)
                os_info[j, 'median_scc'] = median(d_scc[ut], na.rm=TRUE)
                os_info[j, 'sd_scc'] = sd(d_scc[ut], na.rm=TRUE)
                os_info[j, 'expected_mean'] = mean(os_expected_d, na.rm=TRUE)
                os_info[j, 'expected_median'] = median(os_expected_d, na.rm=TRUE)
                os_info[j, 'expected_max'] = max(os_expected_d, na.rm=TRUE)
                os_info[j, 'expected_min'] = min(os_expected_d, na.rm=TRUE)
            }
            if(sum(same_strain_os & ut) > 0) {
                os_info[j, 'mean_dups'] = mean(d[ut & same_strain_os], na.rm=TRUE)
                os_info[j, 'max_dups'] = max(d[ut & same_strain_os], na.rm=TRUE)
                os_info[j, 'min_dups'] = min(d[ut & same_strain_os], na.rm=TRUE)
                os_info[j, 'median_dups'] = median(d[ut & same_strain_os], na.rm=TRUE)
                os_info[j, 'sd_dups'] = sd(d[ut & same_strain_os], na.rm=TRUE)
            }
            os_info[j, 'expected_mean_copynumber'] = mean(os_expected_d_cn, na.rm=TRUE)
            os_info[j, 'expected_median_copynumber'] = median(os_expected_d_cn, na.rm=TRUE)
            os_info[j, 'expected_max_copynumber'] = max(os_expected_d_cn, na.rm=TRUE)
            os_info[j, 'expected_min_copynumber'] = min(os_expected_d_cn, na.rm=TRUE)
            os_info[j, 'n_kaks'] = sum(!is.na(kk[, 'ka_ks']))
            if(os_info[j, 'n_kaks'] > 0) {
                os_info[j, 'mean_kaks'] = mean(kk[, 'ka_ks'], na.rm=TRUE)
                os_info[j, 'median_kaks'] = median(kk[, 'ka_ks'], na.rm=TRUE)
                os_info[j, 'max_kaks'] = max(kk[, 'ka_ks'], na.rm=TRUE)
                os_info[j, 'min_kaks'] = min(kk[, 'ka_ks'], na.rm=TRUE)
                os_info[j, 'n_dups_kaks'] = sum(!is.na(kk[s1_kk == s2_kk, 'ka_ks']))
            } else {
                os_info[j, 'n_dups_kaks'] = 0
            }
            if(os_info[j, 'n_dups_kaks'] > 0) {
                os_info[j, 'mean_dups_kaks'] = mean(kk[s1_kk == s2_kk, 'ka_ks'], na.rm=TRUE)
                os_info[j, 'median_dups_kaks'] = median(kk[s1_kk == s2_kk, 'ka_ks'], na.rm=TRUE)
                os_info[j, 'max_dups_kaks'] = max(kk[s1_kk == s2_kk, 'ka_ks'], na.rm=TRUE)
                os_info[j, 'min_dups_kaks'] = min(kk[s1_kk == s2_kk, 'ka_ks'], na.rm=TRUE)
            }

            os_pairwise_data[kki, 'symbiotic'] = os_info[j, 'symbiotic']
            os_pairwise_data[kki, 'core'] = os_info[j, 'core']
            os_pairwise_data[kki, 'dup'] = as.numeric(s1_kk == s2_kk)
            for(pair_i in 1:nrow(kk)) {
                os_pairwise_data[kki[pair_i], 'pairwise_prot_dist'] = d[kk[pair_i, 'gene1'], kk[pair_i, 'gene2'], drop=TRUE]
                os_pairwise_data[kki[pair_i], 'scc_strain_prot_dist'] = sc_trimmed_mldist[s1_kk[pair_i], s2_kk[pair_i], drop=TRUE]
            }
        }
    }

    ossd = osd[osd[, 'orthogroup'] == og, 'subsubset']
    for(od in ossd) {
        od_genes = unlist(strsplit(osd[osd[, 'subsubset'] == od, 'genes'], ','))
        od_genes = od_genes[!(od_genes %in% xgenes)]
        if(length(od_genes) == 0) next
        kki = which(kaks[, 'gene1'] %in% od_genes | kaks[, 'gene2'] %in% od_genes)
        kk = kaks[kki, ]
        s1_kk = gsub('_.+', '', kaks[kki, 'gene1'])
        s2_kk = gsub('_.+', '', kaks[kki, 'gene2'])
        stopifnot(all(os_pairwise_data[kki, 'strain1'] == s1_kk))
        stopifnot(all(os_pairwise_data[kki, 'strain2'] == s2_kk))
        os_pairwise_data[kki, 'subsubset'] = od

        j = which(osd_info[, 'subsubset'] == od)
        os = osd_info[j, 'subset']
        od_strains = sapply(od_genes, function(x) unlist(strsplit(x, '_'))[1])
        d = dists[od_genes, od_genes, drop=FALSE]
        d_scc = dists_over_scc[od_genes, od_genes, drop=FALSE]
        od_expected_d = sc_trimmed_mldist[unique(od_strains), unique(od_strains)] %>%
        .[upper.tri(.)]
        od_expected_d_cn = sc_trimmed_mldist[od_strains, od_strains] %>%
        .[upper.tri(.)]
        same_strain_od = same_strain[od_genes, od_genes]
        ut = upper.tri(d)
        ns = length(unique(sapply(od_genes, function(x) unlist(strsplit(x, '_'))[1])))
        osd_info[j, 'symbiotic'] = as.numeric(os %in% sym_oss[, 'subset'])
        osd_info[j, 'n_strains'] = ns
        osd_info[j, 'n_genes'] = length(od_genes)
        osd_info[j, 'core'] = as.numeric(ns == 27)
        if(length(od_genes) > 1) {
            osd_info[j, 'mean'] = mean(d[ut], na.rm=TRUE)
            osd_info[j, 'max'] = max(d[ut], na.rm=TRUE)
            osd_info[j, 'median'] = median(d[ut], na.rm=TRUE)
            osd_info[j, 'min'] = min(d[ut], na.rm=TRUE)
            osd_info[j, 'sd'] = sd(d[ut], na.rm=TRUE)
            if(ns > 1) {
                osd_info[j, 'mean_scc'] = mean(d_scc[ut], na.rm=TRUE)
                osd_info[j, 'max_scc'] = max(d_scc[ut], na.rm=TRUE)
                osd_info[j, 'min_scc'] = min(d_scc[ut], na.rm=TRUE)
                osd_info[j, 'median_scc'] = median(d_scc[ut], na.rm=TRUE)
                osd_info[j, 'sd_scc'] = sd(d_scc[ut], na.rm=TRUE)
                osd_info[j, 'expected_mean'] = mean(od_expected_d, na.rm=TRUE)
                osd_info[j, 'expected_median'] = median(od_expected_d, na.rm=TRUE)
                osd_info[j, 'expected_max'] = max(od_expected_d, na.rm=TRUE)
                osd_info[j, 'expected_min'] = min(od_expected_d, na.rm=TRUE)
            }
            if(sum(same_strain_od & ut) > 0) {
                osd_info[j, 'mean_dups'] = mean(d[ut & same_strain_od], na.rm=TRUE)
                osd_info[j, 'max_dups'] = max(d[ut & same_strain_od], na.rm=TRUE)
                osd_info[j, 'min_dups'] = min(d[ut & same_strain_od], na.rm=TRUE)
                osd_info[j, 'median_dups'] = median(d[ut & same_strain_od], na.rm=TRUE)
                osd_info[j, 'sd_dups'] = sd(d[ut & same_strain_od], na.rm=TRUE)
            }
            osd_info[j, 'expected_mean_copynumber'] = mean(od_expected_d_cn, na.rm=TRUE)
            osd_info[j, 'expected_median_copynumber'] = median(od_expected_d_cn, na.rm=TRUE)
            osd_info[j, 'expected_max_copynumber'] = max(od_expected_d_cn, na.rm=TRUE)
            osd_info[j, 'expected_min_copynumber'] = min(od_expected_d_cn, na.rm=TRUE)
            osd_info[j, 'n_kaks'] = sum(!is.na(kk[, 'ka_ks']))
            if(osd_info[j, 'n_kaks'] > 0) {
                osd_info[j, 'mean_kaks'] = mean(kk[, 'ka_ks'], na.rm=TRUE)
                osd_info[j, 'median_kaks'] = median(kk[, 'ka_ks'], na.rm=TRUE)
                osd_info[j, 'max_kaks'] = max(kk[, 'ka_ks'], na.rm=TRUE)
                osd_info[j, 'min_kaks'] = min(kk[, 'ka_ks'], na.rm=TRUE)
                osd_info[j, 'n_dups_kaks'] = sum(!is.na(kk[s1_kk == s2_kk, 'ka_ks']))
            } else {
                osd_info[j, 'n_dups_kaks'] = 0
            }
            if(osd_info[j, 'n_dups_kaks'] > 0) {
                osd_info[j, 'mean_dups_kaks'] = mean(kk[s1_kk == s2_kk, 'ka_ks'], na.rm=TRUE)
                osd_info[j, 'median_dups_kaks'] = median(kk[s1_kk == s2_kk, 'ka_ks'], na.rm=TRUE)
                osd_info[j, 'max_dups_kaks'] = max(kk[s1_kk == s2_kk, 'ka_ks'], na.rm=TRUE)
                osd_info[j, 'min_dups_kaks'] = min(kk[s1_kk == s2_kk, 'ka_ks'], na.rm=TRUE)
            }
        }
    }
    gc()

}

og_info[og_info[, 'n_strains'] == 1, 'expected_mean'] = 0
og_info[og_info[, 'n_strains'] == 1, 'expected_median'] = 0
og_info[og_info[, 'n_strains'] == 1, 'expected_max'] = 0
og_info[og_info[, 'n_strains'] == 1, 'expected_min'] = 0
og_info[og_info[, 'n_strains'] == 1, 'expected_mean_copynumber'] = 0
og_info[og_info[, 'n_strains'] == 1, 'expected_median_copynumber'] = 0
og_info[og_info[, 'n_strains'] == 1, 'expected_max_copynumber'] = 0
og_info[og_info[, 'n_strains'] == 1, 'expected_min_copynumber'] = 0

os_info[os_info[, 'n_strains'] == 1, 'expected_mean'] = 0
os_info[os_info[, 'n_strains'] == 1, 'expected_median'] = 0
os_info[os_info[, 'n_strains'] == 1, 'expected_max'] = 0
os_info[os_info[, 'n_strains'] == 1, 'expected_min'] = 0
os_info[os_info[, 'n_strains'] == 1, 'expected_mean_copynumber'] = 0
os_info[os_info[, 'n_strains'] == 1, 'expected_median_copynumber'] = 0
os_info[os_info[, 'n_strains'] == 1, 'expected_max_copynumber'] = 0
os_info[os_info[, 'n_strains'] == 1, 'expected_min_copynumber'] = 0
os_info[os_info[, 'n_kaks'] < 2, 'mean_kaks'] = NaN
os_info[os_info[, 'n_kaks'] < 2, 'median_kaks'] = NaN
os_info[os_info[, 'n_kaks'] < 2, 'max_kaks'] = NaN
os_info[os_info[, 'n_kaks'] < 2, 'min_kaks'] = NaN
os_info[os_info[, 'n_dups_kaks'] < 2, 'mean_dups_kaks'] = NaN
os_info[os_info[, 'n_dups_kaks'] < 2, 'median__dupskaks'] = NaN
os_info[os_info[, 'n_dups_kaks'] < 2, 'max_dups_kaks'] = NaN
os_info[os_info[, 'n_dups_kaks'] < 2, 'min_dups_kaks'] = NaN

osd_info[osd_info[, 'n_strains'] == 1, 'expected_mean'] = 0
osd_info[osd_info[, 'n_strains'] == 1, 'expected_median'] = 0
osd_info[osd_info[, 'n_strains'] == 1, 'expected_max'] = 0
osd_info[osd_info[, 'n_strains'] == 1, 'expected_min'] = 0
osd_info[osd_info[, 'n_strains'] == 1, 'expected_mean_copynumber'] = 0
osd_info[osd_info[, 'n_strains'] == 1, 'expected_median_copynumber'] = 0
osd_info[osd_info[, 'n_strains'] == 1, 'expected_max_copynumber'] = 0
osd_info[osd_info[, 'n_strains'] == 1, 'expected_min_copynumber'] = 0
osd_info[osd_info[, 'n_kaks'] < 2, 'mean_kaks'] = NaN
osd_info[osd_info[, 'n_kaks'] < 2, 'median_kaks'] = NaN
osd_info[osd_info[, 'n_kaks'] < 2, 'max_kaks'] = NaN
osd_info[osd_info[, 'n_kaks'] < 2, 'min_kaks'] = NaN
osd_info[osd_info[, 'n_dups_kaks'] < 2, 'mean_dups_kaks'] = NaN
osd_info[osd_info[, 'n_dups_kaks'] < 2, 'median__dupskaks'] = NaN
osd_info[osd_info[, 'n_dups_kaks'] < 2, 'max_dups_kaks'] = NaN
osd_info[osd_info[, 'n_dups_kaks'] < 2, 'min_dups_kaks'] = NaN


## Get information on which species are present
osets_spp = sapply(ortho[, 'strains'], function(s) gsub(' ', '_', paste(sort(unique(species[unlist(strsplit(s, ','))])), collapse=',')), USE.NAMES=FALSE)
odsts_spp = sapply(osd[, 'strains'], function(s) gsub(' ', '_', paste(sort(unique(species[unlist(strsplit(s, ','))])), collapse=',')), USE.NAMES=FALSE)
osets_gen = sapply(ortho[, 'strains'], function(s) paste(unique(gsub(' .+', '', sort(unique(species[unlist(strsplit(s, ','))])))), collapse=','), USE.NAMES=FALSE)
odsts_gen = sapply(osd[, 'strains'], function(s) paste(unique(gsub(' .+', '', sort(unique(species[unlist(strsplit(s, ','))])))), collapse=','), USE.NAMES=FALSE)
osets_nspp = sapply(osets_spp, function(x) length(unique(unlist(strsplit(x, ',')))))
odsts_nspp = sapply(odsts_spp, function(x) length(unique(unlist(strsplit(x, ',')))))
osets_ngen = sapply(osets_gen, function(x) length(unique(unlist(strsplit(x, ',')))))
odsts_ngen = sapply(odsts_gen, function(x) length(unique(unlist(strsplit(x, ',')))))

osets_spp_gen = data.frame('subset'=ortho[, 'subset'], 'species'=osets_spp, 'genera'=osets_gen,
                           stringsAsFactors=FALSE)
odsts_spp_gen = data.frame('subsubset'=osd[, 'subsubset'], 'species'=odsts_spp, 'genera'=odsts_gen,
                           stringsAsFactors=FALSE)


## Merge all the information
os = left_join(os_info,
               ortho[, c('subset', 'single_copy', 'strains', 'genes')]) %>%
left_join(annot_by_og) %>%
left_join(osets_spp_gen)
od = left_join(osd_info,
               osd[, c('subsubset', 'single_copy', 'strains', 'genes')]) %>%
left_join(annot_by_og) %>%
left_join(odsts_spp_gen)


f = file(paste0(outpre, '.orthosets.tsv'), open='w')
cat('# orthogroup: OrthoFinder orthogroup\n',
    '# subset: orthoset determined by phylogenetic reconciliation\n',
    '# mean: mean pairwise evolutionary distance among sequences measured by FastME on full protein alignment\n',
    '# max: max pairwise evolutionary distance among sequences measured by FastME on full protein alignment\n',
    '# min: min pairwise evolutionary distance among sequences measured by FastME on full protein alignment\n',
    '# sd: std. deviation of pairwise evolutionary distance among sequences measured by FastME on full protein alignment\n',
    '# median: median pairwise evolutionary distance among sequences measured by FastME on full protein alignment\n',
    '# mean_scc: mean pairwise evol. dist. divided by pairwise evol. distance in concatenated single copy-core trimmed alignment; duplicated copies (divided by zero) are ignored\n',
    '# max_scc: as above, but max value\n',
    '# min_scc: as above, but min\n',
    '# median_scc: as above, but median\n',
    '# sd_scc: as above, but standard deviation\n',
    '# mean_dups: mean pairwise evol. dist. between paralogs (full alignment, protein distance measured by FastME)\n',
    '# max_dups: as above, but max dist. among paralogs\n',
    '# min_dups: as above, but min dist. among paralogs\n',
    '# median_dups: as above, but median\n',
    '# sd_dups: as above, but stdev\n',
    '# n_kaks: number of sequence pairs for which libsequence gestimator calculated a Ka/Ks value (Comeron 1995 method)\n',
    '# mean_kaks: mean pairwise Ka/Ks calculated from full codon alignment using libsequence gestimator\n',
    '# median_kaks: as above, but median\n',
    '# max_kaks: as above, but max\n',
    '# min_kaks: as above, but min\n',
    '# n_kaks: number of paralogous sequence pairs for which libsequence gestimator calculated a Ka/Ks value (Comeron 1995 method)\n',
    '# mean_dups_kaks: mean pairwise Ka/Ks for paralogs\n',
    '# median_dups_kaks: as above, but median\n',
    '# max_dups_kaks: as above, but max\n',
    '# min_dups_kaks: as above, but min\n',
    '# symbiotic: 1 if symbiosis gene, 0 if not\n',
    '# core: 1 if core gene (within 27 strain group), 0 if not\n',
    '# n_strains: number of strain in this orthoset\n',
    '# n_genes: number of total sequences in this orthoset\n',
    '# expected_mean: mean pairwise evol. dist. among the strains in this orthoset measured by FastME on single-copy-core concatenated protein alignment, counting each strain once\n',
    '# expected_median: as above, but median\n',
    '# expected_max: as above, but maximum\n',
    '# expected_min: as above, but minimum\n',
    '# expected_mean_copynumber: mean pairwise evol. dist. among the strains in this orthoset measured by FastME on single-copy-core concatenated protein alignment, counting each strain by the number of sequences it contributes\n',
    '# expected_median_copynumber: as above, but median\n',
    '# expected_max_copynumber: as above, but maximum\n',
    '# expected_min_copynumber: as above, but minimum\n',
    '# single_copy: all genes present in only one copy per strain\n',
    '# strains: comma-separated list of strains present\n',
    '# genes: comma-separated list of genes present\n',
    '# gene: comma-separated list of gene name annotations\n',
    '# descriptions: comma-separated list of annotations\n',
    '# species: comma-separated list of species present\n',
    '# genera: comma-separated list of genera present\n',
    file=f, sep='')
write.table(os[os[, 'n_genes'] > 0, ],
            file=f, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
close(f)


f = file(paste0(outpre, '.orthosets_dist.tsv'), open='w')
cat('# orthogroup: OrthoFinder orthogroup\n',
    '# subset: orthoset determined by phylogenetic reconciliation\n',
    '# subsubset: orthoset divided by phylogenetic distance threshold\n',
    '# mean: mean pairwise evolutionary distance among sequences measured by FastME on full protein alignment\n',
    '# max: max pairwise evolutionary distance among sequences measured by FastME on full protein alignment\n',
    '# min: min pairwise evolutionary distance among sequences measured by FastME on full protein alignment\n',
    '# median: median pairwise evolutionary distance among sequences measured by FastME on full protein alignment\n',
    '# sd: standard deviation of pairwise evolutionary distance among sequences measured by FastME on full protein alignment\n',
    '# mean_scc: mean pairwise evol. dist. divided by pairwise evol. distance in concatenated single copy-core trimmed alignment; duplicated copies (divided by zero) are ignored\n',
    '# max_scc: as above, but max value\n',
    '# min_scc: as above, but min\n',
    '# median_scc: as above, but median\n',
    '# sd_scc: as above, but standard deviation\n',
    '# mean_dups: mean pairwise evol. dist. between paralogs (full alignment, protein distance measured by FastME)\n',
    '# max_dups: as above, but max dist. among paralogs\n',
    '# min_dups: as above, but min dist. among paralogs\n',
    '# median_dups: as above, but median\n',
    '# sd_dups: as above, but stdev\n',
    '# n_kaks: number of sequence pairs for which libsequence gestimator calculated a Ka/Ks value (Comeron 1995 method)\n',
    '# mean_kaks: mean pairwise Ka/Ks calculated from full codon alignment using libsequence gestimator\n',
    '# median_kaks: as above, but median\n',
    '# max_kaks: as above, but max\n',
    '# min_kaks: as above, but min\n',
    '# n_kaks: number of paralogous sequence pairs for which libsequence gestimator calculated a Ka/Ks value (Comeron 1995 method)\n',
    '# mean_dups_kaks: mean pairwise Ka/Ks for paralogs\n',
    '# median_dups_kaks: as above, but median\n',
    '# max_dups_kaks: as above, but max\n',
    '# min_dups_kaks: as above, but min\n',
    '# symbiotic: 1 if symbiosis gene, 0 if not\n',
    '# core: 1 if core gene (within 27 strain group), 0 if not\n',
    '# n_strains: number of strain in this orthoset\n',
    '# n_genes: number of total sequences in this orthoset\n',
    '# expected_mean: mean pairwise evol. dist. among the strains in this orthoset measured by FastME on single-copy-core concatenated protein alignment, counting each strain once\n',
    '# expected_median: as above, but median\n',
    '# expected_max: as above, but maximum\n',
    '# expected_min: as above, but minimum\n',
    '# expected_mean_copynumber: mean pairwise evol. dist. among the strains in this orthoset measured by FastME on single-copy-core concatenated protein alignment, counting each strain by the number of sequences it contributes\n',
    '# expected_median_copynumber: as above, but median\n',
    '# expected_max_copynumber: as above, but maximum\n',
    '# expected_min_copynumber: as above, but minimum\n',
    '# single_copy: all genes present in only one copy per strain\n',
    '# strains: comma-separated list of strains present\n',
    '# genes: comma-separated list of genes present\n',
    '# gene: comma-separated list of gene name annotations\n',
    '# descriptions: comma-separated list of annotations\n',
    '# species: comma-separated list of species present\n',
    '# genera: comma-separated list of genera present\n',
    file=f, sep='')
write.table(od[od[, 'n_genes'] > 0, ],
            file=f, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
close(f)


write.table(og_info[og_info[, 'n_genes'] > 0, ],
            file=paste0(outpre, '.orthogroups.tsv'), sep='\t',
            col.names=TRUE, row.names=FALSE, quote=FALSE)


os_p_d_handle = file(paste0(outpre, '.orthosets_pairwise_sequence.tsv'), 'w')
cat('# orthogroup: OrthoFinder orthogroup\n',
    '# strain1, strain2: the pair of strains in this row\n',
    '# orthoset: subset based on phylogenetic reconcilation\n',
    '# gene1, gene2: the two genes with divergence measured in this row\n',
    '# ka: non-synonymous nucleotide divergence measured by libsequence gestimator\n',
    '# ks: synonymous nucleotide divergence measured by libsequence gestimator\n',
    '# ka_ks: ka / ks ratio\n',
    '# symbiotic: 1 if symbiosis gene, 0 if not\n',
    '# core: 1 if core gene, 0 if not\n',
    '# dup: 1 if genes are paralogs, 0 if not\n',
    file=os_p_d_handle, sep='')
write.table(os_pairwise_data, , file=os_p_d_handle, sep='\t',
            col.names=TRUE, row.names=FALSE, quote=FALSE)
close(os_p_d_handle)
