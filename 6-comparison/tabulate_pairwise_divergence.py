#!/usr/bin/env python3
"""
    Make a table of pairwise divergence stats per gene.

    tabulate_pairwise_divergence.py <input file> <genes table> <strain list> [<gene column name>]
"""

import argparse
import collections
import csv
import itertools
import sys

import numpy

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--strains')
parser.add_argument('--div-gene-column', default='orthoset')
parser.add_argument('--table-gene-column', default='subset')
parser.add_argument('--exclude-genes')
parser.add_argument('pairwise')
parser.add_argument('table')
args = parser.parse_args()


## Summary stats to use for calculations and columns to run the summary
## stats on. Standard deviation is set to use the same calculation as R.
stats = {'mean':numpy.nanmean,
         'median':numpy.nanmedian,
         'min':numpy.nanmin,
         'max':numpy.nanmax,
         'sd':lambda l: numpy.nanstd(l, ddof=1),
         'n':lambda l: len(l)-sum(numpy.isnan(l))}
columns = ['ka', 'ks', 'ka_ks', 'pairwise_prot_dist']

gncol = args.div_gene_column
gncol2 = args.table_gene_column

target_strains = set()
with open(args.strains, 'rt') as h:
    for line in h: target_strains.add(line.strip())

exclude = set()
if args.exclude_genes is not None:
    with open(args.exclude_genes, 'rt') as h:
        for line in h: exclude.add(line.strip())

singleton_genes = set()
with open(args.table, 'rt') as h:
    rdr = csv.DictReader(h, delimiter='\t')
    for row in rdr:
        gene = row[gncol2]
        ng = 0
        for g in row['genes'].split(','):
            if g in exclude: continue
            for s in target_strains:
                if g.startswith(s + '_'):
                    ng += 1
        if ng == 1:
            singleton_genes.add(gene)

results = {}
scc_prot_dist = []
with open(args.pairwise, 'rt') as h:
    hh = filter(lambda l: not l.startswith('#'), h)
    rdr = csv.DictReader(hh, delimiter='\t')
    for gene, rows in itertools.groupby(rdr, lambda r: r[gncol]):

        data = collections.defaultdict(list)
        ortho_data = collections.defaultdict(list)
        para_data = collections.defaultdict(list)
        genes = set()
        strains = set()
        pairwise_prot = {}
        for row in rows:
            if row['strain1'] not in target_strains or row['strain2'] not in target_strains:
                continue
            if row['gene1'] in exclude or row['gene2'] in exclude:
                continue
            strains.add(row['strain1']); strains.add(row['strain2']) 
            genes.add(row['gene1']); genes.add(row['gene2']) 
            pairwise_prot[tuple(sorted([row['strain1'], row['strain2']]))] = float(row['pairwise_prot_dist'])
            for k, v in row.items():
                if k in columns:
                    try:
                        v = float(v)
                    except ValueError:
                        v = float('nan')
                data[k].append(v)
                if row['strain1'] == row['strain2']:
                    para_data[k].append(v)
                else:
                    ortho_data[k].append(v)

        results[gene] = {}
        for column in columns:
            for stat, func in stats.items():
                k = column + '.' + stat
                try:
                    results[gene][k + '.all'] = func(data[column])
                except:
                    results[gene][k + '.all'] = float('nan')
                try:
                    results[gene][k + '.ortho'] = func(ortho_data[column])
                except:
                    results[gene][k + '.ortho'] = float('nan')
                try:
                    results[gene][k + '.para'] = func(para_data[column])
                except:
                    results[gene][k + '.para'] = float('nan')

        results[gene]['n_strains'] = len(strains)
        results[gene]['n_genes'] = len(genes)
        if len(strains) == len(target_strains) and len(genes) == len(target_strains):
            scc_prot_dist.append(pairwise_prot) ## NOTE: this only stores one value per pair; will not work for multi-copy genes

_scc_pairwise = collections.defaultdict(list)
for dists in scc_prot_dist:
    for s, d, in dists.items():
        _scc_pairwise[s].append(d)
scc_pairwise = {}
for s, d in _scc_pairwise.items():
    scc_pairwise[s] = numpy.mean(d)
    scc_pairwise[(s[1], s[0])] = scc_pairwise[s]

## Second pass to calculate relative protein distance
with open(args.pairwise, 'rt') as h:
    hh = filter(lambda l: not l.startswith('#'), h)
    rdr = csv.DictReader(hh, delimiter='\t')
    for gene, rows in itertools.groupby(rdr, lambda r: r[gncol]):

        rel_prot_dist = []
        for row in rows:
            if row['strain1'] not in target_strains or row['strain2'] not in target_strains:
                continue
            if row['strain1'] == row['strain2']: continue
            rel_prot_dist.append(float(row['pairwise_prot_dist']) /
                                 scc_pairwise[(row['strain1'], row['strain2'])])
        for stat, func in stats.items():
            k = 'relative_prot_dist.' + stat
            try:
                results[gene][k + '.all'] = func(rel_prot_dist)
            except:
                results[gene][k + '.all'] = float('nan')

output_column_names = sorted(results[list(results.keys())[0]])
sys.stdout.write('gene\t' + '\t'.join(output_column_names) + '\n')
for gene, result in results.items():
    sys.stdout.write(gene + '\t' + '\t'.join(str(result[c]) for c in output_column_names) + '\n')
for gene in singleton_genes:
    if gene in results: continue
    sys.stdout.write(gene)
    for c in output_column_names:
        if '.n.' in c:
            sys.stdout.write('\t0')
        elif c == 'n_genes' or c == 'n_strains':
            sys.stdout.write('\t1')
        else:
            sys.stdout.write('\tnan')
    sys.stdout.write('\n')
