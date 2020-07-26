#!/usr/bin/env python3
"""
    Compare pairwise divergence stats using pairwise data and a list of
    randomly chosen background genes.
"""

import argparse
import csv

import numpy
import pandas


def read_random_genes(fname):
    ret = []
    with open(fname, 'rt') as h:
        for line in h:
            ret.append(line.strip().split('\t'))
    return ret

def write_comparison(data, targets, randoms, column, gncol, handle, log=False):
    x = []
    genes = []
    for g in targets:
        try:
            x.append(float(data[data[gncol] == g][column]))
            genes.append(g)
        except:
            continue
    if log:
        d = numpy.log2([0.0001 + xx for xx in x])
    else:
        d = x
    handle.write('\t'.join(genes) + '\n')
    handle.write('\t'.join(map(str, d)) + '\n')
    for rand_set in randoms:
        x = []
        genes = []
        for g in rand_set:
            try:
                x.append(float(data[data[gncol] == g][column]))
                genes.append(g)
            except:
                continue
        if log:
            rd = numpy.log2([0.0001 + xx for xx in x])
        else:
            rd = x
        handle.write('\t'.join(genes) + '\n')
        handle.write('\t'.join(map(str, rd)) + '\n')


parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--gene-column', default='gene')
parser.add_argument('--output')
parser.add_argument('--targets')
parser.add_argument('--strains')
parser.add_argument('--strain-count-rand')
parser.add_argument('--copy-number-rand')
parser.add_argument('--ngenes-rand')
parser.add_argument('--kaks-rand')
parser.add_argument('--kaks-rand-3')
parser.add_argument('--protein-divergence-rand')
parser.add_argument('--relative-protein-divergence-rand')
parser.add_argument('--delta-rand')
parser.add_argument('--delta-pav-rand')
parser.add_argument('--delta-pav-rand-3')
parser.add_argument('--fitz-purvis-rand')
parser.add_argument('--r2-scc-rand')
parser.add_argument('--dtl-rand')
parser.add_argument('divtable')
args = parser.parse_args()

target_strains = set()
with open(args.strains, 'rt') as h:
    for line in h: target_strains.add(line.strip())

target_genes = set()
with open(args.targets, 'rt') as h:
    for line in h: target_genes.add(line.strip())

div_data = pandas.read_table(args.divtable, comment='#', sep='\t', header=0)
div_data['mean_copy_number'] = div_data['n_genes'] / div_data['n_strains']
gncol = args.gene_column

random_genes = {} # key = file name, value = lists of random genes

if args.strain_count_rand is not None:
    rnd_file = args.strain_count_rand
    if rnd_file not in random_genes:
        random_genes[rnd_file] = read_random_genes(rnd_file)
    with open(args.output + '.strain_count.tsv', 'wt') as h:
        write_comparison(div_data, target_genes, random_genes[rnd_file],
                         'n_strains', gncol, h)

if args.copy_number_rand is not None:
    rnd_file = args.copy_number_rand
    if rnd_file not in random_genes:
        random_genes[rnd_file] = read_random_genes(rnd_file)
    with open(args.output + '.copy_number.tsv', 'wt') as h:
        write_comparison(div_data, target_genes, random_genes[rnd_file],
                         'mean_copy_number', gncol, h)

if args.ngenes_rand is not None:
    rnd_file = args.ngenes_rand
    if rnd_file not in random_genes:
        random_genes[rnd_file] = read_random_genes(rnd_file)
    with open(args.output + '.n_genes.tsv', 'wt') as h:
        write_comparison(div_data, target_genes, random_genes[rnd_file],
                         'n_genes', gncol, h)

if args.kaks_rand is not None:
    rnd_file = args.kaks_rand
    if rnd_file not in random_genes:
        random_genes[rnd_file] = read_random_genes(rnd_file)
    for stat in ['median', 'max']:
        with open(args.output + '.log2_' + stat + '_kaks.tsv', 'wt') as h:
            write_comparison(div_data, target_genes, random_genes[rnd_file],
                             'ka_ks.' + stat + '.all', gncol, h, log=True)
        with open(args.output + '.log2_' + stat + '_kaks.para.tsv', 'wt') as h:
            write_comparison(div_data, target_genes, random_genes[rnd_file],
                             'ka_ks.' + stat + '.para', gncol, h, log=True)

if args.kaks_rand_3 is not None:
    rnd_file = args.kaks_rand_3
    if rnd_file not in random_genes:
        random_genes[rnd_file] = read_random_genes(rnd_file)
    dd = div_data[div_data['ka_ks.n.all'] >= 3]
    for stat in ['median', 'max']:
        with open(args.output + '.log2_' + stat + '_kaks_3.tsv', 'wt') as h:
            write_comparison(dd, target_genes, random_genes[rnd_file],
                             'ka_ks.' + stat + '.all', gncol, h, log=True)

if args.protein_divergence_rand is not None:
    rnd_file = args.protein_divergence_rand
    if rnd_file not in random_genes:
        random_genes[rnd_file] = read_random_genes(rnd_file)
    for stat in ['median', 'max']:
        with open(args.output + '.' + stat + '_absdiv.tsv', 'wt') as h:
            write_comparison(div_data, target_genes, random_genes[rnd_file],
                             'pairwise_prot_dist.' + stat + '.all', gncol, h)
        with open(args.output + '.' + stat + '_absdiv.para.tsv', 'wt') as h:
            write_comparison(div_data, target_genes, random_genes[rnd_file],
                             'pairwise_prot_dist.' + stat + '.para', gncol, h)

if args.relative_protein_divergence_rand is not None:
    rnd_file = args.relative_protein_divergence_rand
    if rnd_file not in random_genes:
        random_genes[rnd_file] = read_random_genes(rnd_file)
    for stat in ['median', 'max']:
        with open(args.output + '.log2_' + stat + '_reldiv.tsv', 'wt') as h:
            write_comparison(div_data, target_genes, random_genes[rnd_file],
                             'relative_prot_dist.' + stat + '.all', gncol, h,
                             log=True)

if args.delta_rand is not None:
    rnd_file = args.delta_rand
    if rnd_file not in random_genes:
        random_genes[rnd_file] = read_random_genes(rnd_file)
    with open(args.output + '.delta.tsv', 'wt') as h:
        write_comparison(div_data, target_genes, random_genes[rnd_file],
                         'delta', gncol, h)

if args.delta_pav_rand is not None:
    rnd_file = args.delta_pav_rand
    if rnd_file not in random_genes:
        random_genes[rnd_file] = read_random_genes(rnd_file)
    with open(args.output + '.delta_pav.tsv', 'wt') as h:
        write_comparison(div_data, target_genes, random_genes[rnd_file],
                         'delta_pav', gncol, h)

if args.delta_pav_rand_3 is not None:
    rnd_file = args.delta_pav_rand
    dd = div_data[div_data['ka_ks.n.all'] >= 3]
    if rnd_file not in random_genes:
        random_genes[rnd_file] = read_random_genes(rnd_file)
    with open(args.output + '.delta_pav_3.tsv', 'wt') as h:
        write_comparison(dd, target_genes, random_genes[rnd_file],
                         'delta_pav', gncol, h)

if args.fitz_purvis_rand is not None:
    rnd_file = args.fitz_purvis_rand
    if rnd_file not in random_genes:
        random_genes[rnd_file] = read_random_genes(rnd_file)
    with open(args.output + '.fitz_purvis.tsv', 'wt') as h:
        write_comparison(div_data, target_genes, random_genes[rnd_file],
                         'fitz_purvis', gncol, h)

if args.r2_scc_rand is not None:
    rnd_file = args.r2_scc_rand
    if rnd_file not in random_genes:
        random_genes[rnd_file] = read_random_genes(rnd_file)
    with open(args.output + '.r2_scc.tsv', 'wt') as h:
        write_comparison(div_data, target_genes, random_genes[rnd_file],
                         'r2', gncol, h)

if args.dtl_rand is not None:
    rnd_file = args.dtl_rand
    if rnd_file not in random_genes:
        random_genes[rnd_file] = read_random_genes(rnd_file)
    with open(args.output + '.dup_rate.tsv', 'wt') as h:
        write_comparison(div_data, target_genes, random_genes[rnd_file],
                         'duplication', gncol, h)
    with open(args.output + '.transfer_rate.tsv', 'wt') as h:
        write_comparison(div_data, target_genes, random_genes[rnd_file],
                         'transfer', gncol, h)
    with open(args.output + '.transfer_minus_dup_rate.tsv', 'wt') as h:
        write_comparison(div_data, target_genes, random_genes[rnd_file],
                         'transfer_minus_dup', gncol, h)
    with open(args.output + '.loss_rate.tsv', 'wt') as h:
        write_comparison(div_data, target_genes, random_genes[rnd_file],
                         'loss', gncol, h)
