#!/usr/bin/env python3
"""
    Randomly pick sets of genes from a file based on various
    characteristics. File should be tab-delimited with columns
    labeled "gene" (gene name), n_strains, and n_genes.

    randomize_divergence_orthosets.py [OPTIONS] <input file> <input gene list> <N samples>

    OPTIONS
    --output-all    Output file for samples from entire list
    --output-strain-count Output file for samples matched by number of strains
    --output-count-copy Output file for samples matched by number of
        strains and mean copy number

    For each random sample, one gene is randomly chosen to match each
    gene in the input data. That gene is then excluded from further
    random draws.
"""

import argparse
import csv
import random

import pandas

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output-all')
parser.add_argument('--output-strain-count')
parser.add_argument('--output-copy-count')
parser.add_argument('--copy-tol', default=0., type=float)
parser.add_argument('div')
parser.add_argument('genes')
parser.add_argument('N', type=int)
args = parser.parse_args()

N = args.N

target_genes = set()
with open(args.genes, 'rt') as h:
    for line in h:
        l = line.strip()
        if len(l) > 0: target_genes.add(l)

data = pandas.read_table(args.div, header=0, sep='\t', comment='#')

target_data = data[data.gene.isin(target_genes)]
non_target_data = data[~data.gene.isin(target_genes)]
count = target_data['n_strains']
copy_number = target_data['n_genes']

if args.output_all is not None:
    with open(args.output_all, 'wt') as h:
        for i in range(N):
            sampled = list(non_target_data['gene'].sample(n=len(target_genes), replace=False))
            h.write('\t'.join(sampled) + '\n')

if args.output_strain_count is not None:
    with open(args.output_strain_count, 'wt') as h:
        for i in range(N):
            sampled = set()
            for c in count:
                g = non_target_data[(non_target_data.n_strains == c) &
                                    (~non_target_data.gene.isin(sampled))]['gene'].sample(n=1)
                sampled.add(g.values[0])
            h.write('\t'.join(sampled) + '\n')

copy_hi = 1 + args.copy_tol
copy_lo = 1 / copy_hi
if args.output_copy_count is not None:
    with open(args.output_copy_count, 'wt') as h:
        for i in range(N):
            sampled = set()
            for c, cn in zip(count, copy_number):
                g = non_target_data[(non_target_data.n_genes >= cn * copy_lo) &
                                    (non_target_data.n_genes <= cn * copy_hi) &
                                    (non_target_data.n_strains == c) & 
                                    (~non_target_data.gene.isin(sampled))]['gene'].sample(n=1)
                sampled.add(g.values[0])
            h.write('\t'.join(sampled) + '\n')

