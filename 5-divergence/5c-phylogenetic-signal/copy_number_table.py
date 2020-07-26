#!/usr/bin/env python3
"""
    Make a table of strain copy number given a tsv file with a column
    containing comma-separated lists of items. Strain name is
    separated by the first underscore

    copy_number_table.py <input file> <column> <gene_column> [<gene list>]
"""

import argparse
import collections
import csv
import gzip
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('input')
parser.add_argument('column')
parser.add_argument('gene_column')
parser.add_argument('genes', nargs='*', default=[])
args = parser.parse_args()

keep = set()
if len(args.genes):
    with open(args.genes[0], 'rt') as ih:
        for line in ih:
            keep.add(line.strip())

counts = []
strains = set()
op = open
if args.input.endswith('.gz'): op = gzip.open
with op(args.input, 'rt') as ih:
    rdr = csv.DictReader(filter(lambda x: not x.startswith('#'), ih), delimiter='\t')
    for i, row in enumerate(rdr):
        if len(args.genes) > 0 and row[args.gene_column] not in keep:
            continue
        gene_names = row[args.column].split(',')
        strain_count = collections.defaultdict(int)
        for gn in gene_names:
            strain = gn.split('_', 1)[0]
            strain_count[strain] += 1
            strains.add(strain)
        counts.append((row[args.gene_column], strain_count))

s = sorted(strains)
sys.stdout.write('gene\t' + '\t'.join(s) + '\n')
for gene, count in counts:
    sys.stdout.write(gene)
    for strain in s:
        sys.stdout.write('\t' + str(count[strain]))
    sys.stdout.write('\n')
