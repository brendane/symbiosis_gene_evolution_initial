#!/usr/bin/env python3
"""
    Given the output of extract_orthosets_from_orthofinder_trees.py,
    extract unaligned files of sequences.

    extract_orthosets_seqs.py [--ignore-missing] [--gene-list] [--orthogroups]
        <orthosets> <fasta> <output directory>

    --ignore-missing Skip missing sequences rather than raising an exception.
    --gene-list      Text file, one entry per line, listing orthosets or
                     orthogroups to include (default = use all).
    --orthogroups    Collect sequences for orthogroups rather than subsets
                     (default = use orthosets).
"""
 
import argparse
import csv
import itertools
import sys

from Bio import SeqIO

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--gene-lists')
parser.add_argument('--ignore-missing', default=False, action='store_true')
parser.add_argument('--orthogroups', default=False, action='store_true')
parser.add_argument('--subsubsets', default=False, action='store_true')
parser.add_argument('orthosets')
parser.add_argument('fasta')
parser.add_argument('output')
args = parser.parse_args()

key = 'subset'
if args.orthogroups:
    key = 'orthogroup'
elif args.subsubsets:
    key = 'subsubset'

seqs = SeqIO.index(args.fasta, 'fasta')

keep = set()
if args.gene_lists is not None:
    with open(args.gene_lists, 'rt') as ih:
        for line in ih:
            keep.add(line.strip())

with open(args.orthosets) as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for o, rows in itertools.groupby(rdr, lambda x: x[key]):
        seqs_written = set()
        if len(keep) > 0 and o not in keep:
            continue
        with open(args.output + '/' + o + '.fasta', 'w') as oh:
            for row in rows:
                for g in row['genes'].split(','):
                    if g in seqs_written:
                        raise Exception('%s found twice in %s' % (g, o))
                    try:
                        SeqIO.write(seqs[g], oh, 'fasta')
                        seqs_written.add(g)
                    except KeyError:
                        if args.ignore_missing:
                            sys.stderr.write('WARNING: %s not found in %s' % (g, args.fasta))
                            continue
                        else:
                            raise Exception('%s not found in %s' % (g, args.fasta))
