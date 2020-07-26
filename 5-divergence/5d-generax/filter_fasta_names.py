#!/usr/bin/env python3
"""
    filter_fasta_names.py [-x] <fasta input> <names>
"""

import argparse
import sys

from Bio import SeqIO

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('-x', default=False, action='store_true')
parser.add_argument('input')
parser.add_argument('names')
args = parser.parse_args()

exclude = args.x

names = set()
with open(args.names, 'rt') as ih:
    for line in ih:
        names.add(line.strip())

for rec in SeqIO.parse(args.input, 'fasta'):
    if rec.id in names:
        if not exclude:
            SeqIO.write(rec, sys.stdout, 'fasta')
    else:
        if exclude:
            SeqIO.write(rec, sys.stdout, 'fasta')
