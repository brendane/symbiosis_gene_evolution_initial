#!/usr/bin/env python3
"""
    Extract certain pieces of information from a GFF3 file created by NCBI.

    ncbi_gff_processor.py <gff3 file>
"""

import argparse
import collections
import csv
import sys
import urllib.parse as upr

def combine_tags(ts0, ts1):
    ret = collections.defaultdict(lambda: '')
    for k, v in ts0.items():
        ret[k] = v
    for k, v in ts1.items():
        ret[k] = v
    return ret

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('gff')
parser.add_argument('prefix')
args = parser.parse_args()

all_tags = {}

with open(args.gff, 'rt') as gffh:
    rdr = csv.reader(gffh, delimiter='\t')
    for row  in rdr:
        if row[0].startswith('#'):
            continue
        tags = collections.defaultdict(lambda: '')
        for t in row[8].split(';'):
            tags[t.split('=')[0]] = '='.join(t.split('=')[1:])
        if tags['gbkey'] not in ['Gene', 'CDS']:
            continue
        ID = tags['ID']
        parent = tags['Parent']
        if ID in all_tags:
            all_tags[ID] = combine_tags(tags, all_tags[ID])
        elif parent in all_tags:
            all_tags[parent] = combine_tags(tags, all_tags[parent])
        else:
            all_tags[ID] = tags


for tags in all_tags.values():
    lt = tags['locus_tag']
    gene = tags['gene']
    product = upr.unquote(tags['product'])
    partial = False; pseudo = False
    if tags['partial'] == 'true':
        partial = True
    if tags['pseudo'] == 'true':
        pseudo = True
    sys.stdout.writelines([args.prefix + lt, '\t', gene, '\t', product, '\t',
                           str(int(partial)), '\t', str(int(pseudo)), '\n'])
