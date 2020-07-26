#!/usr/bin/env python3
"""
    Given orthosets extracted from OrthoFinder and a table of annotation
    information, provide a summary of gene names and descriptions.

    combine_orthotable_and_annotations.py <orthotable> <annotation>

    --orthogroups3: Input is Orthogroups.tsv from OrthoFinder v2.3 (or higher)
        instead of an orthoset table.
    --no-align: Don't align the descriptions to make a consensus description
"""

import argparse
import collections
import csv
import math
import os
import re
import subprocess
import sys
import tempfile

def code_for_mafft(d):
    return re.sub('[=<>]', ' ', re.sub('-', '##', re.sub(' ', '@@', d)))

def decode_from_mafft(d):
    return re.sub('##', '-', re.sub('-', ' ', re.sub('@', ' ', d)))

def process_mafft_output(barray):
    lines = decode_from_mafft(barray.decode()).split('\n')
    ret = []
    k = -1
    for line in lines:
        if line.startswith('>'):
            if k >= 0:
                ret[k] = decode_from_mafft(ret[k])
            k += 1
            ret.append('')
        else:
            ret[k] += line
    return ret

def get_consensus_description(descriptions, tempin):
    if len(descriptions) > 200:
        ## Reduce to just the top 3, in rough proportion to
        ## their commonality
        dd = collections.Counter(descriptions)
        descriptions = []
        for i, mc in enumerate(dd.most_common()):
            if i > 3:
                break
            for j in range(math.ceil(mc[1] / 10)):
                descriptions.append(mc[0])
    with open(tempin, 'wt') as oh:
        for i, d in enumerate(descriptions):
            oh.write('>' + str(i) + '\n' + code_for_mafft(d) + '\n')
    p = subprocess.Popen(['mafft', '--text', tempin], stdout=subprocess.PIPE)
    if p.wait() != 0:
        raise Exception('MAFFT alignment failed')
    aln = process_mafft_output(p.communicate()[0])
    consensus = ''
    for i in range(max(map(len, aln))):
        count = collections.Counter()
        for j in range(len(aln)):
            try:
                count.update(aln[j][i])
            except IndexError:
                count.update(' ')
        mc = count.most_common()[0] 
        if mc[1] / float(len(aln)) >= 0.5:
            consensus += mc[0]
        else:
            consensus += ''
    return re.sub(' {2,}', ' ', consensus.strip())

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--orthogroups3', default=False, action='store_true')
parser.add_argument('--no-align', default=False, action='store_true')
parser.add_argument('orthotable')
parser.add_argument('annotation')
parser.add_argument('refs', nargs='?')
args = parser.parse_args()


refs = []
if args.refs is not None:
    refs = args.refs.split(',')

annot = {}
with open(args.annotation, 'rt') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        annot[row[0]] = (row[1], row[2])

sys.stdout.write('orthogroup\tsubset\tgene\tdescription')
for ref in refs:
    sys.stdout.write('\t' + ref)
sys.stdout.write('\tdescriptions\n')
with open(args.orthotable, 'rt') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    tin = tempfile.mkstemp(dir='.')[1]
    for row in rdr:
        if args.orthogroups3:
            OG = row['Orthogroup']
            SS = OG + '.0'
            genes = []
            for s in rdr.fieldnames[1:]:
                for g in row[s].split(', '):
                    if len(g) > 0:
                        genes.append(s + '_' + g)
        else:
            OG = row['orthogroup']
            SS = row['subset']
            genes = row['genes'].split(',')
        names = set()
        descriptions = []
        ref_genes = collections.defaultdict(list)
        for gene in genes:
            try:
                a = annot[gene]
            except KeyError:
                sys.stderr.write('WARNING: skipping %s because it is not found in the annotation file\n' % gene)
                continue
            names.add(a[0])
            descriptions.append(a[1])
            for ref in refs:
                if gene.startswith(ref + '_'):
                    ref_genes[ref].append(gene)
        if len(descriptions) > 1 and not args.no_align:
            consensus_description = get_consensus_description(descriptions, tin)
        elif len(descriptions) > 1 and args.no_align:
            consensus_description = collections.Counter(descriptions).most_common()[0][0]
        else:
            consensus_description = descriptions[0]
        sys.stdout.writelines([OG, '\t', SS, '\t',
                               ','.join(n for n in names if n != '' and n != ' '),
                               '\t', consensus_description])
        for ref in refs:
            sys.stdout.write('\t' + ','.join(ref_genes[ref]))
        sys.stdout.write('\t' + '; '.join(set(descriptions)) + '\n')
    os.unlink(tin)
