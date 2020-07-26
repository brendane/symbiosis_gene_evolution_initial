#!/usr/bin/env python3
"""
    Trim alignment positions that are all gap.
"""

import sys

from Bio import AlignIO

aln = AlignIO.read(sys.stdin, 'fasta')
#aln = AlignIO.read(sys.argv[1], 'fasta')

drop = set()
for i in range(len(aln[0])):
    a = set(aln[:,i])
    if len(a) == 1 and list(a)[0] == '-':
        drop.add(i)

for rec in aln:
    sys.stdout.write('>' + rec.description + '\n')
    for i, b in enumerate(rec):
        if i not in drop:
            sys.stdout.write(b)
    sys.stdout.write('\n')
