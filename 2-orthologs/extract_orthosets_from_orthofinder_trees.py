#!/usr/bin/env python3
"""
    Break orthogroups into sets of orthologs given a tree from OrthoFinder.

    extract_orthosets_from_orthofinder_trees.py --output <output file>
        [--dups <directory with .dups files>]
        <list of strains> <trees file> <Orthogroups.csv> <Orthogroups_UnassignedGenes.csv>
"""

import argparse
import csv
import os
import os.path as osp
import sys

import dendropy

def get_leaf_genes(node, targets=None):
    if node.is_leaf():
        genes = {node.taxon.label}
    else:
        genes = {n.taxon.label for n in node.leaf_nodes()}
    if targets:
        genes = genes.intersection(targets)
    return genes

def get_leaf_strains(node, targets=None):
    return {g.split('_', 1)[0] for g in get_leaf_genes(node, targets)}

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--dups')
parser.add_argument('strainlist')
parser.add_argument('trees')
parser.add_argument('orthogroups')
parser.add_argument('unassigned')
args = parser.parse_args()

## List of strains to include
focal_strains = set()
with open(args.strainlist, 'rt') as ih:
    for line in ih:
        focal_strains.add(line.strip())

dups = {}
for fname in os.listdir(args.dups):
    if fname.endswith('.dups'):
        strain = osp.splitext(fname)[0]
        with open(osp.join(args.dups, fname), 'rt') as ih:
            for line in ih:
                genes = line.strip().split('\t')
                if len(genes) > 1:
                    dups[strain + '_' + genes[0]] = [strain + '_' + g for g in genes[1:]]


## Full list of orthogroups
orthos = {}
with open(args.orthogroups, 'rt') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for row in rdr:
        OG = row['']
        orthos[OG] = set()
        for strain, genes in row.items():
            if strain == '':
                continue
            if genes == '':
                continue
            for g in genes.split(', '):
                orthos[OG].add(strain + '_' + g)
with open(args.unassigned, 'rt') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for row in rdr:
        OG = row['']
        orthos[OG] = set()
        for strain, genes in row.items():
            if strain == '':
                continue
            if genes == '':
                continue
            for g in genes.split(', '):
                orthos[OG].add(strain + '_' + g)

## Loop over trees
OGs_done = set()
with open(args.output, 'wt') as oh:
    oh.write('orthogroup\tsubset\tsingle_copy\tn_strains\tcore\tstrains\tgenes\n')
    with open(args.trees, 'rt') as ih:
        for row in ih:
            OG, tree_string = row.strip().split('\t')
            tree = dendropy.Tree.get(data=tree_string, schema='newick',
                                     rooting='default-rooted',
                                     preserve_underscores=True)

            ## List all strains and genes in this tree
            all_genes = set()
            all_strains = set()
            for l in tree.leaf_node_iter():
                name = l.taxon.label
                all_genes.add(name)
                all_strains.add(name.split('_', 1)[0])

            ## Check that the genes in the tree are the same as the genes
            ## in the orthogroups file.
            if all_genes != orthos[OG]:
                raise Exception('Tree does not match Orthogroups.csv for %s' % OG)

            ## Set of focal strains that are actually in this tree --
            ## these are the strains and genes we want to make sure we
            ## include.
            target_strains = set.intersection(all_strains, focal_strains)
            target_genes = {g for g in all_genes if g.split('_', 1)[0] in focal_strains}

            ## Traverse the tree, recording clades that contain complete sets
            ## of strains and the siblings of those clades.
            potential_orthosets = {}
            node_nums = {}
            for i, node in enumerate(tree.preorder_node_iter()):
                node_nums[str(node)] = i
                if target_strains.issubset(get_leaf_strains(node, target_genes)):
                    potential_orthosets[str(node)] = get_leaf_genes(node, target_genes)
                    for n in node.sibling_nodes():
                        g = get_leaf_genes(n, target_genes)
                        if len(g) > 0:
                            potential_orthosets[str(n)] = g

            ## Eliminate redundant sets
            remove = set()
            for nn1, genes1 in potential_orthosets.items():
                i = node_nums[nn1]
                for nn2, genes2 in potential_orthosets.items():
                    j = node_nums[nn2]
                    if i == j:
                        continue
                    if genes1 == genes2:
                        remove.add(min(i, j))
                    elif genes1.issubset(genes2):
                        remove.add(j)

            ## Write output
            for nn, genes in potential_orthosets.items():
                i = node_nums[nn]
                if i in remove:
                    continue
                gg = []
                strains_in_set = set()
                for g in genes:
                    s = g.split('_', 1)[0]
                    if s in focal_strains:
                        strains_in_set.add(s)
                        gg.append(g)
                        if g in dups:
                            gg += dups[g]
                if len(strains_in_set) == 0:
                    continue

                oh.writelines([OG, '\t', OG + '.' + str(i), '\t',
                               str(int(len(gg) == len(strains_in_set))), '\t',
                               str(len(strains_in_set)), '\t',
                               str(int(len(strains_in_set) == len(focal_strains))), '\t',
                               ','.join(sorted(strains_in_set)), '\t',
                               ','.join(sorted(gg)), '\n'])
            OGs_done.add(OG)

    ## Write information for remaining orthogroups -- should be the
    ## singleton "groups" only.
    for OG, genes in orthos.items():
        i = 0
        if OG in OGs_done:
            continue
        gg = []
        strains_in_set = set()
        for g in genes:
            s = g.split('_', 1)[0]
            if s in focal_strains:
                strains_in_set.add(s)
                gg.append(g)
                if g in dups:
                    gg += dups[g]
        if len(strains_in_set) == 0:
            continue
        oh.writelines([OG, '\t', OG + '.' + str(i), '\t',
                       str(int(len(gg) == len(strains_in_set))), '\t',
                       str(len(strains_in_set)), '\t',
                       str(int(len(strains_in_set) == len(focal_strains))), '\t',
                       ','.join(sorted(strains_in_set)), '\t',
                       ','.join(sorted(gg)), '\n'])
