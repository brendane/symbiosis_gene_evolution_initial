#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 09:11:11 2017

@author: david

Perform directed 'reconciliation' first and then apply EggNOG method

1 - root gene trees on outgroup: unique one this time
2 - infer orthologues
"""
import os
import os.path as osp
import re
import csv
import glob
import argparse
import operator
import itertools
import numpy as np
import multiprocessing as mp
from collections import Counter, defaultdict

import tree as tree_lib
import resolve, util

def GeneToSpecies_dash(g):
  return g.split("_", 1)[0]
 
OrthoFinderIDs = GeneToSpecies_dash
 
def GeneToSpecies_secondDash(g):
  return "_".join(g.split("_", 2)[:2])
 
def GeneToSpecies_3rdDash(g):
  return "_".join(g.split("_", 3)[:3])
 
def GeneToSpecies_dot(g):
  return g.split(".", 1)[0]
 
def GeneToSpecies_hyphen(g):
  return g.split("-", 1)[0]  
 
class RootMap(object):
    def __init__(self, setA, setB, GeneToSpecies):
        self.setA = setA
        self.setB = setB
        self.GeneToSpecies = GeneToSpecies
        
    def GeneMap(self, gene_name):
        sp = self.GeneToSpecies(gene_name)
        if sp in self.setA: return True 
        elif sp in self.setB: return False
        else: 
            print(gene_name)
            print(sp)
            raise Exception
 
def StoreSpeciesSets(t, GeneMap, tag="sp_"):
    tag_up = tag + "up"
    tag_down = tag + "down"  
    for node in t.traverse('postorder'):
        if node.is_leaf():
            node.add_feature(tag_down, {GeneMap(node.name)})
        elif node.is_root():
            continue
        else:
            node.add_feature(tag_down, set.union(*[ch.__getattribute__(tag_down) for ch in node.get_children()]))
    for node in t.traverse('preorder'):
        if node.is_root():
            node.add_feature(tag_up, set())
        else:
            parent = node.up
            if parent.is_root():
                others = [ch for ch in parent.get_children() if ch != node]
                node.add_feature(tag_up, set.union(*[other.__getattribute__(tag_down) for other in others]))
            else:
                others = [ch for ch in parent.get_children() if ch != node]
                sp_downs = set.union(*[other.__getattribute__(tag_down) for other in others])
                node.add_feature(tag_up, parent.__getattribute__(tag_up).union(sp_downs))
    t.add_feature(tag_down, set.union(*[ch.__getattribute__(tag_down) for ch in t.get_children()]))

def OutgroupIngroupSeparationScore(sp_up, sp_down, sett1, sett2, N_recip, n1, n2):
    f_dup = len(sp_up.intersection(sett1)) * len(sp_up.intersection(sett2)) * len(sp_down.intersection(sett1)) * len(sp_down.intersection(sett2)) * N_recip
    f_a = len(sp_up.intersection(sett1)) * (n2-len(sp_up.intersection(sett2))) * (n1-len(sp_down.intersection(sett1))) * len(sp_down.intersection(sett2)) * N_recip
    f_b = (n1-len(sp_up.intersection(sett1))) * len(sp_up.intersection(sett2)) * len(sp_down.intersection(sett1)) * (n2-len(sp_down.intersection(sett2))) * N_recip
    choice = (f_dup, f_a, f_b)
#    print(choice)
    return max(choice)
 
def GetRoots(tree, species_tree_rooted, GeneToSpecies):
    """
    Allow non-binary gene or species trees.
    (A,B,C) => consider splits A|BC, B|AC, C|AB - this applies to gene and species tree
    If a clean ingroup/outgroup split cannot be found then score root by geometric mean of fraction of expected species actually 
    observed for the two splits
    """
    speciesObserved = set([GeneToSpecies(g) for g in tree.get_leaf_names()])
    if len(speciesObserved) == 1:
        return [next(n for n in tree)] # arbitrary root if all genes are from the same species
    
    # use species tree to find correct outgroup according to what species are present in the gene tree
    n = species_tree_rooted
    children = n.get_children()
    leaves = [set(ch.get_leaf_names()) for ch in children]
    have = [len(l.intersection(speciesObserved)) != 0 for l in leaves]
    while sum(have) < 2:
        n = children[have.index(True)]
        children = n.get_children()
        leaves = [set(ch.get_leaf_names()) for ch in children]
        have = [len(l.intersection(speciesObserved)) != 0 for l in leaves]

    # Get splits to look for
    roots_list = []
    scores_list = []   # the fraction completeness of the two clades
#    roots_set = set()
    for i in xrange(len(leaves)):
        t1 = leaves[i]
        t2 = set.union(*[l for j,l in enumerate(leaves) if j!=i])
        # G - set of species in gene tree
        # First relevant split in species tree is (A,B), such that A \cap G \neq \emptyset and A \cap G \neq \emptyset
        # label all nodes in gene tree according the whether subsets of A, B or both lie below node
        StoreSpeciesSets(tree, GeneToSpecies)   # sets of species
        root_mapper = RootMap(t1, t2, GeneToSpecies)    
        sett1 = set(t1)
        sett2 = set(t2)
        nt1 = float(len(t1))
        nt2 = float(len(t2))
        N_recip = 1./(nt1*nt1*nt2*nt2)
        GeneMap = root_mapper.GeneMap
        StoreSpeciesSets(tree, GeneMap, "inout_") # ingroup/outgroup identification
        # find all possible locations in the gene tree at which the root should be

        T = {True,}
        F = {False,}
        TF = set([True, False])
        for m in tree.traverse('postorder'):
            if m.is_leaf(): 
                if len(m.inout_up) == 1 and m.inout_up != m.inout_down:
                    # this is the unique root
                    return [m]
            else:
                if len(m.inout_up) == 1 and len(m.inout_down) == 1 and m.inout_up != m.inout_down:
                    # this is the unique root
                    return [m]
                nodes = m.get_children() if m.is_root() else [m] + m.get_children()
                clades = [ch.inout_down for ch in nodes] if m.is_root() else ([m.inout_up] + [ch.inout_down for ch in m.get_children()])
                # do we have the situation A | B or (A,B),S?
                if len(nodes) == 3:
                    if all([len(c) == 1 for c in clades]) and T in clades and F in clades:
                        # unique root
                        if clades.count(T) == 1:
                            return [nodes[clades.index(T)]]
                        else:
                            return [nodes[clades.index(F)]]
                    elif T in clades and F in clades:
                        #AB-(A,B) or B-(AB,A)
                        ab = [c == TF for c in clades]
                        i = ab.index(True)
                        roots_list.append(nodes[i])
                        sp_down = nodes[i].sp_down
                        sp_up = nodes[i].sp_up
#                        print(m)
                        scores_list.append(OutgroupIngroupSeparationScore(sp_up, sp_down, sett1, sett2, N_recip, nt1, nt2))
                    elif clades.count(TF) >= 2:  
                        # (A,A,A)-excluded, (A,A,AB)-ignore as want A to be bigest without including B, (A,AB,AB), (AB,AB,AB) 
                        i = 0
                        roots_list.append(nodes[i])
                        sp_down = nodes[i].sp_down
                        sp_up = nodes[i].sp_up
#                        print(m)
                        scores_list.append(OutgroupIngroupSeparationScore(sp_up, sp_down, sett1, sett2, N_recip, nt1, nt2))
                elif T in clades and F in clades:
                    roots_list.append(m)
                    scores_list.append(0)  # last choice
    # If we haven't found a unique root then use the scores for completeness of ingroup/outgroup to root
    if len(roots_list) == 0: 
        return [] # This shouldn't occur
    return [sorted(zip(scores_list, roots_list), reverse=True)[0][1]]
                
def WriteQfO2(orthologues_list_pairs_list, outfilename, qAppend = True):
    """ takes a list where each entry is a pair, (genes1, genes2), which are orthologues of one another
    """
    with open(outfilename, 'ab' if qAppend else 'wb') as outfile:
        for gs1, gs2, _, _ in orthologues_list_pairs_list:
            for g1 in gs1:
                g1 = g1.split()[0].split("_")[-1]
                for g2 in gs2:
                    g2 = g2.split()[0].split("_")[-1]
                    outfile.write("%s\t%s\n" % (g1, g2))

GeneToSpeciesCustom = GeneToSpecies_dash

def split_sp_seq_names(g):
    ## This function added by Brendan
    spname = GeneToSpeciesCustom(g)
    seqname = re.sub(spname + '_', '', g)
    return spname, seqname

    
def GetGeneToSpeciesMap(args):
    ## line added by Brendan:
    return GeneToSpeciesCustom
    GeneToSpecies = GeneToSpecies_dash
    if args.separator and args.separator == "dot":
        GeneToSpecies = GeneToSpecies_dot  
    elif args.separator and args.separator == "second_dash":
        GeneToSpecies = GeneToSpecies_secondDash  
    elif args.separator and args.separator == "3rd_dash":
        GeneToSpecies = GeneToSpecies_3rdDash  
    elif args.separator and args.separator == "hyphen":
        GeneToSpecies = GeneToSpecies_hyphen  
    return GeneToSpecies
  
def OverlapSize(node, GeneToSpecies):  
    descendents = [{GeneToSpecies(l) for l in n.get_leaf_names()} for n in node.get_children()]
    intersection = descendents[0].intersection(descendents[1])
    return len(intersection), intersection, descendents[0], descendents[1]

def ResolveOverlap(overlap, sp0, sp1, ch, tree, neighbours, relOverlapCutoff=4):
    """
    Is an overlap suspicious and if so can ift be resolved by identifying genes that are out of place?
    Args:
        overlap - the species with genes in both clades
        sp0 - the species below ch[0]
        sp1 - the species below ch[1]
        ch - the two child nodes
        tree - the gene tree
        neighbours - dictionary species->neighbours, where neighbours is a list of the sets of species observed at successive topological distances from the species
    Returns:
        qSuccess - has the overlap been resolved
        genes_removed - the out-of-place genes that have been removed so as to resolve the overlap
    
    Implementation:
        - The number of species in the overlap must be a 5th or less of the number of species in each clade
        - for each species with genes in both clades: the genes in one clade must all be more out of place (according to the 
          species tree) than all the gene from that species in the other tree
    """
    oSize = len(overlap)
    if relOverlapCutoff*oSize >= len(sp0) and relOverlapCutoff*oSize >= len(sp1): return False, []
    # The overlap looks suspect, misplaced genes?
    # for each species, we'd need to be able to determine that all genes from A or all genes from B are misplaced
    genes_removed = []
    nA_removed = 0
    nB_removed = 0
    qResolved = True
    for sp in overlap:
        ## Modified by Brendan
        #A = [g for g in ch[0].get_leaf_names() if g.split("_")[0] == sp]
        #B = [g for g in ch[1].get_leaf_names() if g.split("_")[0] == sp]
        A = [g for g in ch[0].get_leaf_names() if split_sp_seq_names(g)[0] == sp]
        B = [g for g in ch[1].get_leaf_names() if split_sp_seq_names(g)[0] == sp]
        A_levels = []
        B_levels = []
        for X, level in zip((A,B),(A_levels, B_levels)):
            for g in X:
                gene_node = tree & g
                r = gene_node.up
                ## Modified by Brendan
                #nextSpecies = set([gg.split("_")[0] for gg in r.get_leaf_names()])
                nextSpecies = set([split_sp_seq_names(gg)[0] for gg in r.get_leaf_names()])
                while len(nextSpecies) == 1:
                    r = r.up
                    ## Modified by Brendan
                    nextSpecies = set([split_sp_seq_names(gg)[0] for gg in r.get_leaf_names()])
                    #nextSpecies = set([gg.split("_")[0] for gg in r.get_leaf_names()])
#                print((g, sp,nextSpecies))
                nextSpecies.remove(sp)
                # get the level
                # the sum of the closest and furthest expected distance toplological distance for the closest genes in the gene tree (based on species tree topology)
                neigh = neighbours[sp]
                observed = [neigh[nSp] for nSp in nextSpecies]
                level.append(min(observed) + max(observed))
        qRemoveA = max(B_levels) < min(A_levels)                           
        qRemoveB = max(A_levels) < min(B_levels)                            
        if qRemoveA and relOverlapCutoff*oSize < len(sp0):
            nA_removed += len(A_levels)
            genes_removed.extend(A)
        elif qRemoveB and relOverlapCutoff*oSize < len(sp1):
            nB_removed += len(B_levels)
            genes_removed.extend(B)
        else:
            qResolved = False
            break
    if qResolved:
        return True, set(genes_removed)
    else:
        return False, set()
          
def GetRoot(tree, species_tree_rooted, GeneToSpecies):
        roots = GetRoots(tree, species_tree_rooted, GeneToSpecies)
        if len(roots) > 0:
            root_dists = [r.get_closest_leaf()[1] for r in roots]
            i, _ = max(enumerate(root_dists), key=operator.itemgetter(1))
            return roots[i]
        else:
            return None # single species tree

def GetOrthologues_from_tree(iog, treeFN, species_tree_rooted, GeneToSpecies, neighbours, qWrite=False, dupsWriter=None, treeWriter=None, seqIDs=None, spIDs=None, all_stride_dup_genes=None):
    """ if dupsWriter != None then seqIDs and spIDs must also be provided"""
    qPrune=True
    orthologues = []
    if (not os.path.exists(treeFN)) or os.stat(treeFN).st_size == 0: return set(orthologues), treeFN, set()
    try:
        tree = tree_lib.Tree(treeFN)
    except:
        tree = tree_lib.Tree(treeFN, format=3)
    if len(tree) == 1: return set(orthologues), tree, set()
    root = GetRoot(tree, species_tree_rooted, GeneToSpecies)
    if root == None: return set(orthologues), tree, set()
    # Pick the first root for now
    if root != tree:
        tree.set_outgroup(root)

    tree = Resolve(tree, GeneToSpecies)
    if qPrune: tree.prune(tree.get_leaf_names())
    if len(tree) == 1: return set(orthologues), tree, set()
    """ At this point need to label the tree nodes """
    iNode = 1
    tree.name = "n0"
    suspect_genes = set()
    empty_set = set()
    # preorder traverse so that suspect genes can be identified first, before their closer ortholgoues are proposed
    for n in tree.traverse('preorder'):
        if n.is_leaf(): continue
        if not n.is_root():
            n.name = "n%d" % iNode
            iNode += 1
        ch = n.get_children()
        if len(ch) == 2: 
            oSize, overlap, sp0, sp1 = OverlapSize(n, GeneToSpecies)
            if oSize != 0:
                qResolved, misplaced_genes = ResolveOverlap(overlap, sp0, sp1, ch, tree, neighbours)
            else:
                misplaced_genes = empty_set
            if oSize != 0 and not qResolved:
                if dupsWriter != None:
                    sp_present = sp0.union(sp1)
                    if len(sp_present) == 1:
                        stNode = species_tree_rooted & next(sp for sp in sp_present)
                        isSTRIDE = "Terminal"
                    else:
                        stNode = species_tree_rooted.get_common_ancestor(sp_present)
                        isSTRIDE = "Shared" if all_stride_dup_genes == None else "STRIDE" if frozenset(n.get_leaf_names()) in all_stride_dup_genes else ""
                    dupsWriter.writerow(["OG%07d" % iog, spIDs[stNode.name] if len(stNode) == 1 else stNode.name, n.name, float(oSize)/(len(stNode)), isSTRIDE, ", ".join([seqIDs[g] for g in ch[0].get_leaf_names()]), ", ".join([seqIDs[g] for g in ch[1].get_leaf_names()])]) 
            else:
                # sort out bad genes - no orthology for all the misplaced genes at this level (misplaced_genes). 
                # For previous levels, (suspect_genes) have their orthologues written to suspect orthologues file
                d0 = defaultdict(list)
                d0_sus = defaultdict(list)
                for g in [g for g in ch[0].get_leaf_names() if g not in misplaced_genes]:
                    ## Modified by Brendan (don't know why this was working in pipeline)
                    #sp, seq = g.split("_")
                    sp, seq = split_sp_seq_names(g)
                    if g in suspect_genes:
                        d0_sus[sp].append(seq)
                    else:
                        d0[sp].append(seq)
#                if len(d0_sus) > 0: print(d0_sus)
                d1 = defaultdict(list)
                d1_sus = defaultdict(list)
                for g in [g for g in ch[1].get_leaf_names() if g not in misplaced_genes]:
                    ## Modified by Brendan
                    #sp, seq = g.split("_")
                    sp, seq = split_sp_seq_names(g)
                    if g in suspect_genes:
                        d1_sus[sp].append(seq)
                    else:
                        d1[sp].append(seq)
                orthologues.append((d0, d1, d0_sus, d1_sus))
                suspect_genes.update(misplaced_genes)
        elif len(ch) > 2:
            species = [{GeneToSpecies(l) for l in n.get_leaf_names()} for n in ch]
            for (n0, s0), (n1, s1) in itertools.combinations(zip(ch, species), 2):
                if len(s0.intersection(s1)) == 0:
                    d0 = defaultdict(list)
                    d0_sus = defaultdict(list)
                    for g in n0.get_leaf_names():
                        ## Modified by Brendan
                        #sp, seq = g.split("_")
                        sp, seq = split_sp_seq_names(g)
                        if g in suspect_genes:
                            d0_sus[sp].append(seq)
                        else:
                            d0[sp].append(seq)
                    d1 = defaultdict(list)
                    d1_sus = defaultdict(list)
                    for g in n1.get_leaf_names():
                        ## Modified by Brendan
                        #sp, seq = g.split("_")
                        sp, seq = split_sp_seq_names(g)
                        if g in suspect_genes:
                            d1_sus[sp].append(seq)
                        else:
                            d1[sp].append(seq)
                    orthologues.append((d0, d1, d0_sus, d1_sus))
#    raise Exception("WriteQfO2")
    ## Added by Brendan
    ## Write the rooted, resolved tree to a file
    if treeWriter != None:
        nw = tree.write(features=['name'], outfile=None)
        treeWriter.write("OG%07d\t" % iog); treeWriter.write(nw); treeWriter.write('\n')
    if qWrite:
        directory = os.path.split(treeFN)[0]
        WriteQfO2(orthologues, directory + "/../Orthologues_M3/" + os.path.split(treeFN)[1], qAppend=False)
    return orthologues, tree, suspect_genes

def Resolve(tree, GeneToSpecies):
    StoreSpeciesSets(tree, GeneToSpecies)
    for n in tree.traverse("postorder"):
        tree = resolve.resolve(n, GeneToSpecies)
    return tree

def GetSpeciesNeighbours(t):
    """
    Args: t = rooted species tree
    
    Returns:
    dict: species -> species_dict, such that species_dict: other_species -> toplogical_dist 
    """
    species = t.get_leaf_names()
    levels = {s:[] for s in species}
    for n in t.traverse('postorder'):
        if n.is_leaf(): continue
        children = n.get_children()
        leaf_sets = [set(ch.get_leaf_names()) for ch in children]
        not_i = [set.union(*[l for j, l in enumerate(leaf_sets) if j != i]) for i in xrange(len(children))]
        for l,n in zip(leaf_sets, not_i):
            for ll in l:
                levels[ll].append(n)
    neighbours = {sp:{other:n for n,others in enumerate(lev) for other in others} for sp, lev in levels.items()}
    return neighbours
   
def GetOrthologuesStandalone_Serial(trees_file_list, species_tree_rooted_fn, GeneToSpecies, output_dir,
                                    qSingleTree, dup_file, tree_file, spIDs, seqIDs):
    # Brendan: modified so that the first argument is a list of tree
    # files instead of a directory. Added spIDs and seqIDs arguments.
    species_tree_rooted = tree_lib.Tree(species_tree_rooted_fn)
    neighbours = GetSpeciesNeighbours(species_tree_rooted)
    with open(dup_file, 'w') as dh:
        with open(tree_file, 'w') as treeWriter:
            dupsWriter = csv.writer(dh, delimiter='\t')
            dupsWriter.writerow(["Orthogroup", "Species Tree Node", "Gene Tree Node", "Support", "Type",	"Genes 1", "Genes 2"])
            for treeFn in trees_file_list:
                ##print(treeFn)
                ## Brendan:
                ## added dupsWriter handle, set qWrite to False, added seqIDs and spIDs
                iog = int(re.sub('OG([0-9]+).+', '\\1',
                                 osp.splitext(osp.basename(treeFn))[0]))
                GetOrthologues_from_tree(iog, treeFn, species_tree_rooted, GeneToSpecies, neighbours, False,
                                         dupsWriter, treeWriter, seqIDs=seqIDs, spIDs=spIDs)


if __name__ == "__main__":
    ## Brendan: add argument for files needed to get orthogroups list and 
    ## species/strain lists.
    parser = argparse.ArgumentParser()
    parser.add_argument('--output')
    parser.add_argument("-s", "--separator",
                        choices=("dot", "dash", "second_dash", "3rd_dash", "hyphen"),
                        help="Separator between species name and gene name in gene tree taxa")
    parser.add_argument("trees_dir")
    parser.add_argument("rooted_species_tree")
    parser.add_argument('species')
    parser.add_argument('sequences')
    parser.add_argument('nbatches', type=int)
    parser.add_argument('batch', type=int) # Should start with zero
    args = parser.parse_args()

    output_dir = args.output
    n_batches = args.nbatches
    batch = args.batch

    qSingleTree = False
 
    ## Added by Brendan.
    ## Make a list of all the tree files in this batch. We assign
    # orthogroups to batches in this manner to balance the amount of
    # time needed by each batch -- the largest orthogroups are given
    # the lowest numbers by OrthoFinder.
    tree_files = []
    all_tree_files = os.listdir(args.trees_dir)
    all_tree_files.sort()
    i = 0
    for fname in all_tree_files:
        if '_tree.txt' in fname:
            i += 1
            if i % n_batches - batch == 0:
                tree_files.append(args.trees_dir + '/' + fname)
   
    GeneToSpecies = GeneToSpeciesCustom
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)


    ## Added by Brendan.
    spNums = {}
    spIDs = {}
    with open(args.species, 'r') as ih:
        for line in ih:
            i, name = line.strip().split(': ')
            nm = osp.splitext(name)[0]
            spNums[int(i)] = nm
            spIDs[nm] = nm
    seqIDs = {}
    with open(args.sequences, 'r') as ih:
        for line in ih:
            i, name = line.strip().split(': ')
            spi = int(i.split('_')[0])
            seqIDs[spNums[spi] + '_' + name] = spNums[spi] + '_' + name

    ## Brendan:
    ## Modified version of function that takes a list of files and has an
    ## argument for the output file.
    GetOrthologuesStandalone_Serial(tree_files, args.rooted_species_tree,
                                    GeneToSpecies, output_dir, qSingleTree,
                                    output_dir + '/Duplications.' + str(batch) + '.csv',
                                    output_dir + '/RootedTrees.' + str(batch) + '.tsv',
                                    spIDs, seqIDs)
    util.Success()
