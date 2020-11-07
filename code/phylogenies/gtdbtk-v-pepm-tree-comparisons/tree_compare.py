#!/usr/bin/python
import pandas as pd
from ete3 import Tree

##### bacteria
genome_bac = Tree("gtdbtk.bac120.user_msa.trimgappyout.shared.tre")
pepm_bac = Tree("bacteria_PEPM.curated.gappyout.genomes.shared.tre")

bacnames = [leaf.name for leaf in genome_bac]
baclen = len(genome_bac)

bac_tree_d = {}
bac_tree_d['pepm'] = pepm_bac
bac_tree_d['self'] = genome_bac
for k in range(0, 10, 1):
	t=Tree()
	t.populate(baclen, bacnames)
	bac_tree_d["rtree_"+str(k)] = t

bac_result_d={}
for tree, obj in bac_tree_d.items():
	result = obj.compare(genome_bac, unrooted='True')
	bac_result_d[tree] = {
						'rf': result['rf'], 
						'max_rf': result['max_rf'],
						'norm_rf': result['norm_rf'],
						'effective_tree_size': result['effective_tree_size'],
						'ref_edges_in_source': result['ref_edges_in_source'],
						'source_edges_in_ref': result['source_edges_in_ref'],
						'source_subtrees': result['source_subtrees']
						}

##### archaea
genome_arc = Tree("gtdbtk.ar122.user_msa.trimgappyout.shared.tre")
pepm_arc = Tree("archaea_PEPM.curated.gappyout.genomes.shared.tre")

arcnames = [leaf.name for leaf in genome_arc]
arclen = len(genome_arc)

arc_tree_d = {}
arc_tree_d['pepm'] = pepm_arc
arc_tree_d['self'] = genome_arc
for k in range(0, 10, 1):
	t=Tree()
	t.populate(arclen, arcnames)
	arc_tree_d["rtree_"+str(k)] = t

arc_result_d={}
for tree, obj in arc_tree_d.items():
	result = obj.compare(genome_arc, unrooted='True')
	arc_result_d[tree] = {
						'rf': result['rf'], 
						'max_rf': result['max_rf'],
						'norm_rf': result['norm_rf'],
						'effective_tree_size': result['effective_tree_size'],
						'ref_edges_in_source': result['ref_edges_in_source'],
						'source_edges_in_ref': result['source_edges_in_ref'],
						'source_subtrees': result['source_subtrees']
						}

##### combine dicts
b = pd.DataFrame.from_dict(bac_result_d, orient='index')
b['group'] = "bacteria"

a = pd.DataFrame.from_dict(arc_result_d, orient='index')
a['group'] = "archaea"

final = pd.concat([b,a])
final.to_csv("tree_comparison.csv", header=True, index_label="tree")