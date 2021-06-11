# libraries
import argparse
import os
import numpy as np
import pandas as pd
import ete3
import logging
import networkx as nx
from networkx.algorithms import community

# argument parser
arp = argparse.ArgumentParser()

# Add the arguments to the parser
arp.add_argument("-i", "--in", required=True, help="Path to a phylogenetic tree in newick format. Each sequence in the tree must have a prefix indicating the species, separated from gene name with a split character. Default split character is \"_\", see --split for options.", type=str)
arp.add_argument("-o", "--out", required=False, default=None, help="OPTIONAL: String. Path to output folder. Defaults to same directory as input file.", type=str)
arp.add_argument("-p", "--phy",  required=False, default=None, help="OPTIONAL: String. Prefix for output files. Defaults to `basename` of input phylogeny. Default behaviour will never overwrite original files, because it adds suffixes.", type=str)
arp.add_argument("-split", "--split", required=False, default="_", help="OPTIONAL: String to use as species prefix delimiter in gene ids, e.g. \"_\" for gene names formatted as speciesA_geneX. Defaults to \"_\".", type=str)
arp.add_argument("-ogprefix","--ogprefix", required=False, default="OG", help="OPTIONAL: String. Prefix for ortholog clusters. Defaults to \"OG\".", type=str)
arl = vars(arp.parse_args())

# input variables
phy_fn = arl["in"]
out_fn = arl["out"]

# output folder
if out_fn is None:
	out_fn = os.path.dirname(phy_fn)
	
# check if out_fn exists, and create it if it doesn't
if not os.path.exists(out_fn):
    os.makedirs(out_fn)

if arl["phy"] is not None:
	phy_id = arl["phy"]
else:
	phy_id = os.path.basename(phy_fn)

split_ch = arl["split"].replace("\"","")
ogprefix = arl["ogprefix"].replace("\"","")



#########################
####### FUNCTIONS #######
#########################

# logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)-5.5s]\t%(message)s", handlers=[ logging.StreamHandler() ] )

# read in phylogeny and execute event parser to obtain table-like network of 
# orthologous relationships, using the species overlap algorithm
def parse_phylo(phy_fn, phy_id, do_root):

	# load input
	phy = ete3.PhyloTree("%s" % (phy_fn))
	logging.info("%s num nodes = %i" % (phy_id,len(phy)))

	# assign species names to tree
	phy.set_species_naming_function(lambda node: node.name.split(split_ch)[0] )

	# resolve polytomies (randomly)
	phy.resolve_polytomy(recursive=True)
	
	# try to find root if unrooted
	if do_root:
	
		# set outgroup using normal midpoint rooting
		logging.info("%s Midpoint root" % phy_id)
		phy_outgroup = phy.get_midpoint_outgroup()

		# set root
		phy.set_outgroup(phy_outgroup)

	# ignore rooting
	else: 

		pass
		logging.info("%s Skip rooting (assume tree is already rooted)" % phy_id)

	# ladderise phylogeny
	phy.ladderize()

	# parse events
	evs, eva, phy, phy_lis = parse_events(phy=phy, do_allpairs=False)
	clu = clusters_ccs(evs=evs, node_list=phy_lis)

	# output from event parsing
	return evs, eva, phy, phy_lis, clu


# parse phylogenies with ETE to obtain a network-like table defining 
# orthologous relationships, using the species overlap algorithm
def parse_events(phy, do_allpairs, min_support_node=0):

	# list of genes in phylogeny
	phy_lis = phy.get_leaf_names()

	# find evolutionary events (duplications and speciations)
	evev = phy.get_descendant_evol_events(sos_thr=0)

	# speciation events
	evs    = np.empty((len(evev)*len(evev), 5), dtype="object")
	evs[:] = np.nan
	n = 0
	for ev in evev:
		# retrieve list of sps in ingroup and outgroup:
		sps_in = np.unique([ i.split(split_ch)[0] for i in ev.in_seqs ])
		sps_ou = np.unique([ i.split(split_ch)[0] for i in ev.out_seqs ])
		# check if node is a speciation node, or a duplication node where both descendant branches have exactly one species, and this is the same species
		if (ev.etype == "S" or (ev.etype == "D" and len(sps_in) == len(sps_ou) == 1 and sps_in == sps_ou)) and ev.branch_supports[0] >= min_support_node:
			for ii in ev.in_seqs:
				for oi in ev.out_seqs:
					evs[n,0] = ii
					evs[n,1] = oi
					evs[n,2] = ev.branch_supports[0]
					evs[n,3] = ev.etype
					evs[n,4] = ev.sos
					n = n + 1
	evs = pd.DataFrame(evs).dropna()
	evs.columns = ["in_gene","out_gene","branch_support","ev_type","sos"]

	# duplications and speciation events
	if do_allpairs:
		eva    = np.empty((len(evev)*len(evev), 5), dtype="object")
		eva[:] = np.nan
		n = 0
		for ev in evev:
			for ii in ev.in_seqs:
				for oi in ev.out_seqs:
					eva[n,0] = ii
					eva[n,1] = oi
					eva[n,2] = ev.branch_supports[0]
					eva[n,3] = ev.etype
					eva[n,4] = ev.sos
					n = n + 1
		eva = pd.DataFrame(eva).dropna()
		eva.columns = ["in_gene","out_gene","branch_support","ev_type","sos"]
	else:
		eva = []

	return evs, eva, phy, phy_lis



# function to cluster a network-like table of orthologs (from ETE) 
def clusters_lpa(evs, node_list, cluster_label="orthogroup"):

	# clustering: create network
	logging.info("Create network")
	evs_e = evs[["in_gene","out_gene","branch_support"]]
	evs_n = nx.convert_matrix.from_pandas_edgelist(evs_e, source="in_gene", target="out_gene", edge_attr="branch_support")
	evs_n.add_nodes_from(node_list)
	
	# clustering: asynchronous label propagation
	logging.info("Find communities LPA")
	clu_c = community.asyn_lpa_communities(evs_n, seed=11)
	clu_c = { frozenset(c) for c in clu_c }
	logging.info("Find communities LPA num clusters = %i" % len(clu_c))
	clu_c_clu = [ i for i, cluster in enumerate(clu_c) for node in cluster ]
	clu_c_noi = [ node for i, cluster in enumerate(clu_c) for node in cluster ]

	# clustering: save output
	clu = pd.DataFrame( { 
		"node"    :  clu_c_noi,
		cluster_label : clu_c_clu,
	}, columns=["node",cluster_label])
	logging.info("Find communities LPA | num clustered genes = %i" % len(clu))

	return clu



# function to cluster a network-like table of orthologs (from ETE) 
def clusters_ccs(evs, node_list, cluster_label="orthogroup", focus_sps = "Hsap"):

	# clustering: create network
	logging.info("Create network")
	evs_s = evs
	evs_e = evs[["in_gene","out_gene","branch_support"]]
	keep_edges = np.where( np.logical_or([ i.split(split_ch)[0] == focus_sps for i in evs["in_gene"].values ] , [ i.split(split_ch)[0] == focus_sps for i in evs["out_gene"].values ] ))[0]
	evs_e = evs_e.loc[keep_edges,]
	evs_n = nx.convert_matrix.from_pandas_edgelist(evs_e, source="in_gene", target="out_gene", edge_attr="branch_support")
	evs_n.add_nodes_from(node_list)
	
	# clustering: asynchronous label propagation
	logging.info("Find communities connected components")
	clu_c = nx.algorithms.connected_components(evs_n)
	clu_c = { frozenset(c) for c in clu_c }
	logging.info("Find communities connected components num clusters = %i" % len(clu_c))
	clu_c_clu = [ i for i, cluster in enumerate(clu_c) for node in cluster ]
	clu_c_noi = [ node for i, cluster in enumerate(clu_c) for node in cluster ]

	# clustering: save output
	clu = pd.DataFrame( { 
		"node"    :  clu_c_noi,
		cluster_label : clu_c_clu,
	}, columns=["node",cluster_label])
	logging.info("Find communities connected components | num clustered genes = %i" % len(clu))

	return clu




# parse phylogenies with ETE lists of orthologous sequences (descendants from a duplication)
def retrieve_dup_descendants(phy):

	# find evolutionary events (duplications and speciations)
	evev = phy.get_descendant_evol_events(sos_thr=0)

	# speciation events
	evs    = np.empty((len(evev)*len(evev), 5), dtype="object")
	evs[:] = np.nan
	n = 0
	g = 0
	for ev in evev:
		if ev.etype == "D":
			g = g + 1
			for ii in ev.in_seqs:
				evs[n,0] = ii
				evs[n,1] = g
				evs[n,2] = ev.branch_supports[0]
				evs[n,3] = ev.etype
				evs[n,4] = ev.sos
				n = n + 1 
			g = g + 1
			for ii in ev.out_seqs:
				evs[n,0] = ii
				evs[n,1] = g
				evs[n,2] = ev.branch_supports[0]
				evs[n,3] = ev.etype
				evs[n,4] = ev.sos
				n = n + 1 
	evs = pd.DataFrame(evs).dropna()
	evs.columns = ["gene","orthogroup","branch_support","ev_type","sos"]

	return evs

def write_tree(phy, out, evc, attributes, sep="|", do_print=True, cut_gene_names=120):
	
	logging.info("Print tree")

	phy_alter = phy.copy(method="newick-extended")

	for i in phy_alter.get_leaves():
		i_name = i.name
		for attribute in attributes:
			c=evc[evc["node"] == i_name][attribute].values
			if c.size == 0: 
				c="NA"
			else:  
				c=c[0]

			# cut if string is too long
			if cut_gene_names is not None:
				cut_gene_names = int(cut_gene_names)
				if len(c) > cut_gene_names+3:
					c=c[:cut_gene_names] + '...'
					
			i.name = str(i.name) + sep + str(c)
		i.name = str(i.name) + sep

	phy_alter.write(outfile=out)

	# print tree in pdf
	if do_print:
		ts = ete3.TreeStyle()
		ts.show_branch_support = True
		ts.show_leaf_name = False
		ts.complete_branch_lines_when_necessary = False
		ts.scale=120
		phy_alter.render("%s.pdf" % out, tree_style=ts)



if __name__ == '__main__':
	
	# read phylogeny, find speciation events, create network, do clustering
	evs, eva, phy, phy_lis, clu = parse_phylo(phy_fn=phy_fn, phy_id=phy_id, do_root=True)
	
	# retrieve all possible orthogroups (genes descending from duplication node)
	clu["orthogroup"] = ogprefix + clu["orthogroup"].astype(str)

	# save
	clu.to_csv("%s/%s.ortholog_groups_from_ccs.csv" % (out_fn,phy_id), sep="\t", index=None, mode="w")
	print_attributes = ["orthogroup"]
	write_tree(phy=phy, out="%s/%s.ortholog_groups.newick" % (out_fn,phy_id), evc=clu, attributes=print_attributes, sep=" | ", do_print=True, cut_gene_names=120)

	logging.info("%s Done" % phy_id)

