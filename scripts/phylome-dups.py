# libraries
import argparse
import os
import numpy as np
import pandas as pd
import ete3
import logging

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

# read in phylogeny and root if needed
def load_phylo(phy_fn, phy_id, do_root=True):

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

	# output from event parsing
	return phy


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

if __name__ == '__main__':
	
	# read phylogeny, find speciation events, create network, do clustering
	phy = load_phylo(phy_fn=phy_fn, phy_id=phy_id, do_root=True)
	
	# retrieve all possible orthogroups (genes descending from duplication node)
	ogs = retrieve_dup_descendants(phy)
	ogs["orthogroup"] = ogprefix + ogs["orthogroup"].astype(str)

	# save
	ogs.to_csv("%s/%s.ortholog_groups_from_dups.csv" % (out_fn,phy_id), sep="\t", index=None, mode="w")

	logging.info("%s Done" % phy_id)

