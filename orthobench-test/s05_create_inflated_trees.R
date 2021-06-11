# libraries
library(stringr)
library(alluvial)
library(scales)
library(ape)

# input
ref_fn = "refOGs.csv"
out_fn = "results_rooting_inflation"

ort_sets = list(
  list(id="raw", ort_fo="orthobench_trees/raw/")
)

num_iter = 20
edge_fraction = 0.05
inflation_factor_bounds = c(5,20)

# load data
ref = read.table(ref_fn, sep="\t", header = F, col.names = c("refOG", "gene"), stringsAsFactors=F)
fam_list = unique(ref$refOG)


for (ort_set in ort_sets) {
  
  # retrieve tree path info for this tree collection
  ort_fo = ort_set$ort_fo
  id = ort_set$id
  
  for (fam in fam_list) {
    
    # load tree    
    tre_fn = sprintf("%s/%s.tre", ort_fo, fam)
    tre = ape::read.tree(tre_fn)
    
    for (it in 1:num_iter) {
      
      # get tree
      tri = tre
      
      # alter length of a fixed % of edges
      num_edges_to_alter = floor(length(tri$edge.length) * edge_fraction)
      ixs = sample(1:length(tri$edge.length), size = num_edges_to_alter)
      
      for (ix in ixs) {
        inflation_factor = runif(1, min = inflation_factor_bounds[1], max = inflation_factor_bounds[2])
        tri$edge.length [ ix ] = tri$edge.length [ ix ] * inflation_factor
      }
      
      # save tree
      ape::write.tree(tri, file = sprintf("%s/%s_it%i.tre", out_fn, fam, it))
      
    }
    
    
  }
}