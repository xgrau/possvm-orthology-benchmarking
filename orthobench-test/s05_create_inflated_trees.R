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
inflation_factor = 10

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
      
      # inflate one edge
      tri = tre
      ix = sample(1:length(tri$edge.length), size = 1)
      tri$edge.length [ ix ] = tri$edge.length [ ix ] * inflation_factor
      
      # save tree
      ape::write.tree(tri, file = sprintf("%s/%s_it%i.tre", out_fn, fam, it))
      
    }
    
    
  }
}