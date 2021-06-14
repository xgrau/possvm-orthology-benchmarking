# libraries
library(stringr)
library(alluvial)
library(scales)
library(ape)

# input
ref_fn = "refOGs.csv"
out_fn = "results_rooting_inflation"
out_fn = "results_rooting_inflation_one"
# out_fn = "results_rooting_inflation_high"
# out_fn = "results_rooting_inflation_higher"

ort_sets = list(
  list(id="raw", ort_fo="orthobench_trees/raw/")
)

num_iter = 20
edge_fraction = 0.05
# edge_fraction = 0.1  # high
# edge_fraction = 0.2 # higher
inflation_factor_bounds = c(5,20)   # high
# inflation_factor_bounds = c(5,50) # high/higher

# load data
ref = read.table(ref_fn, sep="\t", header = F, col.names = c("refOG", "gene"), stringsAsFactors=F)
fam_list = unique(ref$refOG)

branch_length = data.frame()

for (ort_set in ort_sets) {
  
  # retrieve tree path info for this tree collection
  ort_fo = ort_set$ort_fo
  id = ort_set$id
  
  for (fam in fam_list) {
    
    print(fam)
    
    # load tree    
    tre_fn = sprintf("%s/%s.tre", ort_fo, fam)
    tre = ape::read.tree(tre_fn)
    
    # store original branch lengths
    if (!is.null(tre$edge.length)){
      branch_length = rbind(
        branch_length,
        data.frame(length = tre$edge.length, family = fam)
      )
    }
    
    for (it in 1:num_iter) {
      
      # get tree
      tri = tre
      
      # alter length of a fixed % of edges
      # num_edges_to_alter = floor(length(tri$edge.length) * edge_fraction)
      num_edges_to_alter = 1
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


# distribution of branch lengths
branch_length_nonzero = branch_length [ branch_length$length > 0 , ]
branch_length_nonzero = branch_length
branch_length_nonzero [ branch_length_nonzero == 0, "length" ] = min(branch_length_nonzero$length [ branch_length$length > 0 ])

# plot
pdf("results_evaluation/branch_length_distrubitions.pdf", width = 16, height = 4)
boxplot(length ~ family, data = branch_length_nonzero, log="y", col = "gray90")
# mean and median
brl_mean = mean(branch_length_nonzero$length)
brl_median = median(branch_length_nonzero$length)
abline(h=brl_mean, lty = 2, col = "red")
abline(h=brl_median, lty = 2, col = "purple")
abline(h=brl_mean*5, lty = 2, col = "orange")
abline(h=brl_mean*20, lty = 2, col = "orange")
dev.off()

write.table(branch_length, "results_evaluation/branch_length_distrubitions.csv" , quote = FALSE, row.names = FALSE, sep = "\t")
