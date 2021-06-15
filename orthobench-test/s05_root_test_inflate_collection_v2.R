# libraries
library(stringr)
library(alluvial)
library(scales)
library(ape)

# input
ref_fn = "refOGs.csv"
out_fn = "results_rooting_inflated_trees"
# out_fn = "results_rooting_inflation_high"
# out_fn = "results_rooting_inflation_higher"

id="raw"
ort_fo="orthobench_trees/raw/"

num_iter = 20

condition_list = list(
  c(frac = 0.01, quant = 0.95),
  c(frac = 0.05, quant = 0.95),
  c(frac = 0.10, quant = 0.95),
  c(frac = 0.01, quant = 0.99),
  c(frac = 0.05, quant = 0.99),
  c(frac = 0.10, quant = 0.99)
)

# load data
ref = read.table(ref_fn, sep="\t", header = F, col.names = c("refOG", "gene"), stringsAsFactors=F)
fam_list = unique(ref$refOG)


# evaluate branch length quantiles
branch_length = data.frame()
for (fam in fam_list) {
  
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
}

# distribution of branch lengths
branch_length_nonzero = branch_length [ branch_length$length > 0 , ]
branch_length_nonzero = branch_length
branch_length_nonzero [ branch_length_nonzero == 0, "length" ] = min(branch_length_nonzero$length [ branch_length$length > 0 ])

# plot
pdf(sprintf("%s/branch_length_distrubitions.pdf", out_fn), width = 12, height = 4)
boxplot(length ~ family, data = branch_length_nonzero, log="y", col = "gray90")
# mean and median
branch_length_nonzero_quantiles = quantile(branch_length_nonzero$length, c(0.5, 0.95, 0.99, 0.999))
abline(h=branch_length_nonzero_quantiles, lty = 2, col = c("orange", "purple","deeppink","magenta"))
title(sub=sprintf("quantiles %s = %s", paste(names(branch_length_nonzero_quantiles), collapse = ","), paste(format(branch_length_nonzero_quantiles, digits = 2), collapse = ",")))
dev.off()

for (cond in condition_list) {
  
  frac = cond["frac"]
  quant = cond["quant"]
  
  new_length = quantile(branch_length_nonzero$length, quant)
  
  for (fam in fam_list) {
    
    # load tree    
    tre_fn = sprintf("%s/%s.tre", ort_fo, fam)
    tre = ape::read.tree(tre_fn)
    
    # inflate branches  
    for (it in 1:num_iter) {
      
      # get tree
      tri = tre
      
      # alter length of a fixed % of edges
      num_edges_to_alter = floor(length(tri$edge.length) * frac)
      ixs = sample(1:length(tri$edge.length), size = num_edges_to_alter)
      
      for (ix in ixs) {
        tri$edge.length [ ix ] = new_length
      }
      
      # save tree
      ape::write.tree(tri, file = sprintf("%s/%s_it%i_f%.3f_q%.3f.tre", out_fn, fam, it, frac, quant))
      
    }
  }
}
