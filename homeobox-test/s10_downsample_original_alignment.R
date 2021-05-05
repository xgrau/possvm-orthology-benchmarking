# libraries
library(stringr)
library(scales)
library(ape)


# input
fas_fn = "results_trees/ANTP.genes.lt.fasta"
out_fn = "results_trees/downsampling_alignment/"

niter = 20
fracs = c(0.7)

# read data
fas = ape::read.FASTA(fas_fn, type = "AA")
fas_genes = names(fas)
fas_sps = stringr::str_split(fas_genes, pattern = "_", simplify = T)[,1]
sps_lis = unique(fas_sps)
# always exclude Hsap, to keep at least one reference
sps_lis = sps_lis[sps_lis != "Hsap"]

for (frac in fracs) {
  
  nsps = round(frac*length(sps_lis))
  
  for (it in 1:niter) {
    
    print(sprintf("remove %.2f percent (%i) species from %s", frac*100,nsps, fas_fn))
    # species to remove
    sps_remove = sample(sps_lis, nsps, replace = F)
    fas_it = fas[!fas_sps %in% sps_remove]
    write.FASTA(fas_it, file = sprintf("%s/downsample_%.2f-r%i.genes.lt.fasta", out_fn,frac,it))
    
  }
  
}
