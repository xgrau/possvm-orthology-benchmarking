# libraries
library(stringr)
library(scales)
library(ape)

# input
phy_fn = "results_trees/ANTP.genes.iqtree.treefile"
out_fn = "results_trees/permutations_ref/"
ref_fn = "results_searches/ANTP.seed.diamond.csv"

# load reference sequences
ref = read.table(ref_fn, sep="\t", header = F, col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"), stringsAsFactors=F)
ref = ref[ref$pident == 100,]
ref$refOG = stringr::str_split(ref$qseqid, "_", simplify = T)[,3]
ref$refOG_type = stringr::str_split(ref$qseqid, "_", simplify = T)[,2]
ref$species = stringr::str_split(ref$sseqid, pattern="_", simplify = T)[,1]
ref = ref[!duplicated(ref$sseqid),]
ref = ref[order(ref$refOG),]
ref_seqs = ref$sseqid


fracs = c(0.01,0.02,0.05,0.1,0.2,0.5) # fraction of tips to randomise
reps = 1 # how many different randomisations?
set.seed(11)


# original tree
phy = ape::read.tree(phy_fn)

# for each fraction...
for (fra in fracs) {
  
  nperm = round(length(phy$tip.label[phy$tip.label %in% ref_seqs]) * fra)
  # do this X times  
  for (rep in 1:reps) {

    # init tree that will be permuted
    thy = phy
    
    for (per in 1:nperm) {
      
      # get tips to permute
      tips_to_permute = sample(thy$tip.label[thy$tip.label %in% ref_seqs], 2, replace = F)
      
      # reassign
      new_labels = thy$tip.label
      new_labels[thy$tip.label == tips_to_permute[1]] = tips_to_permute[2]
      new_labels[thy$tip.label == tips_to_permute[2]] = tips_to_permute[1]
      thy$tip.label = new_labels
      
    }
    
    # write permuted tree
    write.tree(thy, file=sprintf("%s/ANTP_permutations_%.2f-r%i.newick", out_fn, fra, rep))
    
  }
  
}
