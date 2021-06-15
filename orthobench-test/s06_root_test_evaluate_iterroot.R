# libraries
library(stringr)
library(alluvial)
library(scales)
library(ape)
library(phangorn)
library(digest)

# input
ref_fn = "refOGs.csv"
ort_fo = "results_rooting_inflated_trees/"
list_ort_files = list.files(path = ort_fo, pattern = "*groups.csv")

# load data
ref = read.table(ref_fn, sep="\t", header = F, col.names = c("refOG", "gene"), stringsAsFactors=F)
ref$species = stringr::str_split(ref$gene, pattern="_", simplify = T)[,1]
fam_list = unique(ref$refOG)
sps_list = unique(ref$species)

num_iter = 20
# fam_list = fam_list[1:60]


#### all trees ####

# loop

dia = data.frame()
for (ort_fn in list_ort_files) {
  
  ort_info   = stringr::str_split(ort_fn, "_")[[1]]
  fam        = ort_info[1]
  ort_it     = ort_info[2]
  ort_frac   = ort_info[3]
  ort_quant  = gsub("\\.possom","",ort_info[4])
  ort_method = gsub("\\.ortholog","",ort_info[5])
  ort_phy_fn = gsub("\\.csv$",".newick", ort_fn)
  
  # get root hash
  phy_fn = sprintf("%s/%s", ort_fo, ort_phy_fn)
  phy = ape::read.tree(phy_fn)
  phy$node.label = NULL
  root_descendants_ix = phangorn::Descendants(phy, phy$edge[1,2])[[1]]
  root_descendants = phy$tip.label [ root_descendants_ix ]
  root_descendants = paste(sort(gsub("\\|.*", "", root_descendants)), collapse = " ")
  root_descendants_hash = sapply(root_descendants, digest, algo="md5")
  
  # read in possvm classification
  ort = read.table(sprintf("%s/%s", ort_fo, ort_fn), sep="\t", header = T, stringsAsFactors = F)
  ort$species = stringr::str_split(ort$gene, pattern = "_", simplify = T)[,1]
  ort = ort[ort$species %in% sps_list,]
  
  # subset ref to refOG of interest
  rei = ref[ref$refOG == fam & ref$gene %in% ort$gene,c("gene","refOG")]
  rei$refOG_bool = T
  # ort = ort[ort$gene %in% ref$gene,]
  
  # add ref annot to possvm classification
  ort = merge(ort[,c("gene","orthogroup")],rei, by.x = "gene",by.y = "gene", all.y = T,all.x = T)
  ort[is.na(ort$refOG),"refOG"] = "other"
  ort[is.na(ort$refOG_bool),"refOG_bool"] = F
  
  # identify equivalent possvm orthogroup
  ort_xtab = as.data.frame(table(ort$refOG,ort$orthogroup), stringsAsFactors = F)
  colnames(ort_xtab) = c("refOG","Possvm","Freq")
  ort_xtab_pos = ort_xtab[ort_xtab$refOG == fam,]
  # this line will retrieve all OGs that have a shared gene with refOG (more inclusive, implies better recall)
  ort_xtab_hits = which(ort_xtab_pos$Freq>0)
  ort_hits = ort_xtab_pos[ort_xtab_hits,"Possvm"]
  ort_hits_string = paste(ort_hits,collapse = ",")
  ort_hits_count = sum(ort_xtab_pos$Freq[ort_xtab_hits])
  
  # this line will the BEST OGs with shared genes with refOG (less inclusive, aims to identify a single hit, therefore lower recall)
  ort_xtab_max = which.max(ort_xtab_pos$Freq)
  ort_hits_best = ort_xtab_pos[ort_xtab_max,"Possvm"]
  ort_hits_best_count = sum(ort_xtab_pos$Freq[ort_xtab_max])
  
  
  # assign bool Possvm OGs
  ort$orthogroup_bool = (ort$orthogroup %in% ort_hits_best)
  ort[is.na(ort$orthogroup_bool),"orthogroup_bool"] = F
  
  # evaluate
  ix_TP = ort$orthogroup_bool & ort$refOG_bool
  ix_TN = !ort$orthogroup_bool & !ort$refOG_bool
  ix_FP = ort$orthogroup_bool & !ort$refOG_bool
  ix_FN = !ort$orthogroup_bool & ort$refOG_bool
  ev_TP = sum(ix_TP)
  ev_TN = sum(ix_TN)
  ev_FP = sum(ix_FP)
  ev_FN = sum(ix_FN)
  ge_TP = paste(ort[ix_TP,"gene"], collapse=",")
  ge_TN = paste(ort[ix_TN,"gene"], collapse=",")
  ge_FP = paste(ort[ix_FP,"gene"], collapse=",")
  ge_FN = paste(ort[ix_FN,"gene"], collapse=",")
  
  # evaluate precision = TP / (TP + FP)
  ev_precision = ev_TP / ( ev_TP + ev_FP )
  # recall = TP / (TP + FN)
  ev_recall = ev_TP / ( ev_TP + ev_FN )
  # F score = (2*Precision*Recall) / sum(Precision, Recall)
  ev_Fscore = (2*ev_precision*ev_recall) / sum(ev_recall, ev_precision)
  
  # store diagnostics
  dia = rbind(
    dia, 
    data.frame(
      refOG = fam, 
      hits=ort_hits_string, 
      ref_size=ort_hits_count,
      best_hit=ort_hits_best,
      num_best_hits=ort_hits_best_count,
      precision=ev_precision, 
      recall=ev_recall, 
      Fscore=ev_Fscore,
      iteration_id = paste(fam, ort_it, ort_frac, ort_quant),
      root = root_descendants_hash,
      method = ort_method,
      inflated_fraction = ort_frac,
      inflated_quantile = ort_quant
      
    ))
  
}




#### Merged analysis ####

dia_mid = dia [ dia$method == "mid", ]
dia_ite = dia [ dia$method == "ite", ]

# merge
dit = merge(dia_mid, dia_ite, by.x = "iteration_id", by.y = "iteration_id", all.x = FALSE, all.y = FALSE, suffixes = c(".mid", ".ite"))

# calculate difs?
dit$recall_diff = dit$recall.ite - dit$recall.mid
dit$precision_diff = dit$precision.ite - dit$precision.mid
dit$Fscore_diff = dit$Fscore.ite - dit$Fscore.mid


# variable tree: there's been a change in fscore?
ixs_variable = which(dit$Fscore_diff != 0)

# variable tree: there's been a change in root hash?
ixs_variable = which(as.character(dit$root.ite) != as.character(dit$root.mid))


# plots
pdf("results_evaluation/difference_iter_mid_inflation.pdf", height = 4, width = 4)

breaks = seq(from = -1, to = 1, length.out = 40)

# change in fscore
hist(dit$Fscore_diff[ixs_variable], breaks = breaks, col="gray90", xlim = c(-1,1), main = "Fscore in variable trees", xlab = "d(ite-mid)", ylim=c(0,100))
title(sub=sprintf(
  "f>0 = %i | f<0 = %i", sum(dit$Fscore_diff > 0 ), sum(dit$Fscore_diff < 0 )
))
abline(v=0, lty=2, col="red")

# change in precision
hist(dit$precision_diff[ixs_variable], breaks = breaks, col="gray90", xlim = c(-1,1), main = "Precision in variable trees", xlab = "d(ite-mid)",ylim=c(0,100))
title(sub=sprintf(
  "p>0 = %i | p<0 = %i", sum(dit$precision_diff > 0 ), sum(dit$precision_diff < 0 )
))
abline(v=0, lty=2, col="red")

# change in recall
hist(dit$recall_diff[ixs_variable], breaks = breaks, col="gray90", xlim = c(-1,1), main = "Recall in variable trees", xlab = "d(ite-mid)", ylim=c(0,100))
title(sub=sprintf(
  "r>0 = %i | r<0 = %i", sum(dit$recall_diff > 0 ), sum(dit$recall_diff < 0 )
))
abline(v=0, lty=2, col="red")

dev.off()

write.table(dit, "results_evaluation/difference_iter_mid_inflation.csv", quote = FALSE, sep = "\t", row.names = FALSE)
