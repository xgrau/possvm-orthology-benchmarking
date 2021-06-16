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


#### Loop: accuracy per tree ####

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
  phy = ape::ladderize(phy)
  phy$node.label = NULL
  root_descendants_ix = phangorn::Descendants(phy, phy$edge[1,2])[[1]]
  root_descendants = phy$tip.label [ root_descendants_ix ]
  root_descendants = paste(sort(gsub("\\|.*", "", root_descendants)), collapse = " ")
  root_descendants_hash = sapply(root_descendants, digest, algo="md5")
  
  # read in possvm classification
  ort = read.table(sprintf("%s/%s", ort_fo, ort_fn), sep="\t", header = TRUE, stringsAsFactors = FALSE)
  ort$species = stringr::str_split(ort$gene, pattern = "_", simplify = TRUE)[,1]
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



#### Parallel analysis: mid & ite ####

dia_mid = dia [ dia$method == "mid", ]
dia_ite = dia [ dia$method == "ite", ]

# merge
dit = merge(dia_mid, dia_ite, by.x = "iteration_id", by.y = "iteration_id", all.x = FALSE, all.y = FALSE, suffixes = c(".mid", ".ite"))
# dit = dit[ dit$inflated_fraction.ite != "f0.100" , ]

# calculate change in recall, precision and Fscore
dit$recall_diff = dit$recall.ite - dit$recall.mid
dit$precision_diff = dit$precision.ite - dit$precision.mid
dit$Fscore_diff = dit$Fscore.ite - dit$Fscore.mid

# variable tree: there's been a change in fscore?
ixs_fscd = abs(dit$Fscore_diff) > 0
# variable tree: there's been a change in root hash?
ixs_variable = as.character(dit$root.ite) != as.character(dit$root.mid)
# ixs_variable = which(as.character(dit$root.ite) != as.character(dit$root.mid) & abs(dit$Fscore_diff) > 0)



# add factor with inflation method
dit$inflation_method = apply(stringr::str_split(dit$iteration_id, " ", simplify = TRUE)[,3:4], 1, function(x) paste(x, collapse = " "))
dit$inflation_method = factor(dit$inflation_method, levels = sort(unique(dit$inflation_method)))



pdf("results_evaluation/difference_iter_mid_inflation.pdf", height = 5, width = 5)
par(mar = c(10.1, 4.1, 4.1, 2.1))

# table: number of recall or precision improvements caused by iterative rooting
dit_sum = data.frame(
  row.names  =    levels(dit$inflation_method),
  d_recall_p =    aggregate(recall_diff ~ inflation_method, data = dit[ixs_variable,], FUN = function(x) sum(x>0))[,2],
  d_recall_n =    aggregate(recall_diff ~ inflation_method, data = dit[ixs_variable,], FUN = function(x) sum(x<0))[,2] * -1,
  d_precision_p = aggregate(precision_diff ~ inflation_method, data = dit[ixs_variable,], FUN = function(x) sum(x>0))[,2],
  d_precision_n = aggregate(precision_diff ~ inflation_method, data = dit[ixs_variable,], FUN = function(x) sum(x<0))[,2] * -1
)
# barplots
barplot(
  t(as.matrix(dit_sum))[c(1,3),], 
  beside = TRUE, las = 2, 
  ylim = c(-100, 160),
  col = c("blue","skyblue"), 
  ylab = "# changes in precision or recall")
barplot(
  t(as.matrix(dit_sum))[c(2,4),], 
  beside = TRUE, axes = FALSE, axisnames = FALSE,
  col = c("deeppink3","pink1"), 
  add = TRUE)
legend(1,170, cex = 0.7, legend = c("ite improves recall", "ite improves precision", "ite worsens recall", "ite worsens precision"), fill=c("blue","skyblue","deeppink3","pink1"), bty="n")
title(main = "Change in recall and precision in inflated trees")


# table: median of recall or precision improvement caused by iterative rooting
dit_mag = data.frame(
  row.names  =    levels(dit$inflation_method),
  d_recall_p =    aggregate(recall_diff ~ inflation_method, data = dit[ixs_variable,], FUN = function(x) median(x[x>0]) )[,2],
  d_recall_n =    aggregate(recall_diff ~ inflation_method, data = dit[ixs_variable,], FUN = function(x) median(x[x<0]) )[,2],
  d_precision_p = aggregate(precision_diff ~ inflation_method, data = dit[ixs_variable,], FUN = function(x) median(x[x>0]) )[,2],
  d_precision_n = aggregate(precision_diff ~ inflation_method, data = dit[ixs_variable,], FUN = function(x) median(x[x<0]) )[,2]
)


# barplots
barplot(
  t(as.matrix(dit_mag))[c(1,3),], 
  beside = TRUE, las = 2, 
  ylim = c(-1, 1),
  col = c("blue","skyblue"), 
  ylab = "sum(ite-mid)")
barplot(
  t(as.matrix(dit_mag))[c(2,4),], 
  beside = TRUE, axes = FALSE, axisnames = FALSE,
  col = c("deeppink3","pink1"), 
  add = TRUE)
legend("topleft", cex = 0.7, legend = c("ite improves recall", "ite improves precision", "ite worsens recall", "ite worsens precision"), fill=c("blue","skyblue","deeppink3","pink1"), bty="n")
title(main = "Change in recall and precision in inflated trees")



# final values
# distributions
vioplot::vioplot(Fscore.mid ~ inflation_method, data = dit[ixs_variable,], side="left",  col = "pink1",  colMed = "deeppink3", rectCol = "deeppink4", lineCol = "deeppink4", las = 2 , ylim = c(0,1), ylab = "F-score")
vioplot::vioplot(Fscore.ite ~ inflation_method, data = dit[ixs_variable,], side="right", col = "skyblue",colMed = "blue",      rectCol = "blue4",     lineCol = "blue4",     add= TRUE)
legend("bottomleft", cex = 0.7, legend = c("ite", "mid"), fill=c("skyblue","pink1"), bty="n")
title("F-score per inflation method")

vioplot::vioplot(precision.mid ~ inflation_method, data = dit[ixs_variable,], side="left",  col = "pink1",  colMed = "deeppink3", rectCol = "deeppink4", lineCol = "deeppink4", las = 2 , ylim = c(0,1), ylab = "Precision")
vioplot::vioplot(precision.ite ~ inflation_method, data = dit[ixs_variable,], side="right", col = "skyblue",colMed = "blue",      rectCol = "blue4",     lineCol = "blue4",     add= TRUE)
legend("bottomleft", cex = 0.7, legend = c("ite", "mid"), fill=c("skyblue","pink1"), bty="n")
title("Precision per inflation method")

vioplot::vioplot(recall.mid ~ inflation_method, data = dit[ixs_variable,], side="left",  col = "pink1",  colMed = "deeppink3", rectCol = "deeppink4", lineCol = "deeppink4", las = 2 , ylim = c(0,1), ylab = "Recall")
vioplot::vioplot(recall.ite ~ inflation_method, data = dit[ixs_variable,], side="right", col = "skyblue",colMed = "blue",      rectCol = "blue4",     lineCol = "blue4",     add= TRUE)
legend("bottomleft", cex = 0.7, legend = c("ite", "mid"), fill=c("skyblue","pink1"), bty="n")
title("Recall per inflation method")



# plots
breaks = seq(from = -1, to = 1, length.out = 40)
# change in fscore
hist(dit$Fscore_diff[ixs_variable], breaks = breaks, col="gray90", xlim = c(-1,1), main = "Fscore in variable trees", xlab = "d(ite-mid)", ylim=c(0,200))
title(sub=sprintf(
  "ite improves = %i | ite worsens = %i | n = %i", sum(dit$Fscore_diff[ixs_variable] > 0), sum(dit$Fscore_diff[ixs_variable] < 0), length(ixs_variable)
))
abline(v=0, lty=2, col="red")

# change in precision
hist(dit$precision_diff[ixs_variable], breaks = breaks, col="gray90", xlim = c(-1,1), main = "Precision in variable trees", xlab = "d(ite-mid)",ylim=c(0,200))
title(sub=sprintf(
  "ite improves = %i | ite worsens = %i | n = %i", sum(dit$precision_diff[ixs_variable] > 0), sum(dit$precision_diff[ixs_variable] < 0), length(ixs_variable)
))
abline(v=0, lty=2, col="red")

# change in recall
hist(dit$recall_diff[ixs_variable], breaks = breaks, col="gray90", xlim = c(-1,1), main = "Recall in variable trees", xlab = "d(ite-mid)", ylim=c(0,200))
title(sub=sprintf(
  "ite improves = %i | ite worsens = %i | n = %i", sum(dit$recall_diff[ixs_variable] > 0), sum(dit$recall_diff[ixs_variable] < 0), length(ixs_variable)
))
abline(v=0, lty=2, col="red")

dev.off()



# save?
write.table(dit, "results_evaluation/difference_iter_mid_inflation.csv", quote = FALSE, sep = "\t", row.names = FALSE)


