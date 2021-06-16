# libraries
library(stringr)
library(alluvial)
library(scales)
library(vioplot)
library(fossil)

# input
out_fo = "results_evaluation/"

# list of datasets
set_list = list(
  list(id = "antp_all_ite_mcl",  ref = "results_searches/ANTP.seed.diamond-edit.csv", ort_fn = "results_trees/ANTP.possom.ortholog_groups.csv",      sps_list = c("Drer","Galgal","Hsap","Mmus","Xentro","Dmel","Apimel","Tcas","Bflo")),
  list(id = "antp_ver_ite_mcl",  ref = "results_searches/ANTP.seed.diamond-edit.csv", ort_fn = "results_trees/ANTP.possom.ortholog_groups.csv",      sps_list = c("Drer","Galgal","Hsap","Mmus","Xentro")),
  list(id = "antp_ins_ite_mcl",  ref = "results_searches/ANTP.seed.diamond-edit.csv", ort_fn = "results_trees/ANTP.possom.ortholog_groups.csv",      sps_list = c("Dmel","Apimel","Tcas")),
  list(id = "antp_all_mid_mcl",  ref = "results_searches/ANTP.seed.diamond-edit.csv", ort_fn = "results_trees/ANTP.possom_mid.ortholog_groups.csv",  sps_list = c("Drer","Galgal","Hsap","Mmus","Xentro","Dmel","Apimel","Tcas","Bflo")),
  list(id = "antp_all_ite_mclw", ref = "results_searches/ANTP.seed.diamond-edit.csv", ort_fn = "results_trees/ANTP.possom_mclw.ortholog_groups.csv", sps_list = c("Drer","Galgal","Hsap","Mmus","Xentro","Dmel","Apimel","Tcas","Bflo")),
  list(id = "antp_all_ite_lou",  ref = "results_searches/ANTP.seed.diamond-edit.csv", ort_fn = "results_trees/ANTP.possom_lou.ortholog_groups.csv",  sps_list = c("Drer","Galgal","Hsap","Mmus","Xentro","Dmel","Apimel","Tcas","Bflo")),
  list(id = "antp_all_ite_lpa",  ref = "results_searches/ANTP.seed.diamond-edit.csv", ort_fn = "results_trees/ANTP.possom_lpa.ortholog_groups.csv",  sps_list = c("Drer","Galgal","Hsap","Mmus","Xentro","Dmel","Apimel","Tcas","Bflo")),
  list(id = "antp_all_phylo_h",  ref = "results_searches/ANTP.seed.diamond-edit.csv", ort_fn = "results_trees/ANTP.ccsh.ortholog_groups_from_ccs.csv",sps_list = c("Drer","Galgal","Hsap","Mmus","Xentro","Dmel","Apimel","Tcas","Bflo")),
  list(id = "antp_all_phylo_d",  ref = "results_searches/ANTP.seed.diamond-edit.csv", ort_fn = "results_trees/ANTP.ccsd.ortholog_groups_from_ccs.csv",sps_list = c("Drer","Galgal","Hsap","Mmus","Xentro","Dmel","Apimel","Tcas","Bflo")),
  list(id = "antp_all_brc_m50",  ref = "results_searches/ANTP.seed.diamond-edit.csv", ort_fn = "results_branchclust/ANTP.bc_clusters_50.csv",        sps_list = c("Drer","Galgal","Hsap","Mmus","Xentro","Dmel","Apimel","Tcas","Bflo")),
  list(id = "antp_all_brc_m60",  ref = "results_searches/ANTP.seed.diamond-edit.csv", ort_fn = "results_branchclust/ANTP.bc_clusters_60.csv",        sps_list = c("Drer","Galgal","Hsap","Mmus","Xentro","Dmel","Apimel","Tcas","Bflo")),
  list(id = "antp_all_brc_m70",  ref = "results_searches/ANTP.seed.diamond-edit.csv", ort_fn = "results_branchclust/ANTP.bc_clusters_70.csv",        sps_list = c("Drer","Galgal","Hsap","Mmus","Xentro","Dmel","Apimel","Tcas","Bflo")),
  list(id = "antp_all_brc_m80",  ref = "results_searches/ANTP.seed.diamond-edit.csv", ort_fn = "results_branchclust/ANTP.bc_clusters_80.csv",        sps_list = c("Drer","Galgal","Hsap","Mmus","Xentro","Dmel","Apimel","Tcas","Bflo")),
  list(id = "prds_all_ite_mcl",  ref = "results_searches/PRD.seed.diamond.csv",       ort_fn = "results_trees/PRD.possom.ortholog_groups.csv",       sps_list = c("Drer","Galgal","Hsap","Mmus","Xentro","Dmel","Apimel","Tcas","Bflo")),
  list(id = "tale_all_ite_mcl",  ref = "results_searches/TALE.seed.diamond.csv",      ort_fn = "results_trees/TALE.possom.ortholog_groups.csv",      sps_list = c("Drer","Galgal","Hsap","Mmus","Xentro","Dmel","Apimel","Tcas","Bflo"))
)


# set_list = set_list[1]
eval_list_ids = unlist(lapply(set_list, function(x) x$id))



#### Global Rand index ####

rand_index_adj = c()
rand_index_raw = c()
message("Rand Index...")
for (set in set_list) {
  
  # input info
  ref_fn = set$ref
  id = set$id
  ort_fn = set$ort_fn
  sps_list =set$sps_list
  prefix =set$prefix
  
  # load data
  ref = read.table(ref_fn, sep="\t", header = F, col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"), stringsAsFactors=F)
  ref = ref[ref$pident == 100,]
  ref$refOG = stringr::str_split(ref$qseqid, "_", simplify = T)[,3]
  ref$refOG_type = stringr::str_split(ref$qseqid, "_", simplify = T)[,2]
  ref$species = stringr::str_split(ref$sseqid, pattern="_", simplify = T)[,1]
  ref = ref[ref$species %in% sps_list,]
  ref = ref[!duplicated(ref$sseqid),]
  ref = ref[order(ref$refOG),]
  ref = ref[,c("sseqid","refOG")]
  
  # ignore single-gene reference orthogroups
  fam_count = table(ref$refOG)
  fam_list = names(fam_count)[fam_count>1]
  
  # read in possvm classification
  ort = read.table(ort_fn, sep="\t", header = TRUE, stringsAsFactors = FALSE)
  colnames(ort)[1:2] = c("gene","orthogroup")
  ort$species = stringr::str_split(ort$gene, pattern = "_", simplify = T)[,1]
  
  # merge
  ort = merge(ort[,c("gene","orthogroup")], ref, by.x = "gene", by.y = "sseqid", all.x = FALSE, all.y = TRUE)

  # calculate indexes (raw and adjusted)
  clu_r = as.numeric(as.factor(ort$refOG))
  clu_q = as.numeric(as.factor(ort$orthogroup))
  i_adj = fossil::adj.rand.index(clu_r, clu_q)
  i_raw = fossil::rand.index(clu_r, clu_q)

  # add to vector
  rand_index_adj = c(rand_index_adj, i_adj)
  rand_index_raw = c(rand_index_raw, i_raw)
  
}

rand_index_df = data.frame(
  method = eval_list_ids, 
  rand_index_adj = rand_index_adj,
  rand_index_raw = rand_index_raw
)

write.table(rand_index_df, sprintf("%s/all_methods.rand_index.csv", out_fo), quote = FALSE, sep = "\t", row.names = FALSE)

# plot rand index
pdf(sprintf("%s/all_methods.rand_index.pdf", out_fo), width = 6, height = 5)
par(mar = c(10.1, 4.1, 4.1, 2.1))
mat = t(as.matrix(rand_index_df[,2:3]))
colnames(mat) = eval_list_ids

# side by side
barplot(mat, beside = TRUE, las = 2, col = c("blue","gray90"), ylim = c(0,1), xlim = c(0,60), ylab = "Rand index")
title(main = "Rand index per dataset")
abline(h=mat[1,1], col = "purple", lty=2)
legend("topright", legend = c("adjusted", "raw"), fill = c("blue","gray90"))

# only adjusted
barplot(rand_index_df$rand_index_adj, names.arg = rand_index_df$method, las = 2, col = "blue", ylim = c(0,1), ylab = "Adjusted Rand index")
title(main = "Rand index per dataset")
abline(h=mat[1,1], col = "purple", lty=2)
dev.off()




#### OG-wise accuracy ####

# loop through datasets and calculate accuracy
dit = data.frame()
for (set in set_list) {
  
  # input info
  ref_fn = set$ref
  id = set$id
  ort_fn = set$ort_fn
  sps_list =set$sps_list
  prefix =set$prefix
  
  # load data
  ref = read.table(ref_fn, sep="\t", header = F, col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"), stringsAsFactors=F)
  ref = ref[ref$pident == 100,]
  ref$refOG = stringr::str_split(ref$qseqid, "_", simplify = T)[,3]
  ref$refOG_type = stringr::str_split(ref$qseqid, "_", simplify = T)[,2]
  ref$species = stringr::str_split(ref$sseqid, pattern="_", simplify = T)[,1]
  ref = ref[ref$species %in% sps_list,]
  ref = ref[!duplicated(ref$sseqid),]
  ref = ref[order(ref$refOG),]
  write.table(ref[,c("sseqid","refOG")], sprintf("%s/reference.%s.csv", out_fo, id), quote = F, sep="\t", row.names = F)
  # write.table(ref[,c("sseqid","refOG_type")], sprintf("reference_%s_%s.type.csv", id, spsi), quote = F, sep="\t", row.names = F)
  ref = ref[,c("sseqid","refOG")]
  
  # ignore single-gene reference orthogroups
  fam_count = table(ref$refOG)
  fam_list = names(fam_count)[fam_count>1]
  
  # store diagnostics
  
  #### per-family evaluation ####
  
  dia = data.frame()
  
  pdf(sprintf("%s/eval.%s.alluvial.pdf",out_fo, id),height = 4, width = 6)
  for (fam in fam_list) {
    
    message(paste(id,fam))
    
    # read in possvm classification
    ort = read.table(ort_fn, sep="\t", header = TRUE, stringsAsFactors = FALSE)
    colnames(ort)[1:2] = c("gene","orthogroup")
    ort$species = stringr::str_split(ort$gene, pattern = "_", simplify = T)[,1]
    
    # subset ref to refOG of interest
    rei = ref[ref$refOG == fam,c("sseqid","refOG")]
    rei$refOG_bool = T
    rei = rei[order(rei$refOG,rei$sseqid),]
    ort = ort[ort$gene %in% ref$sseqid,]
    
    # add ref annot to possvm classification
    ort = merge(ort[,c("gene","orthogroup")], rei, by.x = "gene",by.y = "sseqid", all.y = T, all.x = T)
    ort[is.na(ort$refOG),"refOG"] = "other"
    ort[is.na(ort$refOG_bool),"refOG_bool"] = F
    
    # identify equivalent possvm orthogroup
    ort_xtab = as.data.frame(table(ort$refOG,ort$orthogroup), stringsAsFactors = F)
    colnames(ort_xtab) = c("refOG","Possvm","Freq")
    ort_xtab_pos = ort_xtab[ort_xtab$refOG == fam,]
    ort_xtab_oth = ort_xtab[ort_xtab$refOG == "other",]
    ort_xtab_frac = ort_xtab_pos$Freq / (ort_xtab_pos$Freq+ort_xtab_oth$Freq)
    
    
    # this block will find the BEST OGs with shared genes with refOG 
    # (less inclusive, aims to identify a single hit, therefore lower recall)
    ort_xtab_max = which.max(ort_xtab_pos$Freq)
    ort_hits_best = ort_xtab_pos[ort_xtab_max,"Possvm"]
    ort_hits_best_count = sum(ort_xtab_pos$Freq[ort_xtab_max])
    # assign bool POSSVM OGs
    ort$orthogroup_bool = (ort$orthogroup %in% ort_hits_best)
    ort[is.na(ort$orthogroup_bool),"orthogroup_bool"] = FALSE
    
    # evaluate: based on number of present genes
    ix_TP = ort$orthogroup_bool & ort$refOG_bool
    ix_FP = ort$orthogroup_bool & !ort$refOG_bool
    ix_FN = !ort$orthogroup_bool & ort$refOG_bool
    ev_TP = sum(ix_TP)
    ev_FP = sum(ix_FP)
    ev_FN = sum(ix_FN)
    ge_TP = paste(ort[ix_TP,"gene"], collapse=",")
    ge_FP = paste(ort[ix_FP,"gene"], collapse=",")
    ge_FN = paste(ort[ix_FN,"gene"], collapse=",")
    
    
    # this block will add additional OGs that also align with the refOG (more inclusive)
    # defined as: possvm OG which has at least one gene in the refOG, and at least 50% of its contents are part of the refOG
    ort_xtab_hits = which(ort_xtab_pos$Freq>0 & ort_xtab_frac >= 0.5)
    ort_hits = ort_xtab_pos[c(ort_xtab_max, ort_xtab_hits),"Possvm"]
    ort_hits_string = paste(ort_hits,collapse = ",")
    ort_hits_count = sum(ort_xtab_pos$Freq[ort_xtab_hits])
    # assign bool POSSVM OGs
    ort$orthogroup_bool_inclusive = (ort$orthogroup %in% ort_hits)
    ort[is.na(ort$orthogroup_bool_inclusive),"orthogroup_bool_inclusive"] = FALSE
    
    # evaluate: based on number of present genes
    ix_TPi = ort$orthogroup_bool_inclusive & ort$refOG_bool
    ix_FPi = ort$orthogroup_bool_inclusive & !ort$refOG_bool
    ix_FNi = !ort$orthogroup_bool_inclusive & ort$refOG_bool
    ev_TPi = sum(ix_TPi)
    ev_FPi = sum(ix_FPi)
    ev_FNi = sum(ix_FNi)

        
    # calculate Rand index
    rand_index = fossil::adj.rand.index(as.numeric(ort$refOG_bool), as.numeric(ort$orthogroup_bool))
    
    # evaluate: based on number of present gene pairs
    # get ordered pairs from reference
    # ref_genes = ort$gene [ ort$refOG_bool ]
    # ref_pairs = t(combn(ref_genes, 2))
    # ref_pairs = unique(apply(ref_pairs, 1, function(x) paste(sort(x), collapse = ",")))
    # # get ordered pairs from query
    # ogi_genes = ort$gene [ ort$orthogroup_bool ]
    # ogi_pairs = t(combn(ogi_genes, 2))
    # ogi_pairs = unique(apply(ogi_pairs, 1, function(x) paste(sort(x), collapse = ",")))
    # # count
    # ev_TP = sum(ogi_pairs %in% ref_pairs)
    # ev_FP = sum(!(ogi_pairs %in% ref_pairs))
    # ev_FN = sum(!(ref_pairs %in% ogi_pairs))
    # ge_TP = paste(unique(strsplit(paste(ogi_pairs [ ogi_pairs %in% ref_pairs ],    collapse = ","), ",")[[1]]), collapse = ",")
    # ge_FP = paste(unique(strsplit(paste(ogi_pairs [ !(ogi_pairs %in% ref_pairs) ], collapse = ","), ",")[[1]]), collapse = ",")
    # ge_FN = paste(unique(strsplit(paste(ref_pairs [ !(ref_pairs %in% ogi_pairs) ], collapse = ","), ",")[[1]]), collapse = ",")
    
    
    # evaluate precision
    ev_precision   = ev_TP  / ( ev_TP  + ev_FP )
    ev_precision_i = ev_TPi / ( ev_TPi + ev_FPi )
    # recall
    ev_recall   = ev_TPi / ( ev_TPi + ev_FN )
    ev_recall_i = ev_TPi / ( ev_TPi + ev_FNi )
    # F score
    ev_Fscore   = (2 * ev_precision * ev_recall)     / sum(ev_recall, ev_precision)
    ev_Fscore_i = (2 * ev_precision_i * ev_recall_i) / sum(ev_recall_i, ev_precision_i)
    
    # alluvial plot
    ort_xtab = ort_xtab[ort_xtab$Possvm %in% ort_hits, ]
    if (nrow(ort_xtab)>1){
      alluvial(
        ort_xtab[,1:2], 
        freq = ort_xtab[,3], 
        col = c("lightblue"), gap.width = 0.1, blocks = T, 
        border=NA, cex=0.8)
    } else {
      plot(0,0, xlab = "", ylab = "", frame.plot = F, axes = F, col="blue", pch=19)
    }
    title(
      main=sprintf("%s\nn=%i",fam, length(rei$sseqid)),
      sub=sprintf(
        "%s is %s (%s)\nPrecision = %.2f | Recall = %.2f | F-score = %.2f", 
        fam, ort_hits_best, ort_hits_string, ev_precision, ev_recall, ev_Fscore),
      cex.main=0.8,
      cex.sub=0.7)
    
    # store diagnostics
    dia = rbind(
      dia, 
      data.frame(
        refOG    = fam, 
        best_hit = ort_hits_best,
        num_best_hits = ort_hits_best_count,
        all_hits = ort_hits_string, 
        ref_size = ort_hits_count,
        precision = ev_precision,  # precision and recall (best OG)
        recall    = ev_recall, 
        Fscore    = ev_Fscore,
        TP        = ev_TP,
        FP        = ev_FP,
        FN        = ev_FN,
        precision_i = ev_precision_i, # same, but inclusive
        recall_i    = ev_recall_i, 
        Fscore_i    = ev_Fscore_i,
        TPi         = ev_TPi,
        FPi         = ev_FPi,
        FNi         = ev_FNi,
        rand_index  = rand_index,
        TP_genes = ge_TP, # this only available for best OG
        FP_genes = ge_FP,
        FN_genes = ge_FN,
        method = id
      ))
    
  }
  dev.off()
  
  
  
  
  #### Summary plots ####
  
  # dists of fscore, precision and recall
  pdf(sprintf("%s/eval.%s.summary.pdf",out_fo, id),height = 5, width = 7)
  layout(matrix(1:6, nrow = 2))
  
  # fscore
  hist(dia$Fscore, breaks = 10, xlim = c(0,1),main="F-score", col="blue", ylim=c(0,50),border = "white", xlab = "F-score", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("av = %.3f (inc = %.3f)", weighted.mean(dia$Fscore, dia$ref_size), weighted.mean(dia$Fscore_i, dia$ref_size)))
  abline(v=0.95, lty=2, col="grey")
  # precision ecdf
  plot(sort(dia$Fscore), col="blue", ylab = "Fscore", ylim = c(0,1), main="Fscore", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("n(>=.95) = %i", sum(dia$Fscore >= 0.95)))
  abline(h=0.95, lty=2, col="grey")
  
  # precision
  hist(dia$precision, breaks = 10, xlim = c(0,1),main="Precision", col="blue",ylim=c(0,50), border = "white", xlab = "Precision", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("av = %.3f (inc = %.3f)", weighted.mean(dia$precision, dia$ref_size), weighted.mean(dia$precision_i, dia$ref_size)))
  abline(v=0.95, lty=2, col="grey")
  # precision ecdf
  plot(sort(dia$precision), col="blue", ylab = "Precision", ylim = c(0,1), main="Precision", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("n(>=.95) = %i", sum(dia$precision >= 0.95)))
  abline(h=0.95, lty=2, col="grey")

  # recall
  hist(dia$recall, breaks = 10, xlim = c(0,1),main="Recall", col="blue",ylim=c(0,50), border = "white", xlab = "Recall", cex.axis=0.9, cex.lab=0.9)
  ek_mean = ev_TP / ( ev_TP + ev_FP )
  title(sub=sprintf("av = %.3f (inc = %.3f)", weighted.mean(dia$recall, dia$ref_size), weighted.mean(dia$recall_i, dia$ref_size)))
  abline(v=0.95, lty=2, col="grey")
  # recall ecdf
  plot(sort(dia$recall), col="blue", ylab = "Recall", ylim = c(0,1), main="Recall", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("n(>=0.95) = %i", sum(dia$recall >= 0.95)))
  abline(h=0.95, lty=2, col="grey")
  
  # rand_index
  hist(dia$rand_index, breaks = 10, xlim = c(0,1),main="Rand index", col="blue",ylim=c(0,50), border = "white", xlab = "Rand index", cex.axis=0.9, cex.lab=0.9)
  ek_mean = ev_TP / ( ev_TP + ev_FP )
  title(sub=sprintf("av = %.3f (inc = %.3f)", weighted.mean(dia$rand_index, dia$ref_size), weighted.mean(dia$rand_index, dia$ref_size)))
  abline(v=0.95, lty=2, col="grey")
  # rand_index ecdf
  plot(sort(dia$rand_index), col="blue", ylab = "Rand index", ylim = c(0,1), main="Rand index", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("n(>=0.95) = %i", sum(dia$rand_index >= 0.95)))
  abline(h=0.95, lty=2, col="grey")
  
  
  # compare precision and recall
  plot(dia$precision, dia$recall, xlim = c(0,1), ylim=c(0,1), xlab = "Precision", ylab="Recall",
       col=alpha("blue", 0.6), main="Precision & recall", cex.axis=0.9, cex.lab=0.9)
  text(dia$precision, dia$recall, labels = dia$refOG, col=alpha("lightblue",0.8),cex=0.8)
  
  
  # number of ogs per reference group
  plot(dia$ref_size, dia$num_best_hits, col="blue", xlab = "# genes in ref OG", ylab = "# genes in best OG", cex.axis=0.9, cex.lab=0.9)
  title(main = "size OGs",
    sub=sprintf(
    "n(best >=0.9 all) = %i (%.3f | r=%.3f)\nn(best = all) = %i (p=%.3f | r=%.3f)\nn(best < all) = %i (p=%.3f | r=%.3f)", 
    sum(dia$num_best_hits / dia$ref_size >= 0.9), 
    mean(dia[dia$num_best_hits / dia$ref_size >= 0.9,"precision"]), 
    mean(dia[dia$num_best_hits / dia$ref_size >= 0.9,"recall"]),
    sum(dia$num_best_hits == dia$ref_size), 
    mean(dia[dia$num_best_hits == dia$ref_size,"precision"]), 
    mean(dia[dia$num_best_hits == dia$ref_size,"recall"]),
    sum(dia$num_best_hits < dia$ref_size), 
    mean(dia[dia$num_best_hits < dia$ref_size,"precision"]), 
    mean(dia[dia$num_best_hits < dia$ref_size,"recall"])
  ),
  cex.sub=0.6)
  abline(a=0, b=1, lty=2, col="grey")
  
  
  # size of main groups
  # barplot(dia$ref_size, names.arg = dia$refOG,las=2, cex.names =  0.5, border = "white")
  # barplot(dia$num_best_hits, cex.names =  0.5, add = T, col = "blue", axes = F, border = "white")
  layout(matrix(1:2, nrow = 2))
  dia_count = t(as.matrix(dia[,c("num_best_hits")]))
  dia_count = rbind(dia_count, dia$ref_size - dia$num_best_hits)
  barplot(dia_count, names.arg = dia$refOG, las=2, cex.names =  0.5, col = c("blue","gray"), border = "white",
          main="number of reference genes in main and other OGs")
  legend("topright", fill=c("blue","gray"), border = "white",legend = c("main","other"), bty = "n", cex=0.5)
  
  
  # number of OGs mapping to each refOG
  dia_num_OGs = lapply(as.character(dia$all_hits), function(x) length(strsplit(x, ",")[[1]]))
  names(dia_num_OGs) = dia$refOG
  barplot(as.numeric(dia_num_OGs), border = "white", main="num OGs per refOG",
          names.arg = dia$refOG, las=2, cex.names =  0.5)
  
  dev.off()
  
  # save table
  write.table(dia, file=sprintf("%s/eval.%s.summary.csv",out_fo,id), quote = F, sep="\t",row.names = F)
  dit = rbind(dit, dia)
  
}


#### Compare all datasets ####

eval_list_ids = unlist(lapply(set_list, function(x) x$id))
dit$method = factor(dit$method, levels = eval_list_ids)

pdf(sprintf("%s/all_methods.pdf", out_fo), width = 6, height = 5)
par(mar = c(10.1, 4.1, 4.1, 2.1))

# fscore
vioplot::vioplot(Fscore ~ method, data = dit, las = 2, col = "skyblue", main = "F-score", ylim = c(0,1))
vioplot::vioplot(Fscore ~ method, data = dit, las = 2, col = "skyblue", main = "F-score", ylim = c(0,1), side = "left")
vioplot::vioplot(Fscore_i ~ method, data = dit, las = 2, col = "orange", main = "F-score", ylim = c(0,1), side = "right", add = TRUE)
legend("bottomleft", fill = c("skyblue","orange"), legend = c("best OG", "inclusive OG"), bty = "n")

# precision
vioplot::vioplot(precision ~ method, data = dit, las = 2, col = "skyblue", main = "Precision", ylim = c(0,1))
vioplot::vioplot(precision ~ method, data = dit, las = 2, col = "skyblue", main = "Precision", ylim = c(0,1), side = "left")
vioplot::vioplot(precision_i ~ method, data = dit, las = 2, col = "orange", main = "Precision", ylim = c(0,1), side = "right", add = TRUE)
legend("bottomleft", fill = c("skyblue","orange"), legend = c("best OG", "inclusive OG"), bty = "n")

# recall
vioplot::vioplot(recall ~ method, data = dit, las = 2, col = "skyblue", main = "Recall", ylim = c(0,1))
vioplot::vioplot(recall ~ method, data = dit, las = 2, col = "skyblue", main = "Recall", ylim = c(0,1), side = "left")
vioplot::vioplot(recall_i ~ method, data = dit, las = 2, col = "orange", main = "Recall", ylim = c(0,1), side = "right", add = TRUE)
legend("bottomleft", fill = c("skyblue","orange"), legend = c("best OG", "inclusive OG"), bty = "n")

# rand_index
vioplot::vioplot(rand_index ~ method, data = dit, las = 2, col = "skyblue", main = "Rand index", ylim = c(0,1))

dev.off()
