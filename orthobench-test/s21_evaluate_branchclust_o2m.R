# libraries
library(stringr)
library(alluvial)
library(scales)

# input
ref_fn = "refOGs.csv"

ort_sets = list(
  list(id="raw", ort_fo="orthobench_trees/raw/", many=50),
  list(id="raw", ort_fo="orthobench_trees/raw/", many=60),
  list(id="raw", ort_fo="orthobench_trees/raw/", many=70),
  list(id="raw", ort_fo="orthobench_trees/raw/", many=80)
)


# load data
ref = read.table(ref_fn, sep="\t", header = F, col.names = c("refOG", "gene"), stringsAsFactors=F)
ref$species = stringr::str_split(ref$gene, pattern="_", simplify = T)[,1]
fam_list = unique(ref$refOG)
sps_list = unique(ref$species)

#### per-family evaluation ####

for (ort_set in ort_sets) {
  
  
  ort_fo = ort_set$ort_fo
  id = ort_set$id
  many = ort_set$many
  dia = data.frame()
  
  
  pdf(sprintf("results_evaluation/eval_branchclust-%s_o2m_%s_classification_alluvial.pdf", many, id),height = 8, width = 7)
  for (fam in fam_list) {
    
    ort_fn = sprintf("results_branchclust/%s.bc_clusters_%s.csv", fam, many)
    print(fam)
    
    if( file.exists(ort_fn)) {
      
      # read in possvm classification
      ort = read.table(ort_fn, sep="\t", header = F, stringsAsFactors = F, col.names = c("gene","orthogroup"))
      
      if (nrow(ort)>0) {
      ort$species = stringr::str_split(ort$gene, pattern = "_", simplify = T)[,1]
      ort = ort[ort$species %in% sps_list,]
      
      # subset ref to refOG of interest
      rei = ref[ref$refOG == fam ,c("gene","refOG")]
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
      ort$orthogroup_bool = (ort$orthogroup %in% ort_hits)
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
      
      # alluvial plot
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
        main=sprintf("%s\nn=%i",fam, length(rei$gene)),
        sub=sprintf("%s is %s (%s)\nPrecision = %.2f | Recall = %.2f | F-score = %.2f", fam,ort_hits_best, ort_hits_string,ev_precision, ev_recall, ev_Fscore),
        cex.main=0.8,
        cex.sub=0.7)
      
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
          TP=ev_TP,
          TN=ev_TN,
          FP=ev_FP,
          FN=ev_FN,
          TP_genes=ge_TP,
          FP_genes=ge_FP,
          FN_genes=ge_FN
        ))
      
      }
    }
  }
  dev.off()
  
  
  #### Summary plots ####
  
  # compare precsion and recall
  
  pdf(sprintf("results_evaluation/eval_branchclust-%s_o2m_%s_summary.pdf", many, id),height = 5, width = 7)
  layout(matrix(1:6, nrow = 2))
  plot(dia$precision, dia$recall, xlim = c(0,1), ylim=c(0,1), xlab = "Precision", ylab="Recall",
       col=alpha("blue", 0.6), main="Precision & recall", cex.axis=0.9, cex.lab=0.9)
  text(dia$precision, dia$recall, labels = dia$refOG, col=alpha("lightblue",0.8),cex=0.8)
  
  # dists of fscore, precision and recall
  hist(dia$Fscore, breaks = 10, xlim = c(0,1),main="F-score", col="blue", ylim=c(0,50),border = "white", xlab = "F-score", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("n(>=0.95) = %i | mean = %.3f", sum(dia$Fscore >= 0.95), weighted.mean(dia$Fscore, dia$ref_size)))
  abline(v=0.95, lty=2, col="grey")
  
  hist(dia$precision, breaks = 10, xlim = c(0,1),main="Precision", col="blue",ylim=c(0,50), border = "white", xlab = "Precision", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("n(>=0.95) = %i | mean = %.3f", sum(dia$precision >= 0.95), weighted.mean(dia$precision, dia$ref_size)))
  abline(v=0.95, lty=2, col="grey")
  
  hist(dia$recall, breaks = 10, xlim = c(0,1),main="Recall", col="blue",ylim=c(0,50), border = "white", xlab = "Recall", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("n(>=0.95) = %i | mean = %.3f", sum(dia$recall >= 0.95), weighted.mean(dia$recall, dia$ref_size)))
  abline(v=0.95, lty=2, col="grey")
  
  plot(sort(dia$precision), col="blue", ylab = "Precision", ylim = c(0,1), main="Precision", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("n(>=0.95) = %i", sum(dia$precision >= 0.95)))
  abline(h=0.95, lty=2, col="grey")
  
  plot(sort(dia$recall), col="blue", ylab = "Recall", ylim = c(0,1), main="Recall", cex.axis=0.9, cex.lab=0.9)
  title(sub=sprintf("n(>=0.95) = %i", sum(dia$recall >= 0.95)))
  abline(h=0.95, lty=2, col="grey")
  
  # number of ogs per reference group
  plot(dia$ref_size, dia$num_best_hits, col="blue", xlab = "# genes in ref OG", ylab = "# genes in best OG", cex.axis=0.9, cex.lab=0.9,
       xlim=c(0,100), ylim = c(0,100))
  title(sub=sprintf("n(best >=0.9 all) = %i (%.3f | r=%.3f)\nn(best = all) = %i (p=%.3f | r=%.3f)\nn(best < all) = %i (p=%.3f | r=%.3f)", 
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
  dia_num_OGs = lapply(as.character(dia$hits), function(x) length(strsplit(x, ",")[[1]]))
  names(dia_num_OGs) = dia$refOG
  barplot(as.numeric(dia_num_OGs), border = "white", main="num OGs per refOG",
          names.arg = dia$refOG, las=2, cex.names =  0.5)
  
  dev.off()
  
  # save table
  write.table(dia, file=sprintf("results_evaluation/eval_branchclust-%s_o2m_%s_summary.csv", many, id), quote = F, sep="\t",row.names = F)
  
}