# libraries
library(stringr)
library(scales)
library(vioplot)

# input
ort_fo = "results_trees/downsampling_alignment//"
out_fo = "results_evaluation/"

# list of datasets
set_list = list(
  list(id = "ANTP", ref = "results_searches/ANTP.seed.diamond.csv", eval_sps = c("Hsap","Mmus","Galgal","Xentro","Drer","Bflo","Dmel","Tcas","Apimel"))
)

for (set in set_list) {
  
  ref_fn = set$ref
  id = set$id
  sps_list = set$eval_sps
  
  # load data
  ref = read.table(ref_fn, sep="\t", header = F, col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"), stringsAsFactors=F)
  ref = ref[ref$pident == 100,]
  ref$refOG = stringr::str_split(ref$qseqid, "_", simplify = T)[,3]
  ref$refOG_type = stringr::str_split(ref$qseqid, "_", simplify = T)[,2]
  ref$species = stringr::str_split(ref$sseqid, pattern="_", simplify = T)[,1]
  # ref = ref[ref$species %in% sps_list,]
  # sps_list = unique(ref$species)
  ref = ref[!duplicated(ref$sseqid),]
  ref = ref[order(ref$refOG),]
  write.table(ref[,c("sseqid","refOG")], sprintf("reference_%s.csv", id), quote = F, sep="\t", row.names = F)
  write.table(ref[,c("sseqid","refOG_type")], sprintf("reference_%s.type.csv", id), quote = F, sep="\t", row.names = F)
  ref = ref[,c("sseqid","refOG")]
  
  # ignore single-gene reference orthogroups
  fam_count = table(ref$refOG)
  fam_list = names(fam_count)[fam_count>1]
  
  # store diagnostics
  
  #### per-family evaluation ####
  
  dia_sum = data.frame()
  dia_all = data.frame()
  
  perm_list = list.files("*groups.csv",path = ort_fo, full.names = T)
  
  for (per in perm_list) {
    
    dia = data.frame()
    # load orthology for this permutation
    ort_fn = sprintf(per)
    print(paste(id,per))
    
    # data
    per_fra = stringr::str_split(basename(per), pattern = "_|-", simplify = T)[,c(2)]
    per_rou = gsub("\\..*", "", stringr::str_split(basename(per), pattern = "_|-", simplify = T)[,c(3)])
    
    for (fam in fam_list) {
      
      # read in possvm classification
      ort = read.table(ort_fn, sep="\t", header = T, stringsAsFactors = F)
      ort$species = stringr::str_split(ort$gene, pattern = "_", simplify = T)[,1]
      
      # reference species that are included in this iteration
      # (at least human will always be there)
      sps_list_i = sps_list[sps_list %in% ort$species ]
      
      # subset ref to refOG of interest
      rei = ref[ref$refOG == fam & stringr::str_split(ref$sseqid, "_", simplify = T)[,1] %in% sps_list_i ,c("sseqid","refOG")]
      
      if (nrow(rei)>0) {
        rei$refOG_bool = T
        rei = rei[order(rei$refOG,rei$sseqid),]
        
        ort = ort[ort$gene %in% ref$sseqid,]
        
        # add ref annot to possvm classification
        ort = merge(ort[,c("gene","orthogroup")], rei, by.x = "gene",by.y = "sseqid", all.y = T,all.x = T)
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
        
        # this line will find the BEST OGs with shared genes with refOG (less inclusive, aims to identify a single hit, therefore lower recall)
        ort_xtab_max = which.max(ort_xtab_pos$Freq)
        ort_hits_best = ort_xtab_pos[ort_xtab_max,"Possvm"]
        ort_hits_best_count = sum(ort_xtab_pos$Freq[ort_xtab_max])
        
        # assign bool POSSVM OGs
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
            TP=ev_TP,
            TN=ev_TN,
            FP=ev_FP,
            FN=ev_FN,
            perm_fraction=per_fra,
            perm_round=per_rou
          )
        )
      }
      
    }
    
    # store round-wide diagnostics
    dia_sum = rbind(
      dia_sum,
      data.frame(
        perm_fraction=per_fra,
        perm_round=per_rou,
        precision = weighted.mean(dia$precision, dia$ref_size),
        recall = weighted.mean(dia$recall, dia$ref_size),
        Fscore = weighted.mean(dia$Fscore, dia$ref_size)
      )
    )
    dia_all = rbind(dia_all, dia)
    
  }
  
}



#### Summary plots ####

pdf(sprintf("%s/eval_%s_all_downsampling.pdf",out_fo,id),height = 5, width = 7)
layout(matrix(1:6, nrow = 2))


# precision
vioplot(dia_sum$precision ~ dia_sum$perm_fraction, xlab = "Fraction of removed species", ylab = "Precision", main="Mean precision",
        ylim=c(0.7,1), cex.axis=0.5, border="dodgerblue4", col="powderblue",rectCol="dodgerblue4", lineCol="dodgerblue4", pchMed=19, las=1)
abline(h=0.975, lty=2, col="red")
vioplot(dia_all$precision ~ dia_all$perm_fraction, xlab = "Fraction of removed species", ylab = "Precision", main="Per-gene precision",
        ylim=c(0,1), cex.axis=0.5, border="dodgerblue4", col="powderblue",rectCol="dodgerblue4", lineCol="dodgerblue4", pchMed=19, las=1)
abline(h=0.975, lty=2, col="red")

# recall
vioplot(dia_sum$recall ~ dia_sum$perm_fraction, xlab = "Fraction of removed species", ylab = "Recall", main="Mean recall", 
        ylim=c(0.7,1), cex.axis=0.5, border="dodgerblue4", col="powderblue",rectCol="dodgerblue4", lineCol="dodgerblue4", pchMed=19, las=1)
abline(h=0.936, lty=2, col="red")
vioplot(dia_all$recall ~ dia_all$perm_fraction, xlab = "Fraction of removed species", ylab = "Recall", main="Per-gene recall", 
        ylim=c(0,1), cex.axis=0.5, border="dodgerblue4", col="powderblue",rectCol="dodgerblue4", lineCol="dodgerblue4", pchMed=19, las=1)
abline(h=0.936, lty=2, col="red")

# fscore
vioplot(dia_sum$Fscore ~ dia_sum$perm_fraction, xlab = "Fraction of removed species", ylab = "F1", main="Mean F1 score", 
        ylim=c(0.7,1), cex.axis=0.5, border="dodgerblue4", col="powderblue",rectCol="dodgerblue4", lineCol="dodgerblue4", pchMed=19, las=1)
abline(h=0.941, lty=2, col="red")
vioplot(dia_all$Fscore ~ dia_all$perm_fraction, xlab = "Fraction of removed species", ylab = "F1", main="Per-gene F1 score", 
        ylim=c(0,1), cex.axis=0.5, border="dodgerblue4", col="powderblue",rectCol="dodgerblue4", lineCol="dodgerblue4", pchMed=19, las=1)
abline(h=0.941, lty=2, col="red")
dev.off()

