# libraries
library(stringr)
library(alluvial)
library(scales)

# list of datasets
eval_list = list(
  list(id = "Possvm (best, MCL)",     ref = "results_evaluation/eval_o2o_ANTP_all_summary.csv"),
  list(id = "Possvm (all, MCL)",      ref = "results_evaluation/eval_o2m_ANTP_all_summary.csv"),
  list(id = "Possvm (best, Louvain)", ref = "results_evaluation/eval_o2o_lou_ANTP_all_summary.csv"),
  list(id = "Possvm (best, LPA)",     ref = "results_evaluation/eval_o2o_lpa_ANTP_all_summary.csv"),
  list(id = "Possvm (best, CNM)",     ref = "results_evaluation/eval_o2o_gre_ANTP_all_summary.csv"),
  list(id = "PhylomeDB",              ref = "results_evaluation/eval_ccsfocus_o2o_ANTP_all_summary.csv"),
  list(id = "BranchClust M=50%",      ref = "results_evaluation/eval_branchclust-50_o2o_ANTP_all_summary.csv"),
  list(id = "BranchClust M=60%",      ref = "results_evaluation/eval_branchclust-60_o2o_ANTP_all_summary.csv"),
  list(id = "BranchClust M=70%",      ref = "results_evaluation/eval_branchclust-70_o2o_ANTP_all_summary.csv"),
  list(id = "BranchClust M=80%",      ref = "results_evaluation/eval_branchclust-80_o2o_ANTP_all_summary.csv")
)

tat = data.frame()
for (eval in eval_list) {
  
  tai = read.table(eval$ref, sep="\t", header = TRUE, stringsAsFactors = FALSE)
  tai$method = eval$id
  tat = rbind(tat, tai)
  
}

eval_list_ids = unlist(lapply(eval_list, function(x) x$id))
tat$method = factor(tat$method, levels = eval_list_ids)

pdf("all_methods.pdf", width = 6, height = 5)
par(mar = c(10.1, 4.1, 4.1, 2.1))
boxplot(Fscore ~ method, data = tat, las = 2, col = "skyblue", main = "F-score", ylim = c(0,1))
boxplot(precision ~ method, data = tat, las = 2, col = "skyblue", main = "Precision", ylim = c(0,1))
boxplot(recall ~ method, data = tat, las = 2, col = "skyblue", main = "Recall", ylim = c(0,1))
dev.off()
