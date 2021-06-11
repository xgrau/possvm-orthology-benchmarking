# libraries
library(scales)

# input
mid = read.table("results_evaluation/eval_o2o_raw_summary.csv",      header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names = 1)
ite = read.table("results_evaluation/eval_o2o_iter_raw_summary.csv", header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names = 1)


# ogs_with_changed_root = c("RefOG048","RefOG063","RefOG054","RefOG047","RefOG064","RefOG007","RefOG009","RefOG014","RefOG016","RefOG019","RefOG026","RefOG029","RefOG033","RefOG035","RefOG052","RefOG055","RefOG057")

# subset
# mid = mid [ rownames(mid) %in% ogs_with_changed_root , c("precision","recall","Fscore") ]
# ite = ite [ rownames(ite) %in% ogs_with_changed_root , c("precision","recall","Fscore") ]
mid = mid [ rownames(mid) %in% rownames(ite) , c("precision","recall","Fscore") ]
ite = ite [ rownames(ite) %in% rownames(mid) , c("precision","recall","Fscore") ]

dif = mid 
dif$precision = ite$precision - mid$precision 
dif$recall = ite$recall - mid$recall 
dif$Fscore = ite$Fscore - mid$Fscore
dif = dif[ order(dif$recall), ]


pdf("difference_iter_mid.pdf", height = 12, width = 3.5)
barplot(as.matrix(t(dif)), beside = TRUE, col = c("purple","orange","blue"), border = c("purple4","orange4","blue4"), xlab = "d(iterative-midpoint)", horiz = TRUE, cex.names = 0.5, las=1, xlim = c(-1,1))
legend("bottomright", legend = c("Precision", "Recall","F-score"), fill = c("purple","orange","blue"), border = c("purple4","orange4","blue4"), bty="n")
dev.off()
