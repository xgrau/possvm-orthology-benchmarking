# libraries
library(scales)

# input
mid = read.table("results_evaluation/eval_o2o_raw_summary.csv",      header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names = 1)
ite = read.table("results_evaluation/eval_o2o_iter_raw_summary.csv", header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names = 1)


# subset
# ogs_with_changed_root = c("RefOG048","RefOG063","RefOG054","RefOG047","RefOG064","RefOG007","RefOG009","RefOG014","RefOG016","RefOG019","RefOG026","RefOG029","RefOG033","RefOG035","RefOG052","RefOG055","RefOG057")
ogs_with_changed_root = c("RefOG007","RefOG009","RefOG014","RefOG016","RefOG019","RefOG026","RefOG029","RefOG033","RefOG035","RefOG047","RefOG048","RefOG052","RefOG054","RefOG055","RefOG057","RefOG063","RefOG064")
mid = mid [ ogs_with_changed_root , c("precision","recall","Fscore") ]
ite = ite [ ogs_with_changed_root , c("precision","recall","Fscore") ]

# subset
# mid = mid [ rownames(mid) %in% rownames(ite) , c("precision","recall","Fscore") ]
# ite = ite [ rownames(ite) %in% rownames(mid) , c("precision","recall","Fscore") ]

dif = mid 
dif$precision = ite$precision - mid$precision 
dif$recall = ite$recall - mid$recall 
dif$Fscore = ite$Fscore - mid$Fscore
dif = dif[ order(dif$Fscore), ]
mid = mid[ rownames(dif), ]
ite = ite[ rownames(dif), ]

pdf("difference_iter_mid.pdf", height = 6, width = 3)

# diffs
barplot(as.matrix(t(dif)), beside = TRUE, col = c("purple","orange","blue"), border = c("purple4","orange4","blue4"), xlab = "d(iterative-midpoint)", horiz = TRUE, cex.names = 0.5, las=1, xlim = c(-1,1))
legend("bottomright", legend = c("Precision", "Recall","F-score"), fill = c("purple","orange","blue"), border = c("purple4","orange4","blue4"), bty="n")
title("changes in accuracy")

# diff fscore
colnames(mat) = rownames(mid)
barplot(dif$Fscore, col = c("blue"), border = c("blue4"), xlab = "d(iterative-midpoint)", horiz = TRUE, cex.names = 0.5, las=1, xlim = c(-1,1), names.arg = rownames(dif))
title("Fscore change")

# recall
mat = as.matrix(rbind(
	t(ite$recall),
	t(mid$recall)
))
colnames(mat) = rownames(mid)
barplot(mat, beside = TRUE, col = c("springgreen3","gray80"), border = c("springgreen4","gray70"), xlab = "Recall", horiz = TRUE, cex.names = 0.5, las=1, xlim = c(0,1))
legend("bottomright", legend = c("Iterative", "Midpoint"), fill = c("springgreen3","gray80"), border = c("springgreen4","gray70"), bty="n")
title("absolute recall")


# precision
mat = as.matrix(rbind(
	t(ite$precision),
	t(mid$precision)
))
colnames(mat) = rownames(mid)
barplot(mat, beside = TRUE, col = c("springgreen3","gray80"), border = c("springgreen4","gray70"), xlab = "Precision", horiz = TRUE, cex.names = 0.5, las=1, xlim = c(0,1))
legend("bottomright", legend = c("Iterative", "Midpoint"), fill = c("springgreen3","gray80"), border = c("springgreen4","gray70"), bty="n")
title("absolute precision")


# fscore
mat = as.matrix(rbind(
	t(ite$Fscore),
	t(mid$Fscore)
))
colnames(mat) = rownames(mid)
barplot(mat, beside = TRUE, col = c("springgreen3","gray80"), border = c("springgreen4","gray70"), xlab = "F-score", horiz = TRUE, cex.names = 0.5, las=1, xlim = c(0,1))
legend("bottomright", legend = c("Iterative", "Midpoint"), fill = c("springgreen3","gray80"), border = c("springgreen4","gray70"), bty="n")
title("absolute fscore")



dev.off()
