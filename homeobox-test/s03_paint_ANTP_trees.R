# libraries
library(stringr)
library(scales)
library(viridis)
library(ape)

# input
tre_fn = "results_annotation/ANTP.possom.ortholog_groups.newick2"
ort_fn = "results_annotation/ANTP.possom.ortholog_groups.csv"
# cnidarians
list_sps_cni = c("Aaur","Adig","Chem","Dgig","Exapal","Fspp","Gasp","Gfas","Hvul","Morvir","Nvec","Spis","Xesp")
list_sps_cni_Regexp = paste(list_sps_cni, collapse = "|")
list_sps_pla = c("Hhon","Tadh")

# read
ort = read.table(ort_fn, header = T, sep="\t", stringsAsFactors = F)

# load
tre = ape::read.tree(tre_fn)
# tre = phangorn::midpoint(tre)
tre = ape::ladderize(tre) 

# add taxonomy data to tips
tre_ort = data.frame(node = tre$tip.label, species = stringr::str_split(tre$tip.label, "_", simplify = T)[,1])
tre_ort = merge(tre_ort, ort, by.x="node", by.y="gene", all.x = T, all.y = F)
rownames(tre_ort) = tre_ort$node
tre_ort = tre_ort[tre$tip.label,]

# list of orthologs to highlight in the tree
list_ort = unique(ort$orthogroup)
# keep only reference OGs
list_ort = list_ort [ !grepl(":like:", list_ort) ]

tax_col = viridis::plasma(length(list_ort),begin = 0.05,end = 0.9)
names(tax_col) = list_ort

tre_ort$color = "gray"
tre_ort$color = tax_col[tre_ort$orthogroup]



# plot phylogenies
# identify cnidarians
ix_is_cnid = grep(list_sps_cni_Regexp,tre$tip.label)
ix_is_nema = grep("^Nvec_",tre$tip.label)
ix_is_othe = which(!grepl(list_sps_cni_Regexp,tre$tip.label))

pdf(sprintf("results_evaluation/tree_phylogram-full.pdf"), height = 120, width = 9)
tre_with_og = tre
tre_with_og$tip.label = paste(tre_with_og$tip.label, str_trunc(tre_ort$orthogroup, 30), sep = " | ")
ape::plot.phylo(
    tre_with_og, edge.color = "darkgray", underscore = T, font=1, 
    col="grey" ,type = "p",
    show.node.label = T, 
    show.tip.label = T, cex=0.5, tip.color = tre_ort$color)
ape::tiplabels(tip=ix_is_cnid, col="cyan3", pch=19, cex=0.5)
ape::tiplabels(tip=ix_is_nema, col="blue", pch=19, cex=0.5)
ape::add.scale.bar(x=0, y=+10, lcol = "darkgray", cex=0.6,length = 0.1)
legend("topright", 
       legend = sprintf("%s", str_trunc(list_ort, width = 30)), 
       col = tax_col[list_ort], pch=19,
       cex=0.5, title = "Species:", bty="n")
dev.off()

pdf(sprintf("results_evaluation/tree_phylogram-small.pdf"), height = 18, width = 6)
tre$tip.label 
ape::plot.phylo(
    tre, edge.color = "darkgray", underscore = T, font=1, 
    col="grey" ,type = "p",
    show.node.label = F, 
    show.tip.label = F, cex=0.5, tip.color = tre_ort$color)
ape::tiplabels(tip=ix_is_cnid, col="cyan3", pch=19, cex=0.5)
ape::tiplabels(tip=ix_is_nema, col="blue", pch=19, cex=0.5)
ape::tiplabels(text=tre$tip.label[ix_is_nema],tip=ix_is_nema, col="blue", cex=0.5, bg = NA, frame = "none", adj = -0.05)
ape::add.scale.bar(lcol = "darkgray", cex=0.6,length = 0.1)

# identify ancestral node for each orthogroup
for (ori in list_ort) {
    
    genes_in_ort = ort[ ort$orthogroup == ori, "gene" ]
    root_of_ort_ix = ape::getMRCA(tre, genes_in_ort)
    root_label = paste(str_trunc(ori,30), unique(ort[ort$orthogroup == ori,"orthogroup_support"]), sep="\n")
    ape::nodelabels("", root_of_ort_ix, col=tax_col[ori], cex=0.4, frame ="none", bg = NA, pch=19)
    ape::nodelabels(root_label, root_of_ort_ix, col=tax_col[ori], cex=0.4, frame ="none", bg = NA)
    
}


dev.off()


# unrooted fll
pdf(sprintf("results_evaluation/tree_unrooted-full.pdf"), height = 80, width = 80)
tre_with_og = tre
tre_with_og$tip.label = paste(tre_with_og$tip.label, str_trunc(tre_ort$orthogroup, 30), sep = " | ")
ape::plot.phylo(
    tre_with_og, edge.color = "darkgray", underscore = T, font=1, 
    col="grey" ,type = "u",
    show.node.label = T, 
    show.tip.label = T, cex=0.5, tip.color = tre_ort$color)
ape::tiplabels(tip=ix_is_cnid, col="cyan3", pch=19, cex=0.5)
ape::tiplabels(tip=ix_is_nema, col="blue", pch=19, cex=0.5)
ape::add.scale.bar(x=0, y=+10, lcol = "darkgray", cex=0.6,length = 0.1)
legend("topright", 
       legend = sprintf("%s", str_trunc(list_ort, width = 30)), 
       col = tax_col[list_ort], pch=19,
       cex=0.5, title = "Species:", bty="n")
dev.off()

# unrooted small, with taxa
pdf(sprintf("results_evaluation/tree_unrooted-small-tax.pdf"), height = 6, width = 6)
tre_with_og = tre
tre_with_og$tip.label = paste(tre_with_og$tip.label, str_trunc(tre_ort$orthogroup, 30), sep = " | ")
ape::plot.phylo(
    tre_with_og, edge.color = "darkgray", underscore = T, font=1, 
    col="grey" ,type = "u",
    show.node.label = F, 
    show.tip.label = F, cex=0.5, tip.color = tre_ort$color)
ape::tiplabels(tip=ix_is_cnid, col="cyan3", pch=19, cex=0.4)
ape::tiplabels(tip=ix_is_nema, col="blue", pch=19, cex=0.4)
ape::tiplabels(tip=ix_is_othe, col="red", pch=19, cex=0.3)
ape::add.scale.bar(lcol = "darkgray", cex=0.6,length = 0.1)
dev.off()


# unrooted small, with OGs
pdf(sprintf("results_evaluation/tree_unrooted-small-ogs.pdf"), height = 6, width = 6)
tre_with_og = tre
tre_with_og$tip.label = rep(".",length(tre$tip.label))
ape::plot.phylo(
    tre_with_og, edge.color = "darkgray", underscore = T, font=1, 
    col="grey" ,type = "f",
    show.node.label = F, align.tip.label = T,
    show.tip.label = T, cex=0.5, tip.color = tre_ort$color)
ape::tiplabels("", pch=19, col = tre_ort$color, cex=0.5, bg = NA, frame = "none")
ape::add.scale.bar(lcol = "darkgray", cex=0.6,length = 0.1)
dev.off()


#### summary OG annotations ####

pdf(sprintf("results_evaluation/summary_annotations.pdf"), height = 3, width = 3)

# first for Nvec
ort_nvec = ort[grepl("^Nvec",ort$gene), ]
ort_nvec$annotation_string = NA

ort_nvec [ grepl(":like:",ort_nvec$orthogroup) , "annotation_string"] = "No reference"
ort_nvec [ !grepl(":like:",ort_nvec$orthogroup) & !grepl("/", ort_nvec$orthogroup ) , "annotation_string"] = "Reference, single human ortholog"
ort_nvec [ !grepl(":like:",ort_nvec$orthogroup) & grepl("/", ort_nvec$orthogroup ) , "annotation_string"] = "Reference, multiple human orthologs"
ort_nvec_xtab = table(ort_nvec$annotation_string)

par(mar=c(5.1,6,4.1,2.1))
pie(ort_nvec_xtab, 
        labels = paste(names(ort_nvec_xtab), ort_nvec_xtab), las=1, horiz = T,
        xlab="Frequency", cex = 0.6)
title("annot Nvec",
      sub=sum(ort_nvec_xtab))

# second for all cnidarians
ort$species = stringr::str_split(ort$gene, "_", simplify = T)[,1]
ort_nvec = ort[ort$species %in% list_sps_cni, ]
ort_nvec$annotation_string = NA

ort_nvec [ grepl(":like:",ort_nvec$orthogroup) , "annotation_string"] = "No reference"
ort_nvec [ !grepl(":like:",ort_nvec$orthogroup) & !grepl("/", ort_nvec$orthogroup ) , "annotation_string"] = "Reference, single human ortholog"
ort_nvec [ !grepl(":like:",ort_nvec$orthogroup) & grepl("/", ort_nvec$orthogroup ) , "annotation_string"] = "Reference, multiple human orthologs"
ort_nvec_xtab = table(ort_nvec$annotation_string)

par(mar=c(5.1,6,4.1,2.1))
pie(ort_nvec_xtab, 
    labels = paste(names(ort_nvec_xtab), ort_nvec_xtab), las=1, horiz = T,
    xlab="Frequency", cex = 0.6)
title("annot cnidaria",
      sub=sum(ort_nvec_xtab))


dev.off()

#### evaluate nvec pairs ####

orp_fn = "results_annotation/ANTP.possom.pairs_all.csv"

orp = read.table(orp_fn, header = T)
orp_nvec = orp [ grepl("^Nvec_",orp$in_gene) | grepl("^Nvec_",orp$out_gene), ]
orp_nvec = orp_nvec [ orp_nvec$ev_type != "outparalog_ext",]



pdf(sprintf("results_evaluation/summary_relationships.pdf"), height = 3, width = 3)

orp_nvec_xtab = table(orp[  grepl("Nvec_v1g192469", orp$in_gene) |  grepl("Nvec_v1g192469", orp$out_gene) , ]$ev_type)
pie(orp_nvec_xtab, 
    labels = paste(names(orp_nvec_xtab), orp_nvec_xtab), las=1, horiz = T,
    xlab="Frequency", cex = 0.6)
title("Evolutionary relationships to Nvec_v1g192469", cex.main=0.5)

orp_nvec_xtab = orp_nvec_xtab[names(orp_nvec_xtab) != "outparalog_ext"]
pie(orp_nvec_xtab, 
    labels = paste(names(orp_nvec_xtab), orp_nvec_xtab), las=1, horiz = T,
    xlab="Frequency", cex = 0.6)
title("Evolutionary relationships to Nvec_v1g192469", cex.main=0.5)
dev.off()
