# libraries
library(igraph)

# input
clu_fn = "results_trees/TALE.possom.ortholog_groups.csv"
grd_fn = "results_trees/TALE.possom.pairs_orthologs.csv"

# read data
clu = read.table(clu_fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(clu) = clu$gene
grd = read.table(grd_fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
grd = grd[,1:2]

# create graph
gra = igraph::graph.data.frame(grd, directed = FALSE)

# get communities
list_vertices = V(gra)$name
list_clusters = clu[ list_vertices, "orthogroup" ]
list_clusters = as.numeric(as.factor(list_clusters))

# get colors for communities
color_palette_fun = colorRampPalette(c("magenta4","firebrick1","orange","khaki1","springgreen2","darkgreen","deepskyblue","cadetblue1","mediumblue","darkviolet","violet"))
num_clusters = length(sort(unique(list_clusters)))
color_palette = color_palette_fun(num_clusters)
names(color_palette) = unique(sort(list_clusters))
list_colors = color_palette [ list_clusters ]

# add colors to grap
V(gra)$color = list_colors

# plot
gra_layout = layout_with_fr(gra)
igraph::plot.igraph(gra, layout = gra_layout, vertex.size=2, vertex.frame.color = NA, vertex.label = NA, vertex.color = list_colors)



