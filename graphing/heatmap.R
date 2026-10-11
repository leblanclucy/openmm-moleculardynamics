#for making the heatmap of the similarity matrix of amino acid sequences
#does clustering with euclidean distance and average linkage
#cuts rows into 8 distinct clusters to show the separation more cleanly

install.packages("RColorBrewer")
install.packages("pheatmap")
library(RColorBrewer)
library(pheatmap)

colpal <- colorRampPalette(brewer.pal(8, "RdBu"))(100)
chitinases <- read.csv("percentidentity_chitinases.csv")
data <- as.matrix(chitinases[, -1])
anno <- read.csv("percentidentity_chitinases_anno.csv")

anno_row <- data.frame(Family = factor(anno$Family,  levels = c("GH18", "GH19", "Unknown"), labels = c("GH18", "GH19", "UnknownF")), Class = factor(anno$Class, levels = c("I", "II", "III", "IV", "V", "Unknown"), labels = c("I", "II", "III", "IV", "V", "UnknownC")))

rownames(data) <- rownames(anno_row)

annotation_colors <- list(Family = c(GH18 = "lightgoldenrod1", GH19 = "lightskyblue1", UnknownF = "gray"), Class  = c(I = "turquoise1", II = "darkorchid3", III = "red", IV = "dodgerblue3", V = "darkorange", UnknownC = "gray"))

pheatmap(data, color=colpal, show_colnames=FALSE, show_rownames=FALSE, cutree_rows=8, clustering_method="average", annotation_row=anno_row, annotation_colors=annotation_colors, border_color=NA)
