library(gplots)
library("RColorBrewer")

plot_heat_map <- function(data) {
  norm_conf = sweep(data, 2, colSums(data), '/')
  norm_conf[norm_conf < 0] = 0
  heatmap.2(norm_conf,col=rev(heat.colors(100)), scale = "none",trace="none", Rowv=FALSE, Colv=FALSE,dendrogram='none')
}
