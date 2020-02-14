library(gplots)
library("RColorBrewer")

plot_heat_map <- function(data) {
  norm_conf = sweep(data, 2, colSums(data), '/')
  heatmap.2(norm_conf,col=brewer.pal(11,"RdBu"), trace="none", Rowv=FALSE, Colv=FALSE,dendrogram='none')
}
