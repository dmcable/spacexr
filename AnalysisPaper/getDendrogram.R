library(RCTD)
#gets dendrogram for interneurons
iv <- init_RCTD(gene_list_reg = T, get_proportions = F) #initial variables
#get dendrogram scratch
distance_mat <- Matrix(0, nrow = iv$cell_type_info[[3]], ncol = iv$cell_type_info[[3]])
rownames(distance_mat) <- iv$cell_type_info[[2]]; colnames(distance_mat) <- iv$cell_type_info[[2]]
epsilon = 1e-9
for(type1 in iv$cell_type_info[[2]])
  for(type2 in iv$cell_type_info[[2]]) {
    diff = log(iv$cell_type_info[[1]][iv$gene_list,type1] + epsilon) - log(iv$cell_type_info[[1]][iv$gene_list,type2] + epsilon)
    distance_mat[type1,type2] <- sqrt(mean(diff^2))
    if(type1 == type2)
      distance_mat[type1, type2] = 0
  }
plot_heat_map(pmin(as.matrix(1/distance_mat),1)*1, normalize = F)
norm_conf <- pmin(as.matrix(1/distance_mat),1)
heatmap.2(as.matrix(norm_conf[hc1$order,hc1$order]),col=rev(heat.colors(100)), breaks=seq(0,1,0.01), scale = "none",trace="none")
plot_heat_map(as.matrix(norm_conf[hc1$order,hc1$order]), normalize = F)
heatmap(pmin(as.matrix(1/distance_mat),0.5))
diag(distance_mat) <- 10
partner <- iv$cell_type_info[[2]][apply(distance_mat,1,which.min)]
names(partner) = iv$cell_type_info[[2]]
pairs = list()
for(type in iv$cell_type_info[[2]]) {
  if(partner[partner[type]] == type)
    pairs = rbind(pairs, c(type, partner[type]))
}
d <- dist(log(t(iv$cell_type_info[[1]][iv$gene_list,]) + epsilon), method = "euclidean")
hc1 <- hclust(d, method = "ward.D" )
plot(hc1, cex = 0.6, hang = -1)
saveRDS(d,file.path(refdir,"Subcluster/subclusterdist.RDS"))
saveRDS(hc1,file.path(refdir,"Subcluster/dendrogram.RDS"))
