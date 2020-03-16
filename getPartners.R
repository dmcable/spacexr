library(RCTD)
iv <- init_RCTD(gene_list_reg = T, get_proportions = F) #initial variables
#get dendrogram scratch
distance_mat <- Matrix(0, nrow = iv$cell_type_info[[3]], ncol = iv$cell_type_info[[3]])
rownames(distance_mat) <- iv$cell_type_info[[2]]; colnames(distance_mat) <- iv$cell_type_info[[2]]
epsilon = 1e-9
for(type1 in iv$cell_type_info[[2]])
  for(type2 in iv$cell_type_info[[2]]) {
    diff = log(iv$cell_type_info[[1]][iv$gene_list,type1] + epsilon) - log(iv$cell_type_info[[1]][iv$gene_list,type2] + epsilon)
    distance_mat[type1,type2] <- mean(diff^2)
    if(type1 == type2)
      distance_mat[type1, type2] = 1000
  }
partner <- iv$cell_type_info[[2]][apply(distance_mat,1,which.min)]
names(partner) = iv$cell_type_info[[2]]
pairs = list()
for(type in iv$cell_type_info[[2]]) {
  if(partner[partner[type]] == type)
    pairs = rbind(pairs, c(type, partner[type]))
}
