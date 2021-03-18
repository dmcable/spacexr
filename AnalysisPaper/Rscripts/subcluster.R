refdir <- 'Data/Reference/DropVizHC'
reference <- readRDS('Data/Reference/DropVizHC/scRef.RDS')
inter_ident <- which(reference@meta.data$liger_ident_coarse == "Interneuron")
inter_ident <- rownames(reference@meta.data)[inter_ident]
counts <- reference@assays$RNA@counts[,inter_ident]
nUMI <- reference@meta.data[inter_ident,"nUMI"]; names(nUMI) <- inter_ident
subcluster <- readRDS('Data/Reference/DropVizHC/Subcluster/F_GRCm38.81.P60Hippocampus.cell_cluster_outcomes.RDS')
inter_ident <- inter_ident[subcluster[inter_ident,]$subcluster != "1"]
table(subcluster[inter_ident,]$subcluster)


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


cell_dict_file <- file.path(refdir,"Subcluster/subclusterlabels.csv")
true_type_names <- remap_celltypes(cell_dict_file, subcluster[inter_ident,]$subcluster)
meta_data = as.data.frame(subcluster[inter_ident,]$subcluster)
meta_data$nUMI = nUMI[inter_ident]
meta_data$liger_ident_coarse = true_type_names
rownames(meta_data) = inter_ident
colnames(meta_data)[1] = 'cluster'
reference_interneuron = Seurat::CreateSeuratObject(counts[,inter_ident], meta.data = meta_data)
cell_type_info_inter_all <- get_cell_type_info(reference_interneuron@assays$RNA@counts, reference_interneuron@meta.data$liger_ident_coarse, reference_interneuron@meta.data$nUMI)
saveRDS(reference_interneuron, file.path(refdir,'interneurons.RDS'))

cell_dict_file <- file.path(refdir,"Subcluster/subclusterlabels_coarse.csv")
true_type_names <- remap_celltypes(cell_dict_file, meta_data$cluster)
meta_data$liger_ident_coarse = true_type_names
reference_interneuron_coarse = Seurat::CreateSeuratObject(counts[,inter_ident], meta.data = meta_data)
cell_type_info_inter <- get_cell_type_info(reference_interneuron_coarse@assays$RNA@counts, reference_interneuron_coarse@meta.data$liger_ident_coarse, reference_interneuron_coarse@meta.data$nUMI)
iv <- init_RCTD(gene_list_reg = T, get_proportions = F) #initial variables

gene_list_inter <- iv$gene_list

#run the gather_results for main
iv2 <- init_RCTD(gene_list_reg = T, get_proportions = T)
puck <- readRDS(file.path(iv$slideseqdir, config_data$puckrds))
gene_list <- union(iv2$gene_list, gene_list_inter)
puck = restrict_counts(puck, gene_list, UMI_thresh = iv2$config$UMI_min, UMI_max = iv2$config$UMI_max)
bulk_vec = rowSums(puck@counts)
weight_avg = rowSums(sweep(iv2$cell_type_info[[1]][gene_list,],2,iv2$proportions,'*'))
target_means = bulk_vec[gene_list]/sum(puck@nUMI)
cell_type_means_merged = cbind(iv2$cell_type_info[[1]][gene_list,iv2$cell_type_info[[2]] != "Interneuron"], cell_type_info_inter[[1]][gene_list,])
cell_type_means_renorm = sweep(cell_type_means_merged[gene_list,],1,weight_avg / target_means,'/')
cell_type_names <- c(iv2$cell_type_info[[2]][iv2$cell_type_info[[2]] != "Interneuron"], cell_type_info_inter[[2]])
cell_type_info <- list(cell_type_means_renorm, cell_type_names, length(cell_type_names))
saveRDS(cell_type_info,file.path(refdir,'info_renorm_coarse.RDS'))

#all interneurons
puck <- readRDS(file.path(iv$slideseqdir, config_data$puckrds))
fine_list <- get_de_genes(cell_type_info_inter_all, puck, fc_thresh = 0.25, expr_thresh = .000125, MIN_OBS = 3)
gene_list <- union(iv2$gene_list, fine_list)
puck = restrict_counts(puck, gene_list, UMI_thresh = iv2$config$UMI_min, UMI_max = iv2$config$UMI_max)
bulk_vec = rowSums(puck@counts)
weight_avg = rowSums(sweep(iv2$cell_type_info[[1]][gene_list,],2,iv2$proportions,'*'))
target_means = bulk_vec[gene_list]/sum(puck@nUMI)
cell_type_means_merged = cbind(iv2$cell_type_info[[1]][gene_list,iv2$cell_type_info[[2]] != "Interneuron"], cell_type_info_inter_all[[1]][gene_list,])
cell_type_means_renorm = sweep(cell_type_means_merged[gene_list,],1,weight_avg / target_means,'/')
cell_type_names <- c(iv2$cell_type_info[[2]][iv2$cell_type_info[[2]] != "Interneuron"], cell_type_info_inter_all[[2]])
cell_type_info <- list(cell_type_means_renorm, cell_type_names, length(cell_type_names))
saveRDS(cell_type_info,file.path(refdir,'info_renorm_all.RDS'))
