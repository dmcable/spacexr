#interneuron
refdir <- 'Data/Reference/DropVizHC'
reference <- readRDS('Data/Reference/DropVizHC/scRef.RDS')
inter_ident <- which(reference@meta.data$liger_ident_coarse == "Interneuron")
inter_ident <- rownames(reference@meta.data)[inter_ident]
counts <- reference@assays$RNA@counts[,inter_ident]
nUMI <- reference@meta.data[inter_ident,"nUMI"]; names(nUMI) <- inter_ident
subcluster <- readRDS('Data/Reference/DropVizHC/Subcluster/F_GRCm38.81.P60Hippocampus.cell_cluster_outcomes.RDS')
inter_ident <- inter_ident[subcluster[inter_ident,]$subcluster != "1"]
table(subcluster[inter_ident,]$subcluster)

meta_data = as.data.frame(cluster)
meta_data$nUMI = colSums(raw.data)
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
iv <- init_RCTD(gene_list_reg = F, get_proportions = F) #initial variables
table(reference_interneuron@meta.data$liger_ident_coarse)

#after get dendrogram:
cell_dict_file <- file.path(refdir,"Subcluster/subclusterlabels_coarse.csv")
true_type_names <- remap_celltypes(cell_dict_file, meta_data$cluster)
meta_data$liger_ident_coarse = true_type_names
reference_interneuron_coarse = Seurat::CreateSeuratObject(counts[,inter_ident], meta.data = meta_data)
saveRDS(reference_interneuron_coarse, file.path(refdir,'interneurons_coarse.RDS'))
cell_type_info_inter <- get_cell_type_info(reference_interneuron_coarse@assays$RNA@counts, reference_interneuron_coarse@meta.data$liger_ident_coarse, reference_interneuron_coarse@meta.data$nUMI)
iv <- init_RCTD(gene_list_reg = T, get_proportions = F) #initial variables
saveRDS(iv$gene_list, file.path(refdir,'gene_names_coarse.RDS'))
my_mat <- iv$cell_type_info[[1]][iv$gene_list,]
my_mat$log_fc <- log(apply(my_mat,1,max)) - log(3/2*(rowMeans(my_mat) - 1/3*apply(my_mat,1,max)) + 1e-10)
my_mat$max_type <- apply(my_mat[],1,which.max)
tail(my_mat[order(my_mat$Basket_OLM),])
tail(my_mat[order(my_mat$Neurogliaform_Lacunosum),], 100)
Basket_OLM: 'Sst'
CGE: 'Vip'
Neurogliaform_Lacunosum: 'Id2'
gene_list_inter <- iv$gene_list

#run the gatherDoubletResults for main
#start with singlets
iv2 <- init_RCTD(gene_list_reg = T, get_proportions = T)
cell_type_means_renorm = iv2$cell_type_info
#cell_type_info_renorm[[1]] = get_norm_ref(iv$puck, iv2$cell_type_info[[1]], union(gene_list_inter,iv2$gene_list), proportions)
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
Q_mat <- readRDS(file.path(resultsdir,'Q_mat.RDS'))
N_X = dim(Q_mat)[2]; delta = 1e-5; X_vals = (1:N_X)^1.5*delta
K_val = dim(Q_mat)[1] - 3; use_Q = T
inter_names<- c('CGE', "Basket_OLM" , "Neurogliaform_Lacunosum")
log_l_thresh <- 10
singlet_ind = results_df$first_type == "Interneuron" & results_df$spot_class == "singlet"
singlet_barcodes <- barcodes[singlet_ind]
doublet_barcodes <- c(barcodes[results_df$first_type == "Interneuron" & results_df$spot_class == "doublet_certain"], barcodes[results_df$first_type == "Interneuron" & results_df$spot_class == "doublet_uncertain"],
                      barcodes[results_df$second_type == "Interneuron" & results_df$spot_class == "doublet_certain"])
doub_first <- c(barcodes[results_df$first_type == "Interneuron" & results_df$spot_class == "doublet_certain"], barcodes[results_df$first_type == "Interneuron" & results_df$spot_class == "doublet_uncertain"])
doub_second <- barcodes[results_df$second_type == "Interneuron" & results_df$spot_class == "doublet_certain"]
second_type_list <- unlist(list(results_df[doub_first,]$second_type,results_df[doub_second,]$first_type))
names(second_type_list) <- doublet_barcodes
inter_barcodes <- c(singlet_barcodes, doublet_barcodes)
empty_cell_types = factor(character(N),levels = cell_type_names)
N = length(inter_barcodes)
inter_df <- data.frame(best_type = factor(character(N),levels = inter_names), confident = logical(N), score_diff = numeric(N))
rownames(inter_df) <- inter_barcodes
for(barcode in singlet_barcodes) {
  score_best <- 100000
  score_second <- 100000
  best_type <- NULL
  for (type in inter_names) {
    score <- get_singlet_score(cell_type_info, gene_list, puck@counts[gene_list,barcode], puck@nUMI[barcode], type, F)
    if(score < score_best) {
      score_second <- score_best
      score_best <- score
      best_type <- type
    } else if(score < score_second) {
      score_second <- score
    }
    inter_df[barcode,type] <- score
  }
  inter_df[barcode,"confident"] <- (score_second - score_best) > log_l_thresh
  inter_df[barcode,"score_diff"] <- (score_second - score_best)
  inter_df[barcode,"best_type"] <- best_type
}

for(barcode in doublet_barcodes) {
  score_best <- 100000
  score_second <- 100000
  best_type <- NULL
  for (type in inter_names) {
    score <- decompose_sparse(cell_type_info[[1]], gene_list, puck@nUMI[barcode], puck@counts[gene_list,barcode], type1=type, type2=as.character(second_type_list[barcode]), score_mode = T, constrain = F)
    if(score < score_best) {
      score_second <- score_best
      score_best <- score
      best_type <- type
    } else if(score < score_second) {
      score_second <- score
    }
    inter_df[barcode,type] <- score
  }
  inter_df[barcode,"confident"] <- (score_second - score_best) > log_l_thresh
  inter_df[barcode,"score_diff"] <- (score_second - score_best)
  inter_df[barcode,"best_type"] <- best_type
}

#doublets
plot_weights_doublet(iv$cell_type_info, puck, resultsdir, weights_doublet, results_df)
cell_type = cell_type_info[[2]][i]
all_weights <- weights_doublet[results_df$spot_class == "doublet_certain" & results_df$second_type == cell_type,2, drop=FALSE]
all_weights <- rbind(all_weights, weights_doublet[!(results_df$spot_class == "reject") & results_df$first_type == cell_type,1,drop=FALSE])
all_weights_vec <- as.vector(all_weights); names(all_weights_vec) <- rownames(all_weights)
if(length(all_weights) > 0)
  plots[[i]] <- plot_puck_continuous(puck, rownames(all_weights), all_weights_vec, title = cell_type, ylimit = c(0,1))

cell_type = "Interneuron"
all_weights <- weights_doublet[results_df$spot_class == "doublet_certain" & results_df$second_type == cell_type,2, drop=FALSE]
all_weights <- rbind(all_weights, weights_doublet[!(results_df$spot_class == "reject") & results_df$first_type == cell_type,1,drop=FALSE])
all_weights_vec <- as.vector(all_weights); names(all_weights_vec) <- rownames(all_weights)
all_weights_type <- all_weights_vec[inter_barcodes[inter_df$confident & inter_df$best_type==inter_names[3]]]
plot_puck_continuous(puck, names(all_weights_type), all_weights_type, title = inter_names[3], ylimit = c(0,1))


#plot all
barcodes_cur = inter_barcodes[inter_df$confident]
barcodes_cur = inter_barcodes
my_class <- inter_df[["best_type"]]
names(my_class) = inter_barcodes
plot_class(puck, barcodes_cur, my_class)
plot_class <- function(puck, barcodes_cur, my_class, counter_barcodes = NULL, title = NULL) {
  my_table = puck@coords[barcodes_cur,]
  my_table$class = my_class[barcodes_cur]
  n_levels <- length(unique(my_class))
  if(n_levels > 36)
    cols = rainbow(n_levels)[sample(1:n_levels, n_levels, replace = F)]
  else {
    if(n_levels > 21)
      cols = unname(pals::polychrome(n_levels))[3:(n_levels+2)]
    else
      cols = pals::kelly(n_levels+1)[3:(n_levels+2)]
  }
  plot <- ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) + ggplot2::geom_point(ggplot2::aes(size = 0.5, shape=19,color=class)) +
    ggplot2::scale_shape_identity() + ggplot2::theme_bw() + ggplot2::scale_size_identity() + scale_color_manual(values = cols)
  if(!is.null(counter_barcodes)) {
    my_table = puck@coords[counter_barcodes,]
    plot <- plot + ggplot2::geom_point(data=my_table,aes(x=x, y=y,size=0.15),alpha = 0.1)
  }
  if(!is.null(title))
    plot <- plot + ggplot2::ggtitle(title)
  plot
}


plot_puck_wrapper(iv$puck, iv$puck@counts['Vip',], cell_type = "Interneuron",max_val = 3, maxUMI = 200000)
plot_puck_wrapper(iv$puck, iv$puck@counts['Vip',],max_val = 3, maxUMI = 200000)

probabilities <- inter_df[,inter_names]
probabilities<- exp(-sweep(probabilities,1, apply(probabilities,1,min),'-')/10)
probabilities <- sweep(probabilities,1,rowSums(probabilities),'/')
plot_val <- probabilities[,"CGE"]; names(plot_val) <- inter_barcodes
plot_puck_continuous(iv$puck, inter_barcodes, plot_val, title = "CGE", ylimit = c(0,1))
plot_val <- probabilities[,inter_names[2]]; names(plot_val) <- inter_barcodes
plot_puck_continuous(iv$puck, inter_barcodes, plot_val, title = inter_names[2], ylimit = c(0,1))
counter_barcodes <- barcodes[results_df$spot_class == "singlet" & results_df$first_type %in% c("CA1","CA3","Denate")]
plot_val <- probabilities[,inter_names[1]]; names(plot_val) <- inter_barcodes
plot_puck_continuous(iv$puck, inter_barcodes, plot_val, title = inter_names[1], ylimit = c(0,1), counter_barcodes = counter_barcodes)
saveRDS(probabilities,file.path(resultsdir, 'probabilities_inter_coarse.RDS'))
saveRDS(results_df,file.path(resultsdir,'results_df.RDS'))
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

#
refdir <- 'Data/Reference/DropVizHC'
cell_type_info <- readRDS(file.path(refdir,'info_renorm_all.RDS'))
gene_list <- rownames(cell_type_info[[1]])
puck <- readRDS(file.path(iv$slideseqdir, iv$puckrds))
puck = restrict_counts(puck, gene_list, UMI_thresh = 100, UMI_max = 200000)
Q_mat <- readRDS(file.path(resultsdir,'Q_mat.RDS'))
N_X = dim(Q_mat)[2]; delta = 1e-5; X_vals = (1:N_X)^1.5*delta
K_val = dim(Q_mat)[1] - 3; use_Q = T
log_l_thresh <- 10
singlet_ind = results_df$first_type == "Interneuron" & results_df$spot_class == "singlet"
singlet_barcodes <- barcodes[singlet_ind]
doublet_barcodes <- c(barcodes[results_df$first_type == "Interneuron" & results_df$spot_class == "doublet_certain"], barcodes[results_df$first_type == "Interneuron" & results_df$spot_class == "doublet_uncertain"],
                      barcodes[results_df$second_type == "Interneuron" & results_df$spot_class == "doublet_certain"])
doub_first <- c(barcodes[results_df$first_type == "Interneuron" & results_df$spot_class == "doublet_certain"], barcodes[results_df$first_type == "Interneuron" & results_df$spot_class == "doublet_uncertain"])
doub_second <- barcodes[results_df$second_type == "Interneuron" & results_df$spot_class == "doublet_certain"]
second_type_list <- unlist(list(results_df[doub_first,]$second_type,results_df[doub_second,]$first_type))
names(second_type_list) <- doublet_barcodes
inter_barcodes <- c(singlet_barcodes, doublet_barcodes)
empty_cell_types = factor(character(N),levels = cell_type_names)
N = length(inter_barcodes)
inter_names <- cell_type_info[[2]][17:43]
inter_df <- data.frame(best_type = factor(character(N),levels = inter_names), confident = logical(N), score_diff = numeric(N))
rownames(inter_df) <- inter_barcodes
for(barcode in singlet_barcodes[1]) {
  score_best <- 100000
  score_second <- 100000
  best_type <- NULL
  for (type in inter_names) {
    score <- get_singlet_score(cell_type_info, gene_list, puck@counts[gene_list,barcode], puck@nUMI[barcode], type, F)
    if(score < score_best) {
      score_second <- score_best
      score_best <- score
      best_type <- type
    } else if(score < score_second) {
      score_second <- score
    }
    inter_df[barcode,type] <- score
  }
  inter_df[barcode,"confident"] <- (score_second - score_best) > log_l_thresh
  inter_df[barcode,"score_diff"] <- (score_second - score_best)
  inter_df[barcode,"best_type"] <- best_type
}

for(barcode in doublet_barcodes) {
  score_best <- 100000
  score_second <- 100000
  best_type <- NULL
  for (type in inter_names) {
    score <- decompose_sparse(cell_type_info[[1]], gene_list, puck@nUMI[barcode], puck@counts[gene_list,barcode], type1=type, type2=as.character(second_type_list[barcode]), score_mode = T, constrain = F)
    if(score < score_best) {
      score_second <- score_best
      score_best <- score
      best_type <- type
    } else if(score < score_second) {
      score_second <- score
    }
    inter_df[barcode,type] <- score
  }
  inter_df[barcode,"confident"] <- (score_second - score_best) > log_l_thresh
  inter_df[barcode,"score_diff"] <- (score_second - score_best)
  inter_df[barcode,"best_type"] <- best_type
}
saveRDS(inter_df, file.path(resultsdir,'inter_df_process.RDS'))
inter_df[singlet_barcodes,] <- job_script_results$inter_df[singlet_barcodes,]
inter_df[doublet_barcodes,] <- job_script2_results$inter_df[doublet_barcodes,]
probabilities <- inter_df[,inter_names]
probabilities<- exp(-sweep(probabilities,1, apply(probabilities,1,min),'-')/10)
probabilities <- sweep(probabilities,1,rowSums(probabilities),'/')
type <- "CA1_lacunosum"
plot_val <- probabilities[,type]; names(plot_val) <- inter_barcodes
plot_puck_continuous(iv$puck, inter_barcodes, plot_val, title = type, ylimit = c(0,1))
counter_barcodes <- barcodes[results_df$spot_class == "singlet" & results_df$first_type %in% c("CA1","CA3","Denate")]
ind = 8
plot_val <- probabilities[,inter_names[ind]]; names(plot_val) <- inter_barcodes
plot_puck_continuous(iv$puck, inter_barcodes, plot_val, title = inter_names[ind], ylimit = c(0,1), counter_barcodes = counter_barcodes)
hc1 <- readRDS(file.path(refdir,"Subcluster/dendrogram.RDS"))
cur_types_sub = list(c('CGE_1','CGE_2','CGE_3'),c('CGE_4','CGE_5','CGE_6','CGE_7','CGE_8','CGE_9','CGE_10','CGE_11','CGE_12'),inter_names[c(3,18)], inter_names[19:22],inter_names[c(2,16)],inter_names[c(1,17,23:27)])
cur_types = c('CGE_1','CGE_2','CGE_3','CGE_4','CGE_5','CGE_6','CGE_7','CGE_8','CGE_9','CGE_10','CGE_11','CGE_12')
plot_val <- rowSums(probabilities[,cur_types])
cur_types2 <- c(inter_names[3],inter_names[18:22])
plot_val2 <- rowSums(probabilities[,cur_types2])
plot_puck_continuous(iv$puck, inter_barcodes, plot_val, title = inter_names[ind], ylimit = c(0,1), counter_barcodes = counter_barcodes)
pdf(file.path(resultsdir,'clusters.pdf'))
for(i in (1:num_clusters)) {
  if(sum(my_class==i) > 1) {
    print(plot_puck_continuous(iv$puck, inter_barcodes[my_class==i], plot_val, title = i, ylimit = c(0,1), counter_barcodes = counter_barcodes))
    print(plot_puck_continuous(iv$puck, inter_barcodes[my_class==i], plot_val2, title = i, ylimit = c(0,1), counter_barcodes = counter_barcodes))
    print(plot_puck_continuous(iv$puck, inter_barcodes[my_class==i], 1 - plot_val -plot_val2, title = i, ylimit = c(0,1), counter_barcodes = counter_barcodes))
  }
}
dev.off()
pdf(file.path(resultsdir,'clusters_sub.pdf'))
for(i in c(25,27,38,47,52,53,61,65,72,77,86,98,140,144,169,175)) {
  for(cur_types_inst in cur_types_sub) {
    plot_val_sub <- rowSums(probabilities[,cur_types_inst])
    print(plot_puck_continuous(iv$puck, inter_barcodes[my_class==i], plot_val_sub, title = i, ylimit = c(0,1), counter_barcodes = counter_barcodes))
  }
}
dev.off()
pdf(file.path(resultsdir,'clusters_sub.pdf'))
for(i in c(25,27,38,47,52,53,61,65,72,77,86,98,140,144,169,175)) {
  for(cur_types_inst in cur_types_sub) {
    plot_val_sub <- rowSums(probabilities[,cur_types_inst])
    print(plot_puck_continuous(iv$puck, inter_barcodes[my_class==i], plot_val_sub, title = i, ylimit = c(0,1), counter_barcodes = counter_barcodes))
  }
}
dev.off()
pdf(file.path(resultsdir,'clusters_label.pdf'))
for(i in c(7,12,14,20,22,23,29,39,42,43,60,62,63,68,71,88,96,127,142,162)) {
  print(plot_puck_continuous(iv$puck, inter_barcodes[my_class==i], plot_val, title = i, ylimit = c(0,1), label=T))
  print(plot_puck_continuous(iv$puck, inter_barcodes[my_class==i], plot_val2, title = i, ylimit = c(0,1), label=T))
  print(plot_puck_continuous(iv$puck, inter_barcodes[my_class==i], 1 - plot_val -plot_val2, title = i, ylimit = c(0,1), label=T))
}
dev.off()
pdf(file.path(resultsdir,'clusters_sub_label.pdf'))
for(i in c(25,47,52,61,77,86,175)) {
  for(cur_types_inst in cur_types_sub) {
    plot_val_sub <- rowSums(probabilities[,cur_types_inst])
    print(plot_puck_continuous(iv$puck, inter_barcodes[my_class==i], plot_val_sub, title = i, ylimit = c(0,1),label=T))
  }
}
dev.off()
library(forcats)
d <- dist(puck@coords[inter_barcodes,], method = "euclidean")
hc1 <- hclust(d, method = "average")
plot(hc1, cex = 0.6, hang = -1)
num_clusters = 200
cutree(hc1,k=num_clusters)

barcodes_cur = inter_barcodes
my_class <- as.factor(cutree(hc1,k=200))
my_class <- fct_shuffle(my_class)
plot_class(puck, barcodes_cur, my_class)

my_class
relabel <- read.csv(file.path(iv$slideseqdir,'cluster_relabels.csv'))
#validate
for(cl_label in unique(relabel$Cluster)) {
  print(cl_label)
  set1 = inter_barcodes[my_class==cl_label]
  set2 = barcodes[relabel[relabel$Cluster==cl_label,"Index"]]
  print(length(set1) == length(set2) && length(set1) == length(intersect(set1,set2)))
  #which(barcodes[relabel[relabel$Cluster==cl_label,"Index"]] == setdiff(set2,set1))
}

for(cl_label in unique(relabel$Cluster)) {
  print(cl_label)
  set1 = inter_barcodes[my_class==cl_label]
  set2 = barcodes[relabel[relabel$Cluster==cl_label,"Index"]]
  print(length(set1) == length(set2) && length(set1) == length(intersect(set1,set2)))
  #which(barcodes[relabel[relabel$Cluster==cl_label,"Index"]] == setdiff(set2,set1))
}
#fix dumb mistake
for(i in unique(relabel$Cluster))
  relabel[relabel$Cluster==i,"barcodes"] <- mapvalues(relabel[relabel$Cluster==i,]$Index,from = which(rownames(puck@coords) %in% inter_barcodes[my_class==i]), to=inter_barcodes[my_class==i])
rownames(relabel) <- relabel$barcodes
new_labels <- as.character(my_class)
names(new_labels) <- names(my_class)
new_labels[relabel$barcode] <- apply(relabel,1,function(x) paste(x[3],x[2],sep='_'))
new_labels <- as.factor(new_labels)




pdf(file.path(resultsdir,'clusters_final.pdf'))
for(i in levels(new_labels)) {
  if(sum(new_labels==i) > 3) {
    print(plot_puck_continuous(iv$puck, inter_barcodes[new_labels==i], plot_val, title = i, ylimit = c(0,1), counter_barcodes = counter_barcodes))
    print(plot_puck_continuous(iv$puck, inter_barcodes[new_labels==i], plot_val2, title = i, ylimit = c(0,1), counter_barcodes = counter_barcodes))
    print(plot_puck_continuous(iv$puck, inter_barcodes[new_labels==i], 1 - plot_val -plot_val2, title = i, ylimit = c(0,1), counter_barcodes = counter_barcodes))
  }
}
dev.off()
saveRDS(new_labels, file.path(refdir,'Subcluster/new_labels.RDS'))
library(plyr)
new_labels <- mapvalues(new_labels, from = levels(new_labels), to = sample(1:length(levels(new_labels))))
plot_class(puck, sample(inter_barcodes), new_labels)

probabilities <- inter_df[,inter_names]
probabilities <- exp(-sweep(probabilities,1, apply(probabilities,1,min),'-')/10)
probabilities <- sweep(probabilities,1,rowSums(probabilities),'/')

final_labels = factor(character(length(inter_barcodes)),levels = inter_names)
names(final_labels) = inter_barcodes
n_classes = length(levels(new_labels))
for(i in(1:n_classes)) {
  final_labels[new_labels==i] <- names(which.max(colMeans(probabilities[inter_barcodes[new_labels==i],])))
}
counter_barcodes <- barcodes[results_df$spot_class == "singlet" & results_df$first_type %in% c("CA1","CA3","Denate")]

plot_class(puck, inter_barcodes, final_labels, counter_barcodes = counter_barcodes)
save(my_en)

hc1 <- readRDS(file.path(refdir,"Subcluster/dendrogram.RDS"))
cutree(hc1,k=6)
final_labels_ind <- apply(probabilities,1,function(x) inter_names[which.max(x)])

med_names = c('CGE1_3', 'CGE4_12', 'Lacunosum', 'Neurogliaform', 'Basket_1_Chandelier', 'OLM_Basket_1_GABA')
pdf(file.path(resultsdir, 'med_plots.pdf'))
for(i in 1:length(cur_types_sub)) {
  cur_types_inst = cur_types_sub[[i]]
  plot_val_sub <- rowSums(probabilities[,cur_types_inst])
  plot <- plot_puck_continuous(iv$puck, inter_barcodes, plot_val_sub, title = med_names[i], ylimit = c(0,1),my_pal = pals::brewer.blues(20))
  print(plot)
}
dev.off()
plot_puck_continuous(iv$puck, inter_barcodes, plot_val_sub, title = 10, ylimit = c(0,1),my_pal = pals::brewer.blues(20))
conf_thresh = 0.5
pdf(file.path(resultsdir,'cut_tree_levels.pdf'))
conf_proportion <- numeric(26)
names(conf_proportion) <- as.character(2:27)
correct <- numeric(26)
names(correct) <- as.character(2:27)
for(K_cat in 2:27) {
  agg <- t(aggregate(t(probabilities),by = list(cutree(hc1,k=K_cat)),FUN = sum))[2:(length(inter_barcodes)+1),]
  conf_inter <- apply(agg,1,max) > conf_thresh
  these_labels <- factor(apply(agg,1,which.max))
  #print(plot_class(puck, names(conf_inter[conf_inter]), these_labels[conf_inter], counter_barcodes = counter_barcodes, title = paste("N_levels=",K_cat)))
  conf_proportion[as.character(K_cat)] = sum(conf_inter) / length(conf_inter)
  total = 0
  agree = 0
  for(i in levels(new_labels)) {
    these_barcodes = which(new_labels == i & conf_inter)
    if(length(these_barcodes) >= 2) {
      total = total + length(these_barcodes)*(length(these_barcodes)-1)/2
      agree = agree + sum(unlist(lapply(table(these_labels[these_barcodes]),function(x) x*(x-1)/2)))
    }
  }
  correct[as.character(K_cat)] = agree / total
}
dev.off()

plot(2:27,correct,ylim = c(0,1))
plot(2:27,conf_proportion,ylim = c(0,1))

final_labels = factor(character(length(inter_barcodes)),levels = inter_names)
names(final_labels) = inter_barcodes
n_classes = length(levels(new_labels))
cluster_probabilities = probabilities
for(i in(1:n_classes)) {
  #final_labels[new_labels==i] <- names(which.max(colMeans(probabilities[inter_barcodes[new_labels==i],])))
  N_curr <- sum(new_labels==i)
  my_fact <- sqrt(N_curr)
  row_this <- colMeans(probabilities[inter_barcodes[new_labels==i],])^my_fact
  cluster_probabilities[inter_barcodes[new_labels==i], ] <- matrix(rep(row_this,each=N_curr), N_curr)
}
cluster_probabilities <- sweep(cluster_probabilities, 1, rowSums(cluster_probabilities) ,'/')
colSums(cluster_probabilities)
plot_class(puck, inter_barcodes, final_labels, counter_barcodes = counter_barcodes)
hist(apply(cluster_probabilities,1,max))

pdf(file.path(resultsdir,'cut_tree_levels_cluster.pdf'))
conf_proportion <- numeric(26)
names(conf_proportion) <- as.character(2:27)
for(K_cat in 2:27) {
  agg <- t(aggregate(t(cluster_probabilities),by = list(cutree(hc1,k=K_cat)),FUN = sum))[2:(length(inter_barcodes)+1),]
  conf_inter <- apply(agg,1,max) > conf_thresh
  these_labels <- factor(apply(agg,1,which.max))
  print(plot_class(puck, names(conf_inter[conf_inter]), these_labels[conf_inter], counter_barcodes = counter_barcodes, title = paste("N_levels=",K_cat)))
  conf_proportion[as.character(K_cat)] = sum(conf_inter) / length(conf_inter)
}
dev.off()
plot(2:27,conf_proportion,ylim=c(0,1))

plot_class(puck, names(conf_inter[conf_inter]), final_labels, counter_barcodes = counter_barcodes)
plot_class(puck, names(conf_inter[conf_inter & (final_labels %in% c('OLM1','OLM2','OLM3','OLM4'))]), final_labels, counter_barcodes = counter_barcodes)
plot_class(puck, names(conf_inter[conf_inter & (final_labels %in% c('CGE_1','CGE_2','CGE_3','CGE_4','CGE_5','CGE_6','CGE_7','CGE_8','CGE_9','CGE_10','CGE_11','CGE_12'))]), final_labels, counter_barcodes = counter_barcodes)
plot_class(puck, names(conf_inter[conf_inter & (final_labels %in% inter_names[19:22])]), final_labels, counter_barcodes = counter_barcodes)
plot_class(puck, names(conf_inter[conf_inter & (final_labels %in% inter_names[19:22])]), final_labels, counter_barcodes = counter_barcodes)
plot_class(puck, names(conf_inter[conf_inter & (final_labels %in% inter_names[c(3,18)])]), final_labels, counter_barcodes = counter_barcodes)
plot_class(puck, names(conf_inter[conf_inter & (final_labels %in% inter_names[c(1,2,16,17)])]), final_labels, counter_barcodes = counter_barcodes)
save.image(file=file.path(resultsdir,'inter.RData'))

marker_data = get_marker_data(cell_type_info[[2]], NULL, cell_type_info[[1]], gene_list, marker_provided = TRUE)
my_mat <- cell_type_info[[1]][gene_list[marker_data$cell_type == "OLM1"],c('OLM1','OLM2','OLM3','OLM4')]
my_mat$avg <- apply(my_mat[,2:4],1,mean)
my_mat$log_fc <- log(my_mat$OLM1) - log(my_mat$avg)
saveRDS(probabilities,file.path(resultsdir,'inter_probabilities.RDS'))
saveRDS(cluster_probabilities,file.path(resultsdir,'inter_cluster_probabilities.RDS'))
