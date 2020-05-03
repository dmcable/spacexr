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
gene_list_inter <- iv$gene_list


#assign slide-seq interneurons to probabilities of subclusters

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
hc1 <- readRDS(file.path(refdir,"Subcluster/dendrogram.RDS"))


#spatially cluster interneurons
d <- dist(puck@coords[inter_barcodes,], method = "euclidean")
hc1 <- hclust(d, method = "average")
plot(hc1, cex = 0.6, hang = -1)
num_clusters = 200
cutree(hc1,k=num_clusters)

barcodes_cur = inter_barcodes
my_class <- as.factor(cutree(hc1,k=200))
my_class <- fct_shuffle(my_class)
plot_class(puck, barcodes_cur, my_class)

#manually split doublet spatial clusters
relabel <- read.csv(file.path(iv$slideseqdir,'cluster_relabels.csv'))
for(i in unique(relabel$Cluster))
  relabel[relabel$Cluster==i,"barcodes"] <- mapvalues(relabel[relabel$Cluster==i,]$Index,from = which(rownames(puck@coords) %in% inter_barcodes[my_class==i]), to=inter_barcodes[my_class==i])
rownames(relabel) <- relabel$barcodes
new_labels <- as.character(my_class)
names(new_labels) <- names(my_class)
new_labels[relabel$barcode] <- apply(relabel,1,function(x) paste(x[3],x[2],sep='_'))
new_labels <- as.factor(new_labels)
new_labels <- mapvalues(new_labels, from = levels(new_labels), to = sample(1:length(levels(new_labels))))
final_labels = factor(character(length(inter_barcodes)),levels = inter_names)
names(final_labels) = inter_barcodes
n_classes = length(levels(new_labels))
for(i in(1:n_classes)) {
  final_labels[new_labels==i] <- names(which.max(colMeans(probabilities[inter_barcodes[new_labels==i],])))
}
counter_barcodes <- barcodes[results_df$spot_class == "singlet" & results_df$first_type %in% c("CA1","CA3","Denate")]

plot_class(puck, inter_barcodes, final_labels, counter_barcodes = counter_barcodes)

#new analysis by Dylan to add the likelihoods
#agregate probabilities across each cluster
likelihoods <- inter_df[,4:30]
final_labels = factor(character(length(inter_barcodes)),levels = inter_names)
names(final_labels) = inter_barcodes
n_classes = length(levels(new_labels))
for(i in(1:n_classes)) {
  row_this <- colSums(likelihoods[inter_barcodes[new_labels==i],])
  final_labels[inter_barcodes[new_labels==i]] <- names(which.min(row_this))
  conf_inter[inter_barcodes[new_labels==i]] <- -(min(row_this) - row_this[order(row_this)[2]]) >= 10
}
plot_class(puck, inter_barcodes[conf_inter], final_labels, counter_barcodes = counter_barcodes)
save(puck, inter_barcodes, conf_inter,final_labels, counter_barcodes, file = "AnalysisPaper/Plotting/Results/inter27.RData")
hist(apply(cluster_probabilities,1,max))

#next do 3 main subclasses
inter_names_coarse<- c('CGE', "Basket_OLM" , "Neurogliaform_Lacunosum")
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

