#testing performance of ordinary least squares

library(RCTD)
library(Matrix)
source('Plotting/figure_utils.R')

DropViz <- T
iv <- init_RCTD(gene_list_reg = F, get_proportions = DropViz, load_info = F)
if(DropViz) {
  proportions = iv$proportions
  cell_type_info_unnorm <- iv$cell_type_info
}
puck = iv$puck
iv <- init_RCTD(load_info_renorm = T) #initial variables
if(DropViz) {
  common_cell_types = c("Astrocytes", "Bergmann", "Endothelial", "Fibroblast", "Golgi", "Granule", "MLI1", "MLI2", "Oligodendrocytes", "Polydendrocytes", "Purkinje", "UBCs")
} else {
  common_cell_types <- iv$cell_type_info[[2]]
}
resultsdir <- paste0(iv$slideseqdir,"/results")
metadir <- file.path(iv$slideseqdir,"MetaData")
meta_data <- readRDS(file.path(metadir,"meta_data.RDS"))
meta_df <- meta_data$meta_df
Q_mat <- readRDS(file.path(resultsdir,'Q_mat.RDS'))
N_X = dim(Q_mat)[2]; delta = 1e-5; X_vals = (1:N_X)^1.5*delta
K_val = dim(Q_mat)[1] - 3; use_Q = T

my_barc <- c(rownames(meta_df[meta_df$first_UMI == 0,]),rownames(meta_df[meta_df$first_UMI == 1000,])) #c(0,1000)
true_names <- c(as.character(meta_df[rownames(meta_df[meta_df$first_UMI == 0,]), "second_type"]),as.character(meta_df[rownames(meta_df[meta_df$first_UMI == 1000,]), "first_type"]))
names(true_names) <- my_barc
puck <- restrict_puck(puck, my_barc)
puck@cell_labels <- factor(true_names, levels = iv$cell_type_info[[2]])
#puck_sm <- restrict_puck(puck, rownames(meta_df[meta_df$first_UMI == 0,]))
#Figure 2A: OLS Prediction works on Training data, but not cross-reference
test_results = process_data(puck, iv$gene_list, cell_type_info_unnorm, proportions = NULL, trust_model = F, constrain = F, OLS = T)
#scratch
cell_type_lev = factor(1:iv$cell_type_info[[3]])
cell_type_map = data.frame(cindex = 1:iv$cell_type_info[[3]], row.names = iv$cell_type_info[[2]])
pred_labels = test_results[[3]]
true_labels = lapply(puck@cell_labels, function(x) cell_type_map[as.character(x),"cindex"])
true_labels = as.integer(puck@cell_labels[as.character(colnames(puck@counts))])
conf_mat = caret::confusionMatrix(factor(pred_labels,cell_type_lev),factor(true_labels,cell_type_lev))
norm_conf = sweep(conf_mat$table, 2, colSums(conf_mat$table), '/')
rownames(norm_conf) <- iv$cell_type_info[[2]]; colnames(norm_conf) <- iv$cell_type_info[[2]]
#end
library(reshape2)
data <- melt(norm_conf[,common_cell_types])
saveRDS(data,file="plotting/Results/cross_confusion.RDS")
