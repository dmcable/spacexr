#used for decomposition experiment on known cell types
library(RCTD)
library(Matrix)
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("callDoublets: At least one argument must be supplied: fold_index", call.=FALSE)
}
fold_index = as.integer(args[1])
print(paste0("callDoublets: loading fold index: ",fold_index))
iv <- init_RCTD(puck_file = paste0("SplitPuck/puck",fold_index,".RDS"), MIN_OBS=0, load_info_renorm = T) #initial variables
DropViz <- F
if(DropViz) {
  common_cell_types = c("Astrocytes", "Bergmann", "Endothelial", "Fibroblast", "Golgi", "Granule", "MLI1", "MLI2", "Oligodendrocytes", "Polydendrocytes", "Purkinje", "UBCs")
} else
  common_cell_types <- iv$cell_type_info[[2]]
resultsdir <- paste0(iv$slideseqdir,"/results")
Q_mat <- readRDS(file.path(resultsdir,'Q_mat.RDS'))
N_X = dim(Q_mat)[2]; delta = 1e-5; X_vals = (1:N_X)^1.5*delta
K_val = dim(Q_mat)[1] - 3; use_Q = T
meta_data <- readRDS(file.path(iv$slideseqdir, 'MetaData/meta_data.RDS'))
meta_df <- meta_data$meta_df
n_folds = iv$n_puck_folds; N = meta_data$N_samples
start_index <- (round((fold_index-1)*N/n_folds) + 1)
end_index <- round(fold_index*N/n_folds)
meta_df <- meta_df[start_index:end_index,]
meta_df$first_type <- factor(meta_df$first_type, common_cell_types)
meta_df$second_type <- factor(meta_df$second_type, common_cell_types)

results = process_beads_sparse(iv$cell_type_info, iv$gene_list, iv$puck, meta_df, constrain = F)
saveRDS(results, paste0(iv$slideseqdir,"/DecomposeResults/results",fold_index,".RDS"))
