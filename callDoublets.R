library(RCTD)
library(Matrix)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("callDoublets: At least one argument must be supplied: fold_index", call.=FALSE)
}
fold_index = as.integer(args[1])
print(paste0("callDoublets: loading fold index: ",fold_index))
iv <- init_RCTD(puck_file = paste0("SplitPuck/puck",fold_index,".RDS"), MIN_OBS=0, load_info_renorm = T) #initial variables
resultsdir <- paste0(iv$slideseqdir,"/results")
Q_mat <- readRDS(file.path(resultsdir,'QQ_mat.RDS'))
N_X = dim(Q_mat)[2]; delta = 1e-5; X_vals = (1:N_X)^1.5*delta
K_val = dim(Q_mat)[1] - 3
class_df <- get_class_df(iv$cell_type_info[[2]])
results = process_beads_batch(iv$cell_type_info, iv$gene_list, iv$puck, class_df = class_df, constrain = F, doublet_mode = T)
saveRDS(results, paste0(iv$slideseqdir,"/SplitPuckResults/results",fold_index,".RDS"))
