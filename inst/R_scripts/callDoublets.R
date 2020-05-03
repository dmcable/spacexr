### Running RCTD

# After computing the quadrature, we can run RCTD in parallel on the `SpatialRNA` dataset.
# Because RCTD runs over batches, or folds, of the data, we must provide a `fold_index`
# to run on. This is inputed as a command line parameter.
# RCTD requires the function `set_likelihood_vars` to be run, which sets the quadrature
# as a global variable. RCTD runs with the `process_beads_batch` function, which uses
# parallel processing to estimate cell type proportions on each pixel in the `SpatialRNA` dataset.
# The results are returned and saved for future interpretation.

library(RCTD)
library(Matrix)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("callDoublets: At least one argument must be supplied: fold_index", call.=FALSE)
}
fold_index = as.integer(args[1])
print(paste0("callDoublets: loading fold index: ",fold_index))
iv <- init_RCTD(puck_file = paste0("SplitPuck/puck",fold_index,".RDS"), MIN_OBS=0, load_info_renorm = T) #initial variables
resultsdir <- paste0(iv$SpatialRNAdir,"/results")
Q_mat_loc <- readRDS(file.path(resultsdir,'Q_mat.RDS'))
set_likelihood_vars(Q_mat_loc)
class_df <- data.frame(iv$cell_type_info[[2]], row.names = iv$cell_type_info[[2]]); colnames(class_df)[1] = "class"
results = process_beads_batch(iv$cell_type_info, iv$gene_list, iv$puck, class_df = class_df, constrain = F)
saveRDS(results, paste0(iv$SpatialRNAdir,"/SplitPuckResults/results",fold_index,".RDS"))

