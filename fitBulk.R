library(RCTD)
library(Matrix)
#scratch
iv <- init_RCTD(gene_list_reg = F, get_proportions = F) #initial variables
results <- prepareBulkData(iv$bulkdir, iv$cell_type_info[[1]], iv$puck, iv$gene_list)
resultsdir <- paste0(iv$slideseqdir,"/results")
use_Q = F
X_vals = c(max(results$X+1))
K_val = max(results$b + 1)
sigma = 1
decompose_results <- decompose_full(iv$cell_type_info[[1]], iv$gene_list, sum(iv$puck@nUMI), results$b, verbose = T, constrain = F, MIN_CHANGE = 0.0001, n.iter = 100)
prediction <- as.matrix(results$X) %*% decompose_results$weights
sigma <- sqrt(mean((log(prediction) - log(results$b))^2))
sigma <- chooseSigma(prediction, results$b, resultsdir, sigma_init = sigma, N_epoch = 30, folder_id = "Bulk1")
decompose_results <- decompose_full(iv$cell_type_info[[1]], iv$gene_list, sum(iv$puck@nUMI), results$b, verbose = T, constrain = F, MIN_CHANGE = 0.0001, n.iter = 100)
prediction <- as.matrix(results$X) %*% decompose_results$weights
sigma <- chooseSigma(prediction, results$b, resultsdir, sigma_init = sigma, N_epoch = 30, folder_id = "Bulk2")
decompose_results <- decompose_full(iv$cell_type_info[[1]], iv$gene_list, sum(iv$puck@nUMI), results$b, verbose = T, constrain = F, MIN_CHANGE = 0.0001, n.iter = 100)
proportions = decompose_results$weights
saveRDS(proportions, paste0(resultsdir,'/Bulk/weights.RDS'))
#split_puck(puck, slideseqdir, config$n_puck_folds)

gene_list = get_de_genes(iv$cell_type_info, iv$puck, fc_thresh = iv$config$fc_cutoff_reg, expr_thresh = iv$config$gene_cutoff_reg, MIN_OBS = 0)
metadir <- file.path(iv$slideseqdir,"MetaData")
if(!dir.exists(metadir))
  dir.create(metadir)
#generate the celltypeinfo-renorm
cell_type_info_renorm = iv$cell_type_info
cell_type_info_renorm[[1]] = get_norm_ref(iv$puck, iv$cell_type_info[[1]], gene_list, proportions)
saveRDS(cell_type_info_renorm, file.path(metadir, "cell_type_info_renorm.RDS"))
