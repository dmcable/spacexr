library(RCTD)
library(Matrix)
#puck_file = paste0("SplitPuck/puck",fold_index,".RDS"),
iv <- init_RCTD(MIN_OBS=0, load_info_renorm = T) #initial variables
resultsdir <- paste0(iv$SpatialRNAdir,"/results")
sigma <- choose_sigma_c(iv,resultsdir)
results <- calc_Q_par(K_val, X_vals, sigma, big_params = T)
Q_mat <-t(as.data.frame(matrix(unlist(results), nrow=length(unlist(results[1])))))
saveRDS(Q_mat, file.path(resultsdir,'Q_mat.RDS'))
