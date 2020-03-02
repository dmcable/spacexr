library(RCTD)
library(Matrix)
#puck_file = paste0("SplitPuck/puck",fold_index,".RDS"),
iv <- init_RCTD(MIN_OBS=0, load_info_renorm = T) #initial variables
resultsdir <- paste0(iv$slideseqdir,"/results")
cell_type_means = iv$cell_type_info[[1]]
puck = iv$puck
#get initial classification
MIN_UMI = 300
<<<<<<< HEAD
N_fit = min(500,sum(puck@nUMI > MIN_UMI)) #500
fit_ind = sample(names(puck@nUMI[puck@nUMI > MIN_UMI]), N_fit)
beads = t(as.matrix(puck@counts[iv$gene_list,fit_ind]))
N_X = 50000 # 50000
delta = 1e-5; X_vals = (1:N_X)^1.5*delta
K_val = 100 #100
=======
N_fit = min(iv$config$N_fit,sum(puck@nUMI > MIN_UMI))
fit_ind = sample(names(puck@nUMI[puck@nUMI > MIN_UMI]), N_fit)
beads = t(as.matrix(puck@counts[iv$gene_list,fit_ind]))
N_X = iv$config$N_X # 50000
delta = 1e-5; X_vals = (1:N_X)^1.5*delta
K_val = iv$config$K_val #100
>>>>>>> dev
sigma = 1; use_Q = T
true_predictions = F
if(true_predictions) {
  metadir <- file.path(iv$slideseqdir,"MetaData")
  meta_data <- readRDS(file.path(metadir,"meta_data.RDS"))
  meta_df <- meta_data$meta_df; UMI_tot = meta_data$UMI_tot
  prediction = matrix(0, nrow = length(iv$gene_list), ncol = N_fit)
  colnames(prediction) = fit_ind; rownames(prediction) = iv$gene_list
  for (barcode in fit_ind) {
    type1 = meta_df[barcode, "first_type"]
    type2 = meta_df[barcode, "second_type"]
    p = meta_df[barcode, "first_UMI"] / UMI_tot
    prediction[, barcode] <- get_prediction_sparse(cell_type_means, iv$gene_list, UMI_tot, p, type1, type2)
  }
} else {
  #Q_mat <- readRDS(file.path(resultsdir,'QQ_mat.RDS'))
<<<<<<< HEAD
  results <- calc_Q_par(K_val, X_vals, sigma, big_params = F)
  print(paste('chooseSigma: calculating initial Q_mat with sigma = ',sigma))
=======
  print(paste('chooseSigma: calculating initial Q_mat with sigma = ',sigma))
  results <- calc_Q_par(K_val, X_vals, sigma, big_params = F)
>>>>>>> dev
  Q_mat <-t(as.data.frame(matrix(unlist(results), nrow=length(unlist(results[1])))))
  print(paste('chooseSigma: getting initial weights for #samples: ',N_fit))
  results = decompose_batch(puck@nUMI[fit_ind], iv$cell_type_info[[1]], beads, iv$gene_list, constrain = F)
  weights = Matrix(0, nrow = N_fit, ncol = iv$cell_type_info[[3]])
  rownames(weights) = fit_ind;
  colnames(weights) = iv$cell_type_info[[2]];
  for(i in 1:N_fit)
    weights[i,] = results[[i]]$weights
  prediction <- sweep(as.matrix(cell_type_means[iv$gene_list,]) %*% t(as.matrix(weights)), 2, puck@nUMI[fit_ind], '*')
}
counts <- t(beads)

<<<<<<< HEAD
sigma <- chooseSigma(prediction, counts, resultsdir, sigma_init = sigma)
=======
sigma <- chooseSigma(prediction, counts, resultsdir, sigma_init = sigma, N_epoch = iv$config$N_epoch)
>>>>>>> dev

results <- calc_Q_par(K_val, X_vals, sigma, big_params = T)
Q_mat <-t(as.data.frame(matrix(unlist(results), nrow=length(unlist(results[1])))))
saveRDS(Q_mat, file.path(resultsdir,'Q_mat.RDS'))
#big_params = F; sigma = 1
#if(F) {
#  for(i in 1:20) {
#    Q = numeric(length(Y))
#    for(k in as.numeric(names(table(Y)))) {
#      Q[Y == k] = calc_Q_mat_one(sigma, X[Y==k], k, batch = 100, big_params = big_params)
#    }

#    Q_d = numeric(length(Y))
#    for(k in as.numeric(names(table(Y)))) {
#      Q_d[Y == k] = calc_Q_d_one(sigma, X[Y==k], k, batch = 100, big_params = big_params)
#    }
#   print(sum(log(Q)))
#    sigma = sigma + sum(Q_d/Q)*alpha
#    print(sigma)
#  }
#}

#x = 1
#k = 3
#(get_Q(x, k, sigma, big_params = T) - get_Q(x, k, sigma, big_params = F))/get_Q(x, k, sigma, big_params = T)
#(get_Q_d(x, k, sigma, big_params = T) - get_Q_d(x, k, sigma, big_params = F))/get_Q_d(x, k, sigma, big_params = T)
