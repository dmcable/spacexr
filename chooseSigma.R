library(RCTD)
library(Matrix)
fold_index = 1
iv <- init_RCTD(puck_file = paste0("SplitPuck/puck",fold_index,".RDS"), MIN_OBS=0, load_info_renorm = T) #initial variables
resultsdir <- paste0(iv$slideseqdir,"/results")
cell_type_means = iv$cell_type_info[[1]]
puck = iv$puck
#get initial classification
N_fit = min(500,dim(puck@counts)[2]) #500
fit_ind = sample(colnames(puck@counts), N_fit)
beads = t(as.matrix(puck@counts[iv$gene_list,fit_ind]))
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
  Q_mat <- readRDS(file.path(resultsdir,'QQ_mat.RDS'))
  N_X = dim(Q_mat)[2]; delta = 1e-5; X_vals = (1:N_X)^1.5*delta
  K_val = dim(Q_mat)[1] - 3
  results = decompose_batch(puck@nUMI[fit_ind], iv$cell_type_info[[1]], beads, iv$gene_list, constrain = F)
  weights = Matrix(0, nrow = N_fit, ncol = iv$cell_type_info[[3]])
  rownames(weights) = fit_ind;
  colnames(weights) = iv$cell_type_info[[2]];
  for(i in 1:N_fit)
    weights[i,] = results[[i]]$weights
  prediction <- sweep(as.matrix(cell_type_means) %*% t(as.matrix(weights)), 2, puck@nUMI[fit_ind], '*')
}
counts <- t(beads)
X = as.vector(prediction)
X = pmax(X, 1e-4)
Y = as.vector(counts)
num_sample = min(500000, length(X)) #300000
big_params = F
use_ind = sample(1:length(X), num_sample)
X = X[use_ind]
Y = Y[use_ind]
sigma = 1
alpha_init = 0.0001; batch = 25
X_batch = list()
Y_batch = list()
for(k in as.numeric(names(table(Y)))) {
  X_vals <- X[Y==k]; N_X = length(X_vals)
  for(b in 1:ceiling(N_X/batch)) {
    X_ind = (batch*(b-1) + 1):min((batch*b),N_X)
    curr_X = X_vals[X_ind]
    X_batch[[length(X_batch) + 1]] <- curr_X
    Y_batch[[length(Y_batch) + 1]] <- k
  }
}
ordering <- sample(1:length(X_batch))
sigma_vals = list(sigma)
loss_vals = list()
N_epoch = 15 #30
for(j in 1:N_epoch) {
  alpha = alpha_init / j
  total_loss = 0
  for(i in ordering) {
    Q = get_Q(X_batch[[i]], Y_batch[[i]], sigma, big_params = big_params)
    Q_d = get_Q_d(X_batch[[i]], Y_batch[[i]], sigma, big_params = big_params)
    sigma = sigma + sum(Q_d/Q)*alpha
    if(i%%100 == 0)
      print(sigma)
    total_loss = total_loss + sum(log(Q))
    sigma_vals[[length(sigma_vals) + 1]] <- sigma
  }
  print(total_loss)
  loss_vals[[length(loss_vals) + 1]] <- total_loss
}
sigresdir = file.path(resultsdir, "sigma")
if(!dir.exists(sigresdir))
  dir.create(sigresdir)

saveRDS(sigma_vals, file.path(sigresdir,"sigma_vals.RDS"))
saveRDS(loss_vals, file.path(sigresdir,"loss_vals.RDS"))
saveRDS(sigma, file.path(sigresdir,"sigma.RDS"))

pdf(file.path(sigresdir,"sigma_trace.pdf"))
plot(unlist(sigma_vals),type = 'n')
lines(unlist(sigma_vals))
dev.off()

pdf(file.path(sigresdir,"loss_trace.pdf"))
plot(unlist(loss_vals),type = 'n')
lines(unlist(loss_vals))
dev.off()

#big_params = F; sigma = 1
if(F) {
  for(i in 1:20) {
    Q = numeric(length(Y))
    for(k in as.numeric(names(table(Y)))) {
      Q[Y == k] = calc_Q_mat_one(sigma, X[Y==k], k, batch = 100, big_params = big_params)
    }

    Q_d = numeric(length(Y))
    for(k in as.numeric(names(table(Y)))) {
      Q_d[Y == k] = calc_Q_d_one(sigma, X[Y==k], k, batch = 100, big_params = big_params)
    }
    print(sum(log(Q)))
    sigma = sigma + sum(Q_d/Q)*alpha
    print(sigma)
  }
}

#x = 1
#k = 3
#(get_Q(x, k, sigma, big_params = T) - get_Q(x, k, sigma, big_params = F))/get_Q(x, k, sigma, big_params = T)
#(get_Q_d(x, k, sigma, big_params = T) - get_Q_d(x, k, sigma, big_params = F))/get_Q_d(x, k, sigma, big_params = T)
