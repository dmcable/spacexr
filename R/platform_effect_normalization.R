#' Performs Platform Effect Normalization:
#'
#' Estimates bulk cell type composition and uses this
#' to estimate platform effects and normalize cell type proportions
#'
#' @section Additional Output Files:
#' Saves proportions in results/Bulk/weights.RDS.
#' Saves cell_type_info_renorm in MetaData/cell_type_info_renorm.RDS
#' Additionally, creates trace files of sigma fitting in results/sigmaBulk
#'
#' @param iv Initial Variables: meta data obtained from the \code{\link{init_RCTD}} function
#' @return Returns a list of cell_type_info_renorm (the normalized cell type profiles) and proportions
#' (the predicted bulk cell type proportions)
#' @export
fitBulk <- function(iv) {
  split_puck(iv$puck, iv$SpatialRNAdir, iv$n_puck_folds)
  bulkData <- prepareBulkData(iv$bulkdir, iv$cell_type_info[[1]], iv$puck, iv$gene_list)
  resultsdir <- paste0(iv$SpatialRNAdir,"/results")
  use_Q <<- F
  ABS_MAX = 30000 #prevent numerical issues
  X_vals <<- c(min(max(bulkData$X+1),ABS_MAX))
  K_val <<- min(ABS_MAX,max(bulkData$b + 1))
  sigmavar <<- 1
  print('fitBulk: decomposing bulk')
  decompose_results <- decompose_full(iv$cell_type_info[[1]], iv$gene_list, sum(iv$puck@nUMI), bulkData$b, verbose = T, constrain = F, MIN_CHANGE = iv$config$MIN_CHANGE_BULK*10, n.iter = 100)
  prediction <- as.matrix(bulkData$X) %*% decompose_results$weights
  sigmavar <<- sqrt(mean((log(prediction) - log(bulkData$b))^2))
  sigmavar <<- chooseSigma(prediction, bulkData$b, resultsdir, sigma_init = sigmavar, N_epoch = iv$config$N_epoch_bulk, folder_id = "Bulk1")
  decompose_results <- decompose_full(iv$cell_type_info[[1]], iv$gene_list, sum(iv$puck@nUMI), bulkData$b, verbose = T, constrain = F, MIN_CHANGE = iv$config$MIN_CHANGE_BULK, n.iter = 100)
  prediction <- as.matrix(bulkData$X) %*% decompose_results$weights
  sigma_prev <- sigmavar
  sigmavar <<- chooseSigma(prediction, bulkData$b, resultsdir, sigma_init = sigmavar, N_epoch = iv$config$N_epoch_bulk, folder_id = "Bulk2")
  if(abs(sigmavar - sigma_prev)/(max(0.1,sigmavar)) > 0.02)
    decompose_results <- decompose_full(iv$cell_type_info[[1]], iv$gene_list, sum(iv$puck@nUMI), bulkData$b, verbose = T, constrain = F, MIN_CHANGE = iv$config$MIN_CHANGE_BULK, n.iter = 100)
  proportions = decompose_results$weights
  saveRDS(proportions, paste0(resultsdir,'/Bulk/weights.RDS'))

  gene_list = get_de_genes(iv$cell_type_info, iv$puck, fc_thresh = iv$config$fc_cutoff_reg, expr_thresh = iv$config$gene_cutoff_reg, MIN_OBS = 0)
  metadir <- file.path(iv$SpatialRNAdir,"MetaData")
  if(!dir.exists(metadir))
    dir.create(metadir)
  #generate the celltypeinfo-renorm
  cell_type_info_renorm = iv$cell_type_info
  cell_type_info_renorm[[1]] = get_norm_ref(iv$puck, iv$cell_type_info[[1]], gene_list, proportions)
  saveRDS(cell_type_info_renorm, file.path(metadir, "cell_type_info_renorm.RDS"))
  return(list(cell_type_info_renorm = cell_type_info_renorm, proportions = proportions))
}


chooseSigma <- function(prediction, counts, resultsdir, sigma_init = 1, N_epoch = 15, folder_id = "") {
  X = as.vector(prediction)
  X = pmax(X, 1e-4)
  Y = as.vector(counts)
  num_sample = min(500000, length(X)) #300000
  big_params = F
  use_ind = sample(1:length(X), num_sample)
  X = X[use_ind]
  Y = Y[use_ind]
  sigma = sigma_init
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
  sigresdir = file.path(resultsdir, paste0("sigma",folder_id))
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
  return(sigma)
}

#' Estimates sigma_c by maximum likelihood
#'
#' Estimates sigma_c by stochastic gradient descent
#'
#' @param iv Initial Variables: meta data obtained from the \code{\link{init_RCTD}} function
#' @param resultsdir: directory of results folder
#' @return returns sigma, a numeric representing estimated sigma_c
#' (the predicted bulk cell type proportions)
#' @export
choose_sigma_c <- function(iv,resultsdir) {
  cell_type_means = iv$cell_type_info[[1]]
  puck = iv$puck
  #get initial classification
  MIN_UMI = 300
  N_fit = min(iv$config$N_fit,sum(puck@nUMI > MIN_UMI))
  fit_ind = sample(names(puck@nUMI[puck@nUMI > MIN_UMI]), N_fit)
  beads = t(as.matrix(puck@counts[iv$gene_list,fit_ind]))
  sigma = 1;
  delta <<- 1e-5
  print(paste('chooseSigma: calculating initial Q_mat with sigma = ',sigma))
  results <- calc_Q_par(iv$config$K_val, (1:iv$config$N_X)^1.5*delta, sigma, big_params = F)
  Q_mat_loc <-t(as.data.frame(matrix(unlist(results), nrow=length(unlist(results[1])))))
  set_likelihood_vars(Q_mat_loc)
  print(paste('chooseSigma: getting initial weights for #samples: ',N_fit))
  results = decompose_batch(puck@nUMI[fit_ind], iv$cell_type_info[[1]], beads, iv$gene_list, constrain = F)
  weights = Matrix(0, nrow = N_fit, ncol = iv$cell_type_info[[3]])
  rownames(weights) = fit_ind;
  colnames(weights) = iv$cell_type_info[[2]];
  for(i in 1:N_fit)
    weights[i,] = results[[i]]$weights
  prediction <- sweep(as.matrix(cell_type_means[iv$gene_list,]) %*% t(as.matrix(weights)), 2, puck@nUMI[fit_ind], '*')
  counts <- t(beads)
  sigma <- chooseSigma(prediction, counts, resultsdir, sigma_init = sigma, N_epoch = iv$config$N_epoch)
  return(sigma)
}
