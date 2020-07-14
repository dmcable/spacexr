#' Performs Platform Effect Normalization:
#'
#' Estimates bulk cell type composition and uses this
#' to estimate platform effects and normalize cell type proportions
#'
#' @param RCTD an \code{\linkS4class{RCTD}} object after running the \code{\link{create.RCTD}} function.
#' @return Returns an \code{\linkS4class{RCTD}} object normalized for platform effects.
#' @export
fitBulk <- function(RCTD) {
  bulkData <- prepareBulkData(RCTD@cell_type_info$info[[1]], RCTD@spatialRNA, RCTD@internal_vars$gene_list_bulk)
  print('fitBulk: decomposing bulk')
  decompose_results <- decompose_full(bulkData$X, sum(RCTD@spatialRNA@nUMI),
                                      bulkData$b, verbose = F, constrain = F, MIN_CHANGE = RCTD@config$MIN_CHANGE_BULK,
                                      n.iter = 100, bulk_mode = T)
  RCTD@internal_vars$proportions <- decompose_results$weights
  RCTD@cell_type_info$renorm = RCTD@cell_type_info$info
  RCTD@cell_type_info$renorm[[1]] = get_norm_ref(RCTD@spatialRNA, RCTD@cell_type_info$info[[1]], RCTD@internal_vars$gene_list_bulk, decompose_results$weights)
  return(RCTD)
}

chooseSigma2 <- function(prediction, counts, Q_mat_all, X_vals, sigma) {
  X = as.vector(prediction); X = pmax(X, 1e-4)
  Y = as.vector(counts); num_sample = min(500000, length(X)) #300000
  use_ind = sample(1:length(X), num_sample)
  X = X[use_ind]; Y = Y[use_ind]
  mult_fac_vec <- (7:12 + 0.5)/10; sigma_ind <- c(10:70, (36:100)*2)
  score_vec <- numeric(length(sigma_ind))
  for(i in 1:length(sigma_ind)) {
    sigma <- sigma_ind[i]
    set_likelihood_vars(Q_mat_all[[as.character(sigma)]],X_vals)
    best_val <- calc_log_l_vec(X*mult_fac_vec[1], Y)
    for(mult_fac in mult_fac_vec[2:length(mult_fac_vec)])
      best_val <- min(best_val, calc_log_l_vec(X*mult_fac, Y))
    score_vec[i] <- best_val
  }
  sigma = sigma_ind[which.min(score_vec)]
  print(min(score_vec))
  return(sigma)
}


chooseSigma <- function(prediction, counts, Q_mat_all, X_vals, sigma) {
  X = as.vector(prediction); X = pmax(X, 1e-4)
  Y = as.vector(counts); num_sample = min(1000000, length(X)) #300000
  use_ind = sample(1:length(X), num_sample)
  X = X[use_ind]; Y = Y[use_ind]
  mult_fac_vec <- (8:12)/10; sigma_ind <- c(10:70, (36:100)*2)
  si <- which(sigma_ind == round(sigma))
  sigma_ind <- sigma_ind[max(1,si - 8):min(si+8,length(sigma_ind))]
  score_vec <- numeric(length(sigma_ind))
  for(i in 1:length(sigma_ind)) {
    sigma <- sigma_ind[i]
    set_likelihood_vars(Q_mat_all[[as.character(sigma)]],X_vals)
    best_val <- calc_log_l_vec(X*mult_fac_vec[1], Y)
    for(mult_fac in mult_fac_vec[2:length(mult_fac_vec)])
      best_val <- min(best_val, calc_log_l_vec(X*mult_fac, Y))
    score_vec[i] <- best_val
    #score_vec[i] <- calc_log_l_vec(X, Y)
  }
  sigma = sigma_ind[which.min(score_vec)]
  return(sigma)
}

#' Estimates sigma_c by maximum likelihood
#'
#' @param RCTD an \code{\linkS4class{RCTD}} object after running the \code{\link{fitBulk}} function.
#' @return Returns an \code{\linkS4class{RCTD}} with the estimated \code{sigma_c}.
#' @export
choose_sigma_c <- function(RCTD) {
  puck = RCTD@spatialRNA; MIN_UMI = 300; sigma = 100
  #Q_mat_all <- readRDS('/Users/dcable/Documents/MIT/Research/Rafalab/Projects/slideseq/Cell Demixing/ContentStructure/RCTD/Qmat/Q_mat_c.rds')
  Q1 <- readRDS(system.file("extdata", "Qmat/Q_mat_1.rds", package = "RCTD"))
  Q2 <- readRDS(system.file("extdata", "Qmat/Q_mat_2.rds", package = "RCTD"))
  Q_mat_all <- c(Q1,Q2)
  sigma_vals <- names(Q_mat_all)
  #X_vals <- readRDS('/Users/dcable/Documents/MIT/Research/Rafalab/Projects/slideseq/Cell Demixing/ContentStructure/RCTD/Qmat/X_vals.rds')
  X_vals <- readRDS(system.file("extdata", "Qmat/X_vals.rds", package = "RCTD"))
  #get initial classification
  N_fit = min(RCTD@config$N_fit,sum(puck@nUMI > MIN_UMI))
  fit_ind = sample(names(puck@nUMI[puck@nUMI > MIN_UMI]), N_fit)
  beads = t(as.matrix(puck@counts[RCTD@internal_vars$gene_list_reg,fit_ind]))
  print(paste('chooseSigma: using initial Q_mat with sigma = ',sigma/100))
  for(iter in 1:myRCTD@config$N_epoch) {
    set_likelihood_vars(Q_mat_all[[as.character(sigma)]], X_vals)
    #print(paste('chooseSigma: getting initial weights for #samples: ',N_fit))
    results = decompose_batch(puck@nUMI[fit_ind], RCTD@cell_type_info$renorm[[1]], beads, RCTD@internal_vars$gene_list_reg, constrain = F, max_cores = RCTD@config$max_cores)
    weights = Matrix(0, nrow = N_fit, ncol = RCTD@cell_type_info$renorm[[3]])
    rownames(weights) = fit_ind; colnames(weights) = RCTD@cell_type_info$renorm[[2]];
    for(i in 1:N_fit)
      weights[i,] = results[[i]]$weights
    prediction <- sweep(as.matrix(RCTD@cell_type_info$renorm[[1]][RCTD@internal_vars$gene_list_reg,]) %*% t(as.matrix(weights)), 2, puck@nUMI[fit_ind], '*')
    print(paste('Likelihood value:',calc_log_l_vec(as.vector(prediction), as.vector(t(beads)))))
    sigma_prev <- sigma
    sigma <- chooseSigma(prediction, t(beads), Q_mat_all, X_vals, sigma)
    print(paste('Sigma value: ', sigma/100))
    if(sigma == sigma_prev)
      break
  }
  RCTD@internal_vars$sigma <- sigma/100
  RCTD@internal_vars$Q_mat <- Q_mat_all[[as.character(sigma)]]
  RCTD@internal_vars$X_vals <- X_vals
  return(RCTD)
}
