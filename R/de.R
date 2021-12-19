#' Runs RCTDE on a \code{\linkS4class{RCTD}} object (for a single predictive variable)
#'
#' @param RCTD an \code{\linkS4class{RCTD}} object with annotated cell types e.g. from the \code{\link{run.RCTD}} function.
#' @param explanatory.variable the predictive variable used for detecting differential expression in RCTDE. Names of the variable
#' are the `SpatialRNA` pixel names, and values should be standardized between 0 and 1.
#' @param cell_types the cell types used for RCTDE. If null, cell types will be chosen with aggregate occurences of
#' at least `cell_type_threshold`, as aggregated by \code{\link{choose_cell_types}}
#' @param cell_type_threshold (default 125) min occurence of number of cells for each cell type to be used, as aggregated by \code{\link{choose_cell_types}}
#' @param gene_threshold (default 5e-5) minimum average expression required for selecting genes
#' @param doublet_mode (default TRUE) if TRUE, used RCTD doublet mode weights. Otherwise, uses RCTD full mode weights
#' @param test_mode (default 'individual') if 'individual', tests for DE individually for each parameter. If 'categorical', then tests for differences
#' across multiple categorical parameters
#' @param sigma_gene (default TRUE) if TRUE, fits gene specific overdispersion parameter
#' @param PRECISION.THRESHOLD (default 0.01) for checking for convergence, the maximum parameter change per algorithm step
#' @param cell_types_present cell types (a superset of `cell_types`) to be considered as occuring often enough
#' to consider for gene expression contamination
#' @param fdr (default 0.01) false discovery rate for hypothesis testing
#' @return an \code{\linkS4class{RCTD}} object containing the results of the RCTDE algorithm. Contains objects \code{de_results},
#' which contain the results of the RCTDE algorithm including `gene_fits`, which contains the results of fits on individual genes,
#' in addition `res_gene_list`, a list, for each cell type, of significant genes detected by RCTDE.
#' Additionally, the object contains `internal_vars_de`` a list of variables that are used internally by RCTDE
#' @export
run.RCTDE.single <- function(myRCTD, explanatory.variable,  cell_types = NULL, cell_type_threshold = 125,
                          gene_threshold = 5e-5, doublet_mode = T, test_mode = 'individual', thresh_val = NULL,
                          sigma_gene = T, PRECISION.THRESHOLD = 0.01, cell_types_present = NULL, fdr = .01) {
  X2 <- build.designmatrix.single(explanatory.variable, myRCTD)
  barcodes <- rownames(X2)
  return(run.RCTDE(myRCTD, X2, barcodes, cell_types, gene_threshold = gene_threshold, cell_type_threshold = cell_type_threshold,
                       doublet_mode = doublet_mode, test_mode = test_mode, params_to_test = 2,
                       thresh_val = thresh_val, sigma_gene = sigma_gene,
                       PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                       cell_types_present = cell_types_present, fdr = fdr))
}

# run DE nonparametrically
run.RCTDE.nonparam <- function(myRCTD, df = 15, barcodes = NULL, cell_types = NULL,
                            cell_type_threshold = 125, gene_threshold = 5e-5, doublet_mode = T,
                            test_mode = 'individual', thresh_val = NULL, sigma_gene = T,
                            PRECISION.THRESHOLD = 0.01, cell_types_present = NULL) {
  X2 <- build.designmatrix.nonparam(myRCTD, barcodes = barcodes, df = df)
  barcodes <- rownames(X2)
  return(run.RCTDE(myRCTD, X2, barcodes, cell_types, gene_threshold = gene_threshold,
                       doublet_mode = doublet_mode, test_mode = test_mode, cell_type_threshold = cell_type_threshold,
                       thresh_val = thresh_val, sigma_gene = sigma_gene,
                       PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                       cell_types_present = cell_types_present, params_to_test = 2:df))
}

run.RCTDE.regions <- function(myRCTD, region_list, cell_types = NULL,
                           cell_type_threshold = 125, gene_threshold = 5e-5, doublet_mode = T,
                           test_mode = 'categorical', thresh_val = NULL, sigma_gene = T,
                           PRECISION.THRESHOLD = 0.01, cell_types_present = NULL) {
  X2 <- build.designmatrix.regions(myRCTD, region_list)
  barcodes <- rownames(X2)
  return(run.RCTDE(myRCTD, X2, barcodes, cell_types, cell_type_threshold = cell_type_threshold, gene_threshold = gene_threshold,
                       doublet_mode = doublet_mode, test_mode = test_mode,
                       thresh_val = thresh_val, sigma_gene = sigma_gene, params_to_test = 1:dim(X2)[2],
                       PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                       cell_types_present = cell_types_present))
}

# cell_type_interactions: default to be vector of TRUE's indicating interactiongs for all parameters
run.RCTDE <- function(myRCTD, X, barcodes, cell_types, gene_threshold = 5e-5, cell_type_threshold = 125,
                          doublet_mode = T, test_mode = 'individual', thresh_val = NULL,
                          sigma_gene = T, PRECISION.THRESHOLD = 0.01, cell_types_present = NULL,
                          add_res_genes = T, fdr = .01, cell_type_specific = NULL, params_to_test = NULL) {
  X <- check_designmatrix(X, 'run.RCTDE', require_2d = TRUE)
  if(is.null(cell_type_specific))
    cell_type_specific <- !logical(dim(X)[2])
  check_cell_type_specific(cell_type_specific, dim(X)[2])
  X1 <- X[,!cell_type_specific]; X2 <- X[,cell_type_specific]
  return(run.RCTDE.general(myRCTD, X1, X2, barcodes, cell_types, cell_type_threshold = cell_type_threshold,
                           gene_threshold = gene_threshold,
                       doublet_mode = doublet_mode, test_mode = test_mode, thresh_val = thresh_val,
                       sigma_gene = sigma_gene, PRECISION.THRESHOLD = PRECISION.THRESHOLD, params_to_test = params_to_test,
                       cell_types_present = cell_types_present, add_res_genes = add_res_genes, fdr = fdr))
}

run.RCTDE.general <- function(myRCTD, X1, X2, barcodes, cell_types, gene_threshold = 5e-5, cell_type_threshold = 125,
                          doublet_mode = T, test_mode = 'individual', thresh_val = NULL,
                          sigma_gene = T, PRECISION.THRESHOLD = 0.01, cell_types_present = NULL,
                          add_res_genes = T, fdr = .01, params_to_test = NULL) {
  cell_types <- choose_cell_types(myRCTD, barcodes, doublet_mode, cell_type_threshold, cell_types)
  X1 <- check_designmatrix(X1, 'run.RCTDE.general')
  X2 <- check_designmatrix(X2, 'run.RCTDE.general', require_2d = TRUE)
  if(!(test_mode %in% c('individual', 'categorical')))
    stop(c('run.RCTDE.general: not valid test_mode = ',test_mode,'. Please set test_mode = "categorical" or "individual".'))
  if(is.null(params_to_test))
    if(test_mode == 'individual')
      params_to_test <- 2
    else
      params_to_test <- 1:dim(X2)[2]
  print(paste0("run.RCTDE.general: configure params_to_test = ",params_to_test))
  if(any(!(params_to_test %in% 1:dim(X2)[2])))
    stop(c('run.RCTDE.general: params_to_test must be a vector of integers from 1 to dim(X2)[2] = ', dim(X2)[2],
           'please make sure that tested parameters are in the required range.'))
  if(test_mode == 'categorical' && any(!(X2[,params_to_test] %in% c(0,1))))
    stop(c('run.RCTDE.general: for test_mode = categorical, colums params_to_test, ',params_to_test,', must have values 0 or 1.'))
  if(is.null(cell_types_present))
    cell_types_present <- cell_types
  if(doublet_mode && myRCTD@config$RCTDmode != 'doublet')
    stop('run.RCTDE.general: attempted to run RCTDE in doublet mode, but RCTD was not run in doublet mode. Please run RCTDE in full mode (doublet_mode = F) or run RCTD in doublet mode.')
  puck = myRCTD@originalSpatialRNA
  gene_list_tot <- filter_genes(puck, threshold = gene_threshold)
  nUMI <- puck@nUMI[barcodes]
  cell_type_info <- myRCTD@cell_type_info$info
  if(doublet_mode) {
    my_beta <- get_beta_doublet(barcodes, cell_type_info[[2]], myRCTD@results$results_df, myRCTD@results$weights_doublet)
    thresh = 0.999
  } else {
    my_beta <- as.matrix(sweep(myRCTD@results$weights, 1, rowSums(myRCTD@results$weights), '/'))
    thresh = 0.95
  }
  if(!is.null(thresh_val))
    thresh = thresh_val
  res <- filter_barcodes_cell_types(barcodes, cell_types, my_beta, thresh = thresh)
  barcodes <- res$barcodes; my_beta <- res$my_beta
  set_likelihood_vars(myRCTD@internal_vars$Q_mat, myRCTD@internal_vars$X_vals)
  if(sigma_gene)
    set_global_Q_all()
  sigma_init <- as.character(100*myRCTD@internal_vars$sigma)
  gene_fits <- get_de_gene_fits(X1[barcodes,],X2[barcodes,],my_beta, nUMI[barcodes], gene_list_tot,
                                cell_types, restrict_puck(puck, barcodes), barcodes, sigma_init,
                                test_mode, numCores = myRCTD@config$max_cores, sigma_gene = sigma_gene,
                                PRECISION.THRESHOLD = PRECISION.THRESHOLD, params_to_test = params_to_test)
  if(add_res_genes)
    res_gene_list <- get_res_genes(puck, myRCTD, gene_list_tot, cell_types, my_beta, barcodes, nUMI,
                            gene_fits, cell_types_present, X2, test_mode, fdr = fdr, params_to_test = params_to_test)
  else
    res_gene_list <- NULL
  myRCTD@internal_vars_de <- list(barcodes = barcodes, cell_types = cell_types, doublet_mode = doublet_mode,
                                  cell_types_present = cell_types_present,
                                  my_beta = my_beta, X1 = X1, X2 = X2,
                                  test_mode = test_mode, params_to_test = params_to_test)
  myRCTD@de_results <- list(gene_fits = gene_fits, res_gene_list = res_gene_list)
  return(myRCTD)
}

get_res_genes <- function(puck, myRCTD, gene_list_tot, cell_types, my_beta, barcodes, nUMI,
                          gene_fits, cell_types_present, X2, test_mode, params_to_test = 2,
                          fdr = .01, p_thresh = 1, log_fc_thresh = 0.4) {
  cti_renorm <- get_norm_ref(puck, myRCTD@cell_type_info$info[[1]], intersect(gene_list_tot,rownames(myRCTD@cell_type_info$info[[1]])), myRCTD@internal_vars$proportions)
  res_gene_list <- list()
  for(cell_type in cell_types) {
    gene_list_type <- get_gene_list_type(my_beta, barcodes, cell_type, nUMI, gene_list_tot,
                                         cti_renorm, cell_types_present, gene_fits, test_mode = test_mode)
    if(test_mode == 'individual')
      res_genes <- find_sig_genes_individual(cell_type, cell_types, gene_fits, gene_list_type, X2,
                                         params_to_test = params_to_test, fdr = fdr, p_thresh = p_thresh, log_fc_thresh = log_fc_thresh)
    else if(test_mode == 'categorical') {
      res_genes <- find_sig_genes_categorical(cell_type, cell_types, gene_fits, gene_list_type, X2,
                                        p_thresh = p_thresh, log_fc_thresh = log_fc_thresh,
                                        params_to_test = params_to_test)
    }
    res_gene_list[[cell_type]] <- res_genes
  }
  return(res_gene_list)
}

add_res_genes <- function(myRCTD,  params_to_test = NULL, fdr = .01, p_thresh = 1, log_fc_thresh = 0.4) {
  puck <- myRCTD@originalSpatialRNA
  gene_list_tot <- rownames(myRCTD@de_results$gene_fits$s_mat)
  cell_types <- myRCTD@internal_vars_de$cell_types
  my_beta <- myRCTD@internal_vars_de$my_beta
  X2 <- myRCTD@internal_vars_de$X2
  barcodes <- myRCTD@internal_vars_de$barcodes
  nUMI <- puck@nUMI[barcodes]
  test_mode <- myRCTD@internal_vars_de$test_mode
  cell_types_present <- myRCTD@internal_vars_de$cell_types_present
  gene_fits <- myRCTD@de_results$gene_fits
  if(is.null(params_to_test))
    params_to_test <- myRCTD@internal_vars_de$params_to_test
  res_gene_list <- get_res_genes(puck, myRCTD, gene_list_tot, cell_types, my_beta, barcodes, nUMI,
                                  gene_fits, cell_types_present, X2, test_mode,
                                  params_to_test = params_to_test, fdr = fdr,
                                 p_thresh = p_thresh, log_fc_thresh = log_fc_thresh)
  myRCTD@de_results$res_gene_list <- res_gene_list
  return(myRCTD)
}

find_sig_genes_categorical <- function(cell_type, cell_types, gene_fits, gene_list_type, X2, fdr = 0.01,
                                p_thresh = 1, log_fc_thresh = 0.4, params_to_test = NULL) {
  if(is.null(params_to_test))
    params_to_test <- 1:dim(X2)[2]
  n_regions <- length(params_to_test); n_cell_types <- length(cell_types)
  cell_ind = (which(cell_types == cell_type))
  s_mat_ind <- (1:dim(X2)[2]) + (n_regions*(cell_ind - 1))
  p_val_sig_pair <- numeric(length(gene_list_type)); names(p_val_sig_pair) <- gene_list_type
  log_fc_best_pair <- numeric(length(gene_list_type)); names(log_fc_best_pair) <- gene_list_type
  sd_vec <- numeric(length(gene_list_type)); names(sd_vec) <- gene_list_type
  sd_lfc_vec <- numeric(length(gene_list_type)); names(sd_lfc_vec) <- gene_list_type
  i1_vec <- numeric(length(gene_list_type)); names(i1_vec) <- gene_list_type
  i2_vec <- numeric(length(gene_list_type)); names(i2_vec) <- gene_list_type
  for(gene in gene_list_type) {
    con_regions <- get_con_regions(gene_fits, gene, X2, cell_ind, n_cell_types) &
      (params_to_test %in% 1:dim(X2)[2])
    n_regions_con <- sum(con_regions)
    x <- gene_fits$all_vals[gene, con_regions,cell_ind]
    s_mat_ind_cur <- s_mat_ind[con_regions]
    var_vals <- (gene_fits$s_mat[gene, s_mat_ind_cur])^2
    ovr_best_p_val <- 1
    best_log_fc <- 0; best_sd <- 0
    best_i1 <- 0; best_i2 <- 0
    for(i1 in 1:(n_regions_con-1))
      for(i2 in (i1+1):n_regions_con) {
        log_fc <- abs(x[i1] - x[i2])
        sd_cur <- sqrt(var_vals[i1] + var_vals[i2])
        z_score <- (log_fc) / sd_cur
        p_val <- 2*(1-pnorm(z_score))
        if(p_val < ovr_best_p_val) {
          ovr_best_p_val <- p_val
          best_log_fc <- log_fc
          best_sd <- sd_cur
          best_i1 <- which(con_regions)[i1]
          best_i2 <- which(con_regions)[i2]
        }
      }
    p_val_sig_pair[gene] <- min(1, ovr_best_p_val * choose(n_regions_con, 2))
    log_fc_best_pair[gene] <- best_log_fc
    sd_vec[gene] <- best_sd
    sd_lfc_vec[gene] <- sd(x)
    i1_vec[gene] <- best_i1
    i2_vec[gene] <- best_i2
  }
  gene_list_sig <- fdr_sig_genes(gene_list_type, p_val_sig_pair, fdr)
  res_genes <- data.frame(sd_lfc_vec[gene_list_sig], i1_vec[gene_list_sig], i2_vec[gene_list_sig],
                          sd_vec[gene_list_sig], p_val_sig_pair[gene_list_sig],
                          log_fc_best_pair[gene_list_sig])
  rownames(res_genes) <- gene_list_sig
  custom_names <- c('sd_lfc','sd_best','p_val_best','log_fc_best', 'paramindex1_best', 'paramindex2_best')
  colnames(res_genes) <- custom_names
  res_genes <- res_genes[abs(res_genes$p_val < p_thresh) & abs(res_genes$log_fc) >= log_fc_thresh, ]
  if(length(gene_list_sig) > 1)
    res_genes <- data.frame(res_genes, gene_fits$all_vals[rownames(res_genes),params_to_test,cell_ind]) # add on the means
  else {
    if(length(gene_list_sig) == 1) {
      res_genes <- data.frame(t(unlist((c(res_genes, gene_fits$all_vals[rownames(res_genes),params_to_test,cell_ind],
                                          gene_fits$s_mat[rownames(res_genes),s_mat_ind[params_to_test]])))))
      rownames(res_genes) <- gene_list_sig
      colnames(res_genes)[(length(custom_names)+1):length(res_genes)] <-
        c(lapply(params_to_test,function(x) paste('mean_',x)), lapply(params_to_test,function(x) paste('sd_',x)))
    } else
      res_genes <- list()
  }
  return(res_genes)
}

find_sig_genes_individual <- function(cell_type, cell_types, gene_fits, gene_list_type, X2, params_to_test = 2, fdr = 0.01, p_thresh = 1, log_fc_thresh = 0.4) {
  ct_ind <- which(cell_types == cell_type)
  I_ind = dim(X2)[2]*(ct_ind - 1) + params_to_test
  log_fc <- gene_fits$all_vals[gene_list_type,params_to_test, ct_ind]
  s_vec <- gene_fits$s_mat[gene_list_type,I_ind]
  z_score <- abs(log_fc) / s_vec
  if(length(params_to_test) > 1)
    p_val <- pmin(apply(p_val, 1, min)*length(params_to_test),1)
  else
    p_val <- 2*(1-pnorm(z_score))
  names(p_val) <- gene_list_type
  gene_list_sig <- fdr_sig_genes(gene_list_type, p_val, fdr)
  if(length(gene_list_sig) > 0)
    p_thresh <- min(p_thresh, max(p_val[gene_list_sig]))
  if(length(params_to_test) > 1) {
    best_mat <- function(gene) {
      index <- which(p_val[gene,] < p_thresh)
      if(length(index) > 0) {
        best_ind <- which.max(abs(z_score[gene,index]))
        lfc <- log_fc[gene,index][best_ind]
        sd <- s_vec[gene,index][best_ind]
        z <- z_score[gene,index][best_ind]
        return(best_ind, lfc,sd,z)
      } else {
        return(c(0,0,0,0))
      }
    }
    best_ind <- function(gene) {
      best_mat(gene)[1]
    }
    best_log_fc <- function(gene) {
      best_mat(gene)[2]
    }
    best_sd <- function(gene) {
      best_mat(gene)[3]
    }
    best_Z <- function(gene) {
      best_mat(gene)[4]
    }
    best_ind <- unlist(lapply(gene_list_type, best_ind))
    names(best_ind) <- gene_list_type
    log_fc <- unlist(lapply(gene_list_type, best_log_fc))
    names(log_fc) <- gene_list_type
    z_score <- unlist(lapply(gene_list_type, best_Z))
    names(z_score) <- gene_list_type
    s_vec <- unlist(lapply(gene_list_type, best_sd))
    names(s_vec) <- gene_list_type
    p_val <- pmin(apply(p_val, 1, min)*length(params_to_test),1)
  } else {
    best_ind <- rep(params_to_test, length(gene_list_type))
    names(best_ind) <- gene_list_type
  }
  res_genes <- data.frame(z_score[gene_list_sig], log_fc[gene_list_sig], s_vec[gene_list_sig], best_ind[gene_list_sig])
  names(res_genes) <- c('Z_score','log_fc', 'se', 'paramindex_best')
  res_genes$conv <- gene_fits$con_mat[gene_list_sig, cell_type]
  res_genes$p_val <- p_val[gene_list_sig]
  res_genes <- res_genes[abs(res_genes$p_val < p_thresh) & abs(res_genes$log_fc) >= log_fc_thresh, ]
  return(res_genes)
}

get_de_gene_fits <- function(X1,X2,my_beta, nUMI, gene_list, cell_types, puck, barcodes, sigma_init, test_mode,
                             numCores = 4, sigma_gene = T, PRECISION.THRESHOLD = 0.01, params_to_test = 2) {
  results_list <- fit_de_genes(X1,X2,my_beta, nUMI, gene_list, puck, barcodes, sigma_init, test_mode, numCores = numCores, sigma_gene = sigma_gene, PRECISION.THRESHOLD = PRECISION.THRESHOLD)
  N_genes <- length(results_list)
  intercept_val <- matrix(0,nrow = N_genes, ncol = length(cell_types))
  mean_val <- matrix(0,nrow = N_genes, ncol = length(cell_types))
  all_vals <- array(0, dim = c(N_genes, dim(X2)[2],length(cell_types)))
  dimnames(all_vals)[[1]] <- gene_list
  con_val <- logical(N_genes)
  ll_val <- numeric(N_genes)
  n_val <- numeric(N_genes)
  sigma_g <- numeric(N_genes)
  names(sigma_g) <- gene_list
  I_val <- list()
  names(n_val) <- gene_list
  names(con_val) <- gene_list
  names(ll_val) <- gene_list
  rownames(mean_val) <- gene_list; colnames(mean_val) <- cell_types
  rownames(intercept_val) <- gene_list; colnames(intercept_val) <- cell_types
  d_vals <- matrix(0,nrow=N_genes,ncol=dim(X2)[2]*length(cell_types))
  s_mat <- matrix(0, nrow = N_genes, ncol = dim(X2)[2]*length(cell_types))
  precision_mat <- matrix(0, nrow = N_genes, ncol = dim(X2)[2]*length(cell_types))
  con_all <- matrix(FALSE, nrow = N_genes, ncol = dim(X2)[2]*length(cell_types))
  con_mat <- matrix(FALSE, nrow = N_genes, ncol = length(cell_types))
  error_mat <- matrix(FALSE, nrow = N_genes, ncol = length(cell_types))
  rownames(precision_mat) <- gene_list; rownames(con_all) <- gene_list
  rownames(s_mat) <- gene_list; rownames(con_mat) <- gene_list; rownames(error_mat) <- gene_list
  colnames(s_mat) <- get_param_names(X1,X2, cell_types)
  colnames(precision_mat) <- get_param_names(X1,X2, cell_types)
  colnames(con_all) <- get_param_names(X1,X2, cell_types)
  colnames(con_mat) <- cell_types
  colnames(error_mat) <- cell_types
  rownames(d_vals) <- gene_list
  for(i in 1:N_genes) {
    sigma_g[i] <- results_list[[i]]$sigma_s_best
    res <- results_list[[i]]$res
    d_vals[i,] <- res$d
    mean_val[i,] <- res$alpha2[params_to_test[1],]
    intercept_val[i,] <- res$alpha2[1,]
    all_vals[i, ,] <- res$alpha2
    con_val[i] <- res$converged
    precision_mat[i,] <- res$precision
    ll_val[i] <- res$log_l
    n_val[i] <- res$n.iter
    I_val[[i]] <- res$I
    s_mat[i,] <- sqrt(diag(I_val[[i]]))
    con_mat[i,] <- res$converged_vec
    con_all[i,] <- res$precision < PRECISION.THRESHOLD
    error_mat[i,] <- res$error_vec
  }
  return(list(mean_val = mean_val, con_val = con_val, ll_val = ll_val, I_val = I_val, s_mat = s_mat,
              n.iter = n_val,d_vals = d_vals, intercept_val = intercept_val, all_vals = all_vals,
              precision_mat = precision_mat, sigma_g = sigma_g, con_mat = con_mat, con_all = con_all, error_mat = error_mat))
}

fit_de_genes <- function(X1,X2,my_beta, nUMI, gene_list, puck, barcodes, sigma_init, test_mode, numCores = 4, sigma_gene = T, PRECISION.THRESHOLD = 0.01) {
  results_list <- list()
  if(numCores == 1) {
    for(i in 1:length(gene_list)) {
      print(i)
      gene <- gene_list[i]
      Y <- puck@counts[gene, barcodes]
      results_list[[i]] <- estimate_gene_wrapper(Y,X1,X2,my_beta, nUMI, sigma_init, test_mode, verbose = F, n.iter = 200, MIN_CHANGE = 1e-3, sigma_gene = sigma_gene, PRECISION.THRESHOLD = PRECISION.THRESHOLD)
    }
  } else {
    cl <- parallel::makeCluster(numCores,setup_strategy = "sequential",outfile="") #makeForkCluster
    doParallel::registerDoParallel(cl)
    environ = c('estimate_effects_trust', 'solveIRWLS.effects_trust', 'Q_mat', 'K_val','X_vals',
                'calc_log_l_vec', 'calc_Q_k','get_d1_d2', 'calc_Q_all','psd','construct_hess_fast',
                'choose_sigma_gene', 'estimate_gene_wrapper', 'check_converged_vec')
    if(sigma_gene)
      environ <- c(environ, 'Q_mat_all')
    out_file = "logs/de_log.txt"
    if(!dir.exists('logs'))
      dir.create('logs')
    if(file.exists(out_file))
      file.remove(out_file)
    results_list <- foreach::foreach(i = 1:length(gene_list), .packages = c("quadprog", "RCTD"), .export = environ) %dopar% {
      if(i %% 10 == 0) {
        cat(paste0("Finished sample: ",i," gene ", gene_list[i],"\n"), file=out_file, append=TRUE)
      }
      assign("Q_mat",Q_mat, envir = globalenv()); assign("X_vals",X_vals, envir = globalenv())
      assign("K_val",K_val, envir = globalenv());
      if(sigma_gene)
        assign("Q_mat_all",Q_mat_all, envir = globalenv());
      gene <- gene_list[i]
      Y <- puck@counts[gene, barcodes]
      res <- estimate_gene_wrapper(Y,X1,X2,my_beta, nUMI, sigma_init, test_mode, verbose = F, n.iter = 200, MIN_CHANGE = 1e-3, sigma_gene = sigma_gene)
    }
    parallel::stopCluster(cl)
  }
  return(results_list)
}
