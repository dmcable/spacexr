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
#' @param test_mode (default 'direct') if 'direct', tests for DE along a single parameter. If 'multi', then tests for differences
#' across multiple parameters
#' @param delta (default 0) non-centrality parameter for differential expression tests
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
run.de.single <- function(myRCTD, explanatory.variable, datadir = './', cell_types = NULL, cell_type_threshold = 125,
                          gene_threshold = 5e-5, doublet_mode = T, test_mode = 'direct', delta = 0, thresh_val = NULL,
                          sigma_gene = T, PRECISION.THRESHOLD = 0.01, cell_types_present = NULL, fdr = .01) {
  check_vector(explanatory.variable, 'explanatory.variable', 'run.de.single')
  all_barc = intersect(names(explanatory.variable), colnames(myRCTD@spatialRNA@counts))
  if(length(all_barc) <= 1)
    stop(paste0('run.de.single: ', length(all_barc),
    ' common barcode names found between explanatory.variable and myRCTD@spatialRNA. Please ensure that more common barcodes are found'))
  explanatory.variable <- explanatory.variable[all_barc]
  if(max(explanatory.variable) - min(explanatory.variable) < 1e-9)
    stop('run.de.single: values of explanatory.variable are constant. Please increase the range of this variable.')
  explanatory.variable <- explanatory.variable - min(explanatory.variable)
  explanatory.variable <- explanatory.variable / max(explanatory.variable) # standardize
  X2 <- cbind(1,explanatory.variable)
  X1 <- matrix(0,nrow = length(all_barc),ncol=0)
  rownames(X2) <- all_barc; rownames(X1) <- all_barc
  cell_types <- choose_cell_types(myRCTD, all_barc, doublet_mode, cell_type_threshold, cell_types)
  return(find_de_genes(myRCTD, X1, X2, all_barc, cell_types, datadir, gene_threshold = gene_threshold,
                       doublet_mode = doublet_mode, test_mode = test_mode, delta = delta,
                       thresh_val = thresh_val, sigma_gene = sigma_gene,
                       PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                       cell_types_present = cell_types_present, fdr = fdr))
}

# run DE nonparametrically
run.de.nonparam <- function(myRCTD, df = 15, datadir = './', all_barc = NULL, cell_types = NULL,
                            cell_type_threshold = 125, gene_threshold = 5e-5, doublet_mode = T,
                            test_mode = 'direct', delta = 0, thresh_val = NULL, sigma_gene = T,
                            PRECISION.THRESHOLD = 0.01, cell_types_present = NULL) {
  if(is.null(all_barc))
    all_barc <- colnames(myRCTD@spatialRNA@counts)
  else
    all_barc = intersect(all_barc, colnames(myRCTD@spatialRNA@counts))
  if(length(all_barc) <= 1)
    stop(paste0('run.de.nonparam: ', length(all_barc),
                ' common barcode names found between all_barc and myRCTD@spatialRNA. Please ensure that more common barcodes are found'))
  X2 <- get_spline_matrix(myRCTD@spatialRNA, df = df)
  X1 <- matrix(0,nrow = length(all_barc),ncol=0)
  rownames(X2) <- all_barc; rownames(X1) <- all_barc
  cell_types <- choose_cell_types(myRCTD, all_barc, doublet_mode, cell_type_threshold, cell_types)
  return(find_de_genes(myRCTD, X1, X2, all_barc, cell_types, datadir, gene_threshold = gene_threshold,
                       doublet_mode = doublet_mode, test_mode = test_mode, delta = delta,
                       thresh_val = thresh_val, sigma_gene = sigma_gene,
                       PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                       cell_types_present = cell_types_present, add_res_genes = F))
}

run.de.regions <- function(myRCTD, region_list, datadir = './', cell_types = NULL,
                           cell_type_threshold = 125, gene_threshold = 5e-5, doublet_mode = T,
                           test_mode = 'multi', delta = 0, thresh_val = NULL, sigma_gene = T,
                           PRECISION.THRESHOLD = 0.01, cell_types_present = NULL) {
  if(!is.list(region_list))
    stop('run.de.regions: error, region_list must be a list')
  n_regions <- length(region_list)
  if(n_regions < 3)
    stop('run.de.regions: length(region_list) <= 2. Must be at least 3 to continue.')
  for(i in 1:n_regions) {
    barcodes <- region_list[[i]]
    if(!is.character(barcodes) || !is.atomic(barcodes))
      stop('run.de.regions: error, region_list must be a list of atomic character vectors')
    shorter_barcodes <- intersect(barcodes, colnames(myRCTD@spatialRNA@counts))
    if(length(barcodes) > length(shorter_barcodes))
      warning('run.de.regions: some barcodes from region_list not found in RCTD object')
    region_list[[i]] <- shorter_barcodes
    if(length(shorter_barcodes) < 2)
      stop('run.de.regions: error, region_list must be a list of atomic character vectors of length at least 2.')
  }
  all_barc <- Reduce(union, region_list)
  X2 <- matrix(0, nrow = length(all_barc), ncol = n_regions)
  X1 <- matrix(0,nrow = length(all_barc),ncol=0)
  rownames(X2) <- all_barc; rownames(X1) <- all_barc
  for(i in 1:n_regions)
    X2[region_list[[i]],i] <- 1
  cell_types <- choose_cell_types(myRCTD, all_barc, doublet_mode, cell_type_threshold, cell_types)
  return(find_de_genes(myRCTD, X1, X2, all_barc, cell_types, datadir, gene_threshold = gene_threshold,
                       doublet_mode = doublet_mode, test_mode = test_mode, delta = delta,
                       thresh_val = thresh_val, sigma_gene = sigma_gene,
                       PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                       cell_types_present = cell_types_present))

}

find_de_genes <- function(myRCTD, X1, X2, all_barc, cur_cell_types, datadir = './', gene_threshold = 5e-5,
                          doublet_mode = T, test_mode = 'direct', delta = 0, thresh_val = NULL,
                          sigma_gene = T, PRECISION.THRESHOLD = 0.01, cell_types_present = NULL,
                          add_res_genes = T, fdr = .01) {
  if(!(test_mode %in% c('direct', 'multi')))
    stop(c('find_de_genes: not implemented for test_mode = ',test_mode,'. Please set test_mode = "multi" or "direct".'))
  if(is.null(cell_types_present))
    cell_types_present <- cur_cell_types
  puck = myRCTD@originalSpatialRNA
  gene_list_tot <- filter_genes(puck, threshold = gene_threshold)
  nUMI <- puck@nUMI[all_barc]
  cell_type_info <- myRCTD@cell_type_info$info
  if(doublet_mode) {
    my_beta <- get_beta_doublet(all_barc, cell_type_info[[2]], myRCTD@results$results_df, myRCTD@results$weights_doublet)
    thresh = 0.999
  } else {
    my_beta <- as.matrix(sweep(myRCTD@results$weights, 1, rowSums(myRCTD@results$weights), '/'))
    thresh = 0.95
  }
  if(!is.null(thresh_val))
    thresh = thresh_val
  res <- filter_barcodes_cell_types(all_barc, cur_cell_types, my_beta, thresh = thresh)
  all_barc <- res$all_barc; my_beta <- res$my_beta
  set_likelihood_vars(myRCTD@internal_vars$Q_mat, myRCTD@internal_vars$X_vals)
  if(sigma_gene)
    set_global_Q_all()
  sigma_init <- as.character(100*myRCTD@internal_vars$sigma)
  gene_fits <- get_de_gene_fits(X1[all_barc,],X2[all_barc,],my_beta, nUMI[all_barc], gene_list_tot,
                                cur_cell_types, restrict_puck(puck, all_barc), all_barc, sigma_init,
                                test_mode, numCores = myRCTD@config$max_cores, sigma_gene = sigma_gene,
                                PRECISION.THRESHOLD = PRECISION.THRESHOLD)
  if(add_res_genes)
    res_gene_list <- get_res_genes(puck, myRCTD, gene_list_tot, cur_cell_types, my_beta, all_barc, nUMI,
                            gene_fits, cell_types_present, X2, test_mode, delta, datadir, fdr = fdr)
  else
    res_gene_list <- NULL
  myRCTD@internal_vars_de <- list(all_barc = all_barc, cell_types = cur_cell_types,
                                  cell_types_present = cell_types_present,
                                  my_beta = my_beta, X1 = X1, X2 = X2, delta = delta,
                                  test_mode = test_mode)
  myRCTD@de_results <- list(gene_fits = gene_fits, res_gene_list = res_gene_list)
  return(myRCTD)
}

get_res_genes <- function(puck, myRCTD, gene_list_tot, cur_cell_types, my_beta, all_barc, nUMI,
                          gene_fits, cell_types_present, X2, test_mode, delta, datadir,
                          plot_genes = T, param_position = 2, fdr = .01, p_thresh = 1, log_fc_thresh = 0.4) {
  cti_renorm <- get_norm_ref(puck, myRCTD@cell_type_info$info[[1]], intersect(gene_list_tot,rownames(myRCTD@cell_type_info$info[[1]])), myRCTD@internal_vars$proportions)
  res_gene_list <- list()
  for(cell_type in cur_cell_types) {
    gene_list_type <- get_gene_list_type(my_beta, all_barc, cell_type, nUMI, gene_list_tot,
                                         cti_renorm, cell_types_present, gene_fits, test_mode = test_mode)
    if(test_mode == 'direct')
      res_genes <- find_sig_genes_direct(cell_type, cur_cell_types, gene_fits, gene_list_type, X2,
                                         param_position = param_position, fdr = fdr, p_thresh = p_thresh, log_fc_thresh = log_fc_thresh)
    else if(test_mode == 'multi') {
      res_genes <- find_sig_genes_multi(cell_type, cur_cell_types, gene_fits, gene_list_type, X2,
                                        delta = delta, p_thresh = p_thresh, log_fc_thresh = log_fc_thresh)
    }
    if(plot_genes)
      plot_sig_genes(cell_type,all_barc,my_beta,puck,res_genes, test_mode, datadir)
    res_gene_list[[cell_type]] <- res_genes
  }
  return(res_gene_list)
}

add_res_genes <- function(myRCTD, datadir = './', plot_genes = T, param_position = 2, fdr = .01,
                          p_thresh = 1, log_fc_thresh = 0.4) {
  puck <- myRCTD@originalSpatialRNA
  gene_list_tot <- rownames(myRCTD@de_results$gene_fits$I_mat)
  cur_cell_types <- myRCTD@internal_vars_de$cell_types
  my_beta <- myRCTD@internal_vars_de$my_beta
  X2 <- myRCTD@internal_vars_de$X2
  all_barc <- myRCTD@internal_vars_de$all_barc
  nUMI <- puck@nUMI[all_barc]
  delta <- myRCTD@internal_vars_de$delta
  test_mode <- myRCTD@internal_vars_de$test_mode
  cell_types_present <- myRCTD@internal_vars_de$cell_types_present
  gene_fits <- myRCTD@de_results$gene_fits
  res_gene_list <- get_res_genes(puck, myRCTD, gene_list_tot, cur_cell_types, my_beta, all_barc, nUMI,
                                  gene_fits, cell_types_present, X2, test_mode, delta, datadir,
                                 plot_genes = plot_genes, param_position = param_position, fdr = fdr,
                                 p_thresh = p_thresh, log_fc_thresh = log_fc_thresh)
  myRCTD@de_results$res_gene_list <- res_gene_list
  return(myRCTD)
}

find_sig_genes_multi <- function(cell_type, cur_cell_types, gene_fits, gene_list_type, X2, fdr = 0.01,
                                p_thresh = 1, log_fc_thresh = 0.4, delta = 0) {
  n_regions <- dim(X2)[2]; n_cell_types <- length(cur_cell_types)
  cell_ind = (which(cur_cell_types == cell_type))
  I_mat_ind <- (1:n_regions) + (n_regions*(cell_ind - 1))
  p_val_sig_pair <- numeric(length(gene_list_type)); names(p_val_sig_pair) <- gene_list_type
  log_fc_best_pair <- numeric(length(gene_list_type)); names(log_fc_best_pair) <- gene_list_type
  sd_vec <- numeric(length(gene_list_type)); names(sd_vec) <- gene_list_type
  sd_lfc_vec <- numeric(length(gene_list_type)); names(sd_lfc_vec) <- gene_list_type
  for(gene in gene_list_type) {
    con_regions <- get_con_regions(gene_fits, gene, n_regions, cell_ind, n_cell_types)
    n_regions_con <- sum(con_regions)
    x <- gene_fits$all_vals[gene, con_regions,cell_ind]
    I_mat_ind_cur <- I_mat_ind[con_regions]
    var_vals <- (gene_fits$I_mat[gene, I_mat_ind_cur])^2
    ovr_best_p_val <- 1
    best_log_fc <- 0
    best_sd <- 0
    for(i1 in 1:(n_regions_con-1))
      for(i2 in (i1+1):n_regions_con) {
        log_fc <- abs(x[i1] - x[i2])
        sd_cur <- sqrt(var_vals[i1] + var_vals[i2])
        z_score <- (log_fc - delta) / sd_cur
        p_val <- 2*(1-pnorm(z_score))
        if(p_val < ovr_best_p_val) {
          ovr_best_p_val <- p_val
          best_log_fc <- log_fc
          best_sd <- sd_cur
        }
      }
    p_val_sig_pair[gene] <- min(1, ovr_best_p_val * choose(n_regions_con, 2))
    log_fc_best_pair[gene] <- best_log_fc
    sd_vec[gene] <- best_sd
    sd_lfc_vec[gene] <- sd(x)
  }
  gene_list_sig <- fdr_sig_genes(gene_list_type, p_val_sig_pair, fdr)
  res_genes <- data.frame(sd_vec[gene_list_sig], sd_lfc_vec[gene_list_sig],
                          p_val_sig_pair[gene_list_sig], log_fc_best_pair[gene_list_sig])
  rownames(res_genes) <- gene_list_sig
  colnames(res_genes) <- c('sd','sd_lfc','p_val','log_fc')
  res_genes <- res_genes[abs(res_genes$p_val < p_thresh) & abs(res_genes$log_fc) >= log_fc_thresh, ]
  if(length(gene_list_sig) > 1)
    res_genes <- data.frame(res_genes, gene_fits$all_vals[rownames(res_genes),,cell_ind]) # add on the means
  else {
    if(length(gene_list_sig) == 1) {
      res_genes <- data.frame(t(unlist((c(res_genes, gene_fits$all_vals[rownames(res_genes),,cell_ind])))))
      rownames(res_genes) <- gene_list_sig
      colnames(res_genes)[5:length(res_genes)] <- 1:n_regions
    } else
      res_genes <- list()
  }
  return(res_genes)
}

find_sig_genes_direct <- function(cell_type, cur_cell_types, gene_fits, gene_list_type, X2, param_position = 2, fdr = 0.01, p_thresh = 1, log_fc_thresh = 0.4) {
  ct_ind <- which(cur_cell_types == cell_type)
  I_ind = dim(X2)[2]*(ct_ind - 1) + param_position
  log_fc <- gene_fits$all_vals[gene_list_type,param_position, ct_ind]
  z_score <- abs(log_fc) / gene_fits$I_mat[gene_list_type,I_ind]
  p_val <- 2*(1-pnorm(z_score))
  if(length(param_position) > 1) {
    best_mat <- function(log_fc_mat, gene) {
      index <- which(p_val[gene,] < p_thresh)
      if(length(index) > 0)
        log_fc <- max(abs(log_fc_mat[gene,index]))
      else
        log_fc <- 0
      return(log_fc)
    }
    best_log_fc <- function(gene) {
      best_mat(log_fc, gene)
    }
    best_Z <- function(gene) {
      best_mat(z_score, gene)
    }
    log_fc <- unlist(lapply(gene_list_type, best_log_fc))
    names(log_fc) <- gene_list_type
    z_score <- unlist(lapply(gene_list_type, best_Z))
    names(z_score) <- gene_list_type
    p_val <- pmin(apply(p_val, 1, min)*length(param_position),1)
  }
  names(p_val) <- gene_list_type
  gene_list_sig <- fdr_sig_genes(gene_list_type, p_val, fdr)
  res_genes <- data.frame(z_score[gene_list_sig], log_fc[gene_list_sig])
  names(res_genes) <- c('Z_score','log_fc')
  res_genes$conv <- gene_fits$con_mat[gene_list_sig, cell_type]
  #res_genes$d <- gene_fits$d_vals[gene_list_sig,I_ind]
  #res_genes$precision <- gene_fits$precision[gene_list_sig]
  res_genes$p_val <- p_val[gene_list_sig]
  res_genes <- res_genes[abs(res_genes$p_val < p_thresh) & abs(res_genes$log_fc) >= log_fc_thresh, ]
  return(res_genes)
}

get_de_gene_fits <- function(X1,X2,my_beta, nUMI, gene_list, cur_cell_types, puck, all_barc, sigma_init, test_mode, numCores = 4, sigma_gene = T, PRECISION.THRESHOLD = 0.01) {
  results_list <- fit_de_genes(X1,X2,my_beta, nUMI, gene_list, puck, all_barc, sigma_init, test_mode, numCores = numCores, sigma_gene = sigma_gene, PRECISION.THRESHOLD = PRECISION.THRESHOLD)
  N_genes <- length(results_list)
  base_val <- matrix(0,nrow = N_genes, ncol = length(cur_cell_types))
  mean_val <- matrix(0,nrow = N_genes, ncol = length(cur_cell_types))
  all_vals <- array(0, dim = c(N_genes, dim(X2)[2],length(cur_cell_types)))
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
  rownames(mean_val) <- gene_list; colnames(mean_val) <- cur_cell_types
  rownames(base_val) <- gene_list; colnames(base_val) <- cur_cell_types
  d_vals <- matrix(0,nrow=N_genes,ncol=dim(X2)[2]*length(cur_cell_types))
  I_mat <- matrix(0, nrow = N_genes, ncol = dim(X2)[2]*length(cur_cell_types))
  p_mat <- matrix(0, nrow = N_genes, ncol = dim(X2)[2]*length(cur_cell_types))
  con_all <- matrix(FALSE, nrow = N_genes, ncol = dim(X2)[2]*length(cur_cell_types))
  con_mat <- matrix(FALSE, nrow = N_genes, ncol = length(cur_cell_types))
  error_mat <- matrix(FALSE, nrow = N_genes, ncol = length(cur_cell_types))
  rownames(p_mat) <- gene_list; rownames(con_all) <- gene_list
  rownames(I_mat) <- gene_list; rownames(con_mat) <- gene_list; rownames(error_mat) <- gene_list
  colnames(I_mat) <- get_param_names(X1,X2, cur_cell_types)
  colnames(p_mat) <- get_param_names(X1,X2, cur_cell_types)
  colnames(con_all) <- get_param_names(X1,X2, cur_cell_types)
  colnames(con_mat) <- cur_cell_types
  colnames(error_mat) <- cur_cell_types
  rownames(d_vals) <- gene_list
  for(i in 1:N_genes) {
    sigma_g[i] <- results_list[[i]]$sigma_s_best
    res <- results_list[[i]]$res
    d_vals[i,] <- res$d
    mean_val[i,] <- res$alpha2[2,]
    base_val[i,] <- res$alpha2[1,]
    all_vals[i, ,] <- res$alpha2
    con_val[i] <- res$converged
    p_mat[i,] <- res$precision
    ll_val[i] <- res$log_l
    n_val[i] <- res$n.iter
    I_val[[i]] <- res$I
    I_mat[i,] <- sqrt(diag(I_val[[i]]))
    con_mat[i,] <- res$converged_vec
    con_all[i,] <- res$precision < PRECISION.THRESHOLD
    error_mat[i,] <- res$error_vec
  }
  return(list(mean_val = mean_val, con_val = con_val, ll_val = ll_val, I_val = I_val, I_mat = I_mat,
              n.iter = n_val,d_vals = d_vals, base_val = base_val, all_vals = all_vals,
              p_mat = p_mat, sigma_g = sigma_g, con_mat = con_mat, con_all = con_all, error_mat = error_mat))
}

fit_de_genes <- function(X1,X2,my_beta, nUMI, gene_list, puck, all_barc, sigma_init, test_mode, numCores = 4, sigma_gene = T, PRECISION.THRESHOLD = 0.01) {
  results_list <- list()
  if(numCores == 1) {
    for(i in 1:length(gene_list)) {
      print(i)
      gene <- gene_list[i]
      Y <- puck@counts[gene, all_barc]
      results_list[[i]] <- estimate_gene_wrapper(Y,X1,X2,my_beta, nUMI, sigma_init, test_mode, verbose = F, n.iter = 200, MIN_CHANGE = 1e-3, sigma_gene = sigma_gene, PRECISION.THRESHOLD = PRECISION.THRESHOLD)
    }
  } else {
    cl <- parallel::makeCluster(numCores,setup_strategy = "sequential",outfile="") #makeForkCluster
    doParallel::registerDoParallel(cl)
    environ = c('estimate_effects_trust', 'solveIRWLS.effects_trust', 'Q_mat', 'K_val','X_vals','delta',
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
      Y <- puck@counts[gene, all_barc]
      res <- estimate_gene_wrapper(Y,X1,X2,my_beta, nUMI, sigma_init, test_mode, verbose = F, n.iter = 200, MIN_CHANGE = 1e-3, sigma_gene = sigma_gene)
    }
    parallel::stopCluster(cl)
  }
  return(results_list)
}
