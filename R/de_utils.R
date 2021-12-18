filter_genes <- function(puck, threshold = 5e-5) {
  norm_counts <- sweep(puck@counts,2,puck@nUMI,'/')
  gene_list_tot <- names(which(rowMeans(norm_counts) > threshold))
  if(length(grep("mt-",gene_list_tot)) > 0)
    gene_list_tot <- gene_list_tot[-grep("mt-",gene_list_tot)]
  return(gene_list_tot)
}

get_beta_doublet <- function(barcodes, cell_type_names, results_df, weights_doublet) {
  my_beta <- matrix(0, nrow = length(barcodes), ncol = length(cell_type_names))
  rownames(my_beta) <- barcodes
  colnames(my_beta) <- cell_type_names
  for(barcode in barcodes) {
    if(results_df[barcode, 'spot_class'] %in% c('singlet', 'doublet_certain')) {
      my_beta[barcode, results_df[barcode,'first_type']] <- weights_doublet[barcode,1]
      if(results_df[barcode, 'spot_class'] == "doublet_certain")
        my_beta[barcode, results_df[barcode,'second_type']] <- weights_doublet[barcode,2]
      else
        my_beta[barcode, results_df[barcode,'first_type']] <- 1
    }
  }
  return(my_beta)
}

filter_barcodes_cell_types <- function(barcodes, cell_types, my_beta, thresh = 0.9999) {
  barcodes <- barcodes[(rowSums(my_beta[barcodes, cell_types]) >= thresh)]
  my_beta <- my_beta[barcodes,cell_types]
  return (list(barcodes = barcodes, my_beta = my_beta))
}

get_p_chi <- function(x2, d, S_inv, points_list, delta = 0 ) {
  lambda <- get_l_chi(x2, d, S_inv, points_list, delta = delta)
  return(1 - pchisq(x2, d, ncp = lambda))
}

get_l_chi <- function(x2, d, S_inv, points_list, delta = 0 ) {
  lambda <- 0
  if(delta < 1e-6)
    lambda <- 0
  else {
    max_found <- 0
    for(ind in 1:(2^d - 1)) {
      new_val <- t(points_list[,ind]) %*% (S_inv %*% points_list[,ind])
      if(new_val > max_found)
        max_found <- new_val
    }
    lambda <- (max_found * (delta ^ 2))
  }
  return(lambda)
}

get_gene_list_type <- function(my_beta, barcodes, cell_type, nUMI, gene_list_type, cti_renorm,
                               cell_types_present, gene_fits, test_mode = 'individual') {
  C = 15
  N_cells <- colSums(my_beta[barcodes,])[cell_type]
  UMI_list <- nUMI[names(which(my_beta[barcodes,cell_type] >= .99))]
  if(length(UMI_list) < 10)
    UMI_list <- nUMI[names(which(my_beta[barcodes,cell_type] >= .80))]
  UMI_m <- median(UMI_list)
  expr_thresh <-  C / (N_cells * UMI_m)
  gene_list_type <- setdiff(gene_list_type,gene_list_type[which(cti_renorm[gene_list_type,cell_type] < expr_thresh)])
  cell_type_means <- cti_renorm[gene_list_type,cell_types_present]
  cell_prop <- sweep(cell_type_means,1,apply(cell_type_means,1,max),'/')
  gene_list_type <- gene_list_type[which(cell_prop[gene_list_type,cell_type] > 0.5)]
  if(test_mode == 'categorical') {
    n_cell_types <- dim(my_beta)[2]
    n_regions <- dim(gene_fits$con_all)[2] / n_cell_types
    cell_type_ind <- which(colnames(my_beta) == cell_type)
    my_filter <- unlist(lapply(gene_list_type, function(gene)
      (sum(get_con_regions(gene_fits, gene, n_regions, cell_type_ind, n_cell_types)) >= 2)))
    gene_list_type <- gene_list_type[my_filter]
  } else
    gene_list_type <- intersect(gene_list_type, names(which(gene_fits$con_mat[,cell_type]))) #only take converged genes
  return(gene_list_type)
}

get_con_regions <- function(gene_fits, gene, X2, cell_type_ind, n_cell_types) {
  matrix(gene_fits$con_all[gene,], nrow = dim(X2)[2], ncol = n_cell_types)[,cell_type_ind]
}

get_gene_list_type_wrapper <- function(myRCTD, cell_type, cell_types_present) {
  gene_list <- rownames(myRCTD@de_results$gene_fits$mean_val)
  cti_renorm <- get_norm_ref(myRCTD@originalSpatialRNA, myRCTD@cell_type_info$info[[1]],
                             gene_list, myRCTD@internal_vars$proportions)
  return(get_gene_list_type(myRCTD@internal_vars_de$my_beta, myRCTD@internal_vars_de$barcodes, cell_type,
                   myRCTD@spatialRNA@nUMI, gene_list, cti_renorm, cell_types_present,
                   myRCTD@de_results$gene_fits))
}
# Aggregate cell types
#' @export
aggregate_cell_types <- function(myRCTD, barcodes, doublet_mode = T) {
  if(doublet_mode) {
    return(table(myRCTD@results$results_df[barcodes,]$first_type[myRCTD@results$results_df[barcodes,]$spot_class %in% c('singlet','doublet_certain')]) +
             +     table(myRCTD@results$results_df[barcodes,]$second_type[myRCTD@results$results_df[barcodes,]$spot_class %in% c('doublet_certain')]))
  } else {
    return(colSums(myRCTD@results$weights[barcodes,]))
  }
}

get_param_names <- function(X1,X2, cell_types) {
  cnames <- c()
  if(dim(X1)[2] > 0)
    cnames <- unlist(lapply(1:dim(X1)[2], function(x) paste0('1_',x)))
  for(k in 1:length(cell_types))
    cnames <- c(cnames, unlist(lapply(1:dim(X2)[2], function(x) paste0('2_',x,'_',cell_types[k]))))
  return(cnames)
}

get_cell_type_ind <- function(X1,X2, n_cell_types) {
  cnames <- c()
  if(dim(X1)[2] > 0)
    cnames <- rep(0, dim(X1)[2])
  cnames <- c(cnames, unlist(lapply(1:n_cell_types, function(x) rep(x,dim(X2)[2]))))
  return(cnames)
}

choose_cell_types <- function(myRCTD, barcodes, doublet_mode, cell_type_threshold, cell_types) {
  cell_type_count <- aggregate_cell_types(myRCTD, barcodes, doublet_mode = doublet_mode)
  cell_types_default <- names(which(cell_type_count >= cell_type_threshold))
  passed_cell_types = !is.null(cell_types)
  if(passed_cell_types) {
    diff_types <- setdiff(cell_types, cell_types_default)
    if(length(diff_types) > 0)
      stop(paste0('choose_cell_types: cell types: ',paste(diff_types,collapse = ', '),
                  ' detected using aggregate_cell_types to have less than the minimum cell_type_threshold of ',
                  cell_type_threshold,
                  '. To fix this issue, please remove these cell types or reduce the cell_type_threshold'))
    diff_types <- setdiff(cell_types, myRCTD@cell_type_info$info[[2]])
    if(length(diff_types) > 0)
      stop(paste0('choose_cell_types: cell types: ',paste(diff_types,collapse = ', '),
                  ' are not valid cell types in this RCTD object (myRCTD@cell_type_info$info[[2]]). Please check that cell_types only has valid cell types.'))

  } else {
    cell_types <- cell_types_default
  }
  if(length(cell_types) == 0) {
    if(passed_cell_types)
      stop('choose_cell_types: length(cell_types) is 0. Please pass in at least one cell type in the list cell_types')
    else
      stop(paste0('choose_cell_types: length(cell_types) is 0. According to the aggregate_cell_types fn, no cell types occured greater than cell_type_threshold of ',
                  cell_type_threshold, '. Please check that all data is present and consider reducing cell_type_threshold.'))
  }
  if(length(cell_types) == 1) {
    stop('choose_cell_types: length(cell_types) is 1. This is currently not supported. Please consider adding another cell type or contact the developers to have us add in this capability.')
  }
  print(paste0("choose_cell_types: running de with cell types ",paste(cell_types, collapse = ', ')))
  return(cell_types)
}

fdr_sig_genes <- function(gene_list_type, p_val, fdr) {
  N_genes_type <- length(gene_list_type)
  thresh <- (1:N_genes_type)/N_genes_type * fdr
  N_sig <- max(which(p_val[order(p_val)] < thresh))
  if(N_sig >= 1)
    gene_list_sig <- gene_list_type[order(p_val)[1:N_sig]]
  else
    gene_list_sig <- c()
  return(gene_list_sig)
}

get_spline_matrix <- function(puck, df = 15) {
  center_coords <- puck@coords
  center_coords <- sweep(center_coords, 2, apply(center_coords,2, mean), '-')
  center_coords <- center_coords / sd(as.matrix(center_coords))
  sm <- smoothCon(s(x,y,k=df,fx=T,bs='tp'),data=center_coords)[[1]]
  mm <- as.matrix(data.frame(sm$X))
  X2 <- cbind(mm[,(df - 2):df], mm[,1:(df-3)]) #swap intercept, x, and y
  X2[,2:df] <- sweep(X2[,2:df], 2, apply(X2[,2:df],2, mean), '-')
  X2[,2:df] <- sweep(X2[,2:df], 2, apply(X2[,2:df],2, sd), '/') #standardize
  rownames(X2) <- names(puck@nUMI)
  return(X2)
}

check_converged_vec <- function(X1,X2,my_beta, itera, n.iter, error_vec, precision, PRECISION.THRESHOLD) {
  cell_type_ind <- get_cell_type_ind(X1,X2, dim(my_beta)[2])
  converged_vec <- (1:dim(my_beta)[2]) == 0
  if(itera < n.iter) {
    converged_vec <- !converged_vec
    converged_vec[unique(cell_type_ind[precision > PRECISION.THRESHOLD])] <- FALSE
  }
  converged_vec <- converged_vec & (!error_vec)
  names(converged_vec) <- colnames(my_beta)
  return(converged_vec)
}
