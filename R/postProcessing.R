#functions for processing data after RCTD is fit to the data

# Collects RCTD results
gather_results <- function(RCTD, results) {

  message('Step 4/4: Gather Results')
  pb <- txtProgressBar(max = 3, style = 3)

  cell_type_names <- RCTD@cell_type_info$renorm[[2]]
  barcodes <- colnames(RCTD@spatialRNA@counts)
  N <- length(results)

  empty_cell_types <- factor(character(N),levels = cell_type_names)
  spot_levels <- c("reject", "singlet", "doublet_certain", "doublet_uncertain")

  setTxtProgressBar(pb, 1)
  
  results_df <- data.frame(spot_class = factor(sapply(results,function(X){return(X$spot_class)}),levels=spot_levels),
                           first_type = sapply(results,function(X){return(X$first_type)}),
                           scond_type = sapply(results,function(X){return(X$second_type)}),
                           first_class = sapply(results,function(X){return(X$first_class)}),
                           second_class = sapply(results,function(X){return(X$second_class)}),
                           min_score = sapply(results,function(X){return(X$min_score)}),
                           singlet_score = sapply(results,function(X){return(X$singlet_score)}),
                           conv_all = sapply(results,function(X){return(X$conv_all)}),
                           conv_doublet = sapply(results,function(X){return(X$conv_doublet)}))

  setTxtProgressBar(pb, 2)

  weights_doublet <- do.call(rbind,lapply(results,function(X){return(X$doublet_weights)}))
  weights <- do.call(rbind,lapply(results,function(X){return(X$all_weights)}))

  rownames(weights) <- barcodes
  rownames(weights_doublet) <- barcodes
  colnames(weights) <- cell_type_names
  colnames(weights_doublet) <- c('first_type', 'second_type')

  score_mat <- lapply(results,function(X){return(X$score_mat)})
  singlet_scores <- lapply(results,function(X){return(X$singlet_scores)})

  setTxtProgressBar(pb, 3)

  rownames(results_df) = barcodes
  RCTD@results <- list(results_df = results_df, weights = weights, weights_doublet = weights_doublet,
                       score_mat = score_mat, singlet_scores = singlet_scores)
  return(RCTD)
}

get_decomposed_data_full_doublet <- function(gene_list, puck, weights, ct_info) {
  first_DGE <- Matrix(0, nrow = dim(weights)[1], ncol = length(gene_list))
  second_DGE <- Matrix(0, nrow = dim(weights)[1], ncol = length(gene_list))
  rownames(first_DGE) = rownames(weights); rownames(second_DGE) = rownames(weights)
  colnames(first_DGE) = gene_list; colnames(second_DGE) = gene_list
  for(ind in 1:dim(weights)[1]) {
    barcode = rownames(weights)[ind]
    doub_res <- decompose_doublet_fast(puck@counts[gene_list,barcode], weights[barcode,], gene_list, ct_info, colnames(weights)[1],colnames(weights)[2])
    first_DGE[barcode,] <- doub_res$expect_1; second_DGE[barcode,] <- doub_res$expect_2
  }
  norm1 <- sweep(first_DGE, 1, weights[rownames(weights),1] * puck@nUMI[rownames(first_DGE)], '/')
  norm2 <- sweep(second_DGE, 1, weights[rownames(weights),2] * puck@nUMI[rownames(second_DGE)], '/')
  all_DGE <- rbind(norm1, norm2)
  cell_type_labels <- unlist(list(rep(colnames(weights)[1], dim(weights)[1]), rep(colnames(weights)[2], dim(weights)[1])))
  nUMI <- c(weights[rownames(weights),1] *puck@nUMI[rownames(first_DGE)], weights[rownames(weights),2]*puck@nUMI[rownames(second_DGE)])
  rownames(all_DGE) = 1:dim(all_DGE)[1]; names(cell_type_labels) <- 1:dim(all_DGE)[1]; names(nUMI) <- 1:dim(all_DGE)[1]
  ref_d <- Reference(t(all_DGE), factor(cell_type_labels), nUMI, require_int = F)
  return(ref_d)
}

#' Decomposes SpatialRNA data into individual cells
#'
#' Warning: in the current RCTD version, this function is deprecated, and is no longer supported.
#' For differential expression tasks, we instead recommend the RCTDE method.
#'
#' Applied to the output of \code{\link{gather_results}}.
#' Singlet pixels are left unchanged, and doublet_certain conditions are
#' decomposed into single cells.
#'
#' @param results_df a dataframe of RCTD results (see \code{\link{gather_results}})
#' @param puck an object of type \linkS4class{SpatialRNA}
#' @param weights_doublet a dataframe of predicted weights in doublet mode
#' @param results_df a dataframe of RCTD results
#' @param gene_list a list of genes to be used for the decomposition
#' @param cell_type_info cell type information and profiles of each cell, calculated from the scRNA-seq
#' reference (see \code{\link{get_cell_type_info}})
#' @return An object of type \linkS4class{SpatialRNA} representing the decomposed cells
get_decomposed_data <- function(results_df, gene_list, puck, weights_doublet, cell_type_info) {
  warning("In the current RCTD version, get_decomposed_data is deprecated, and is no longer supported. For differential expression tasks, we instead recommend the RCTDE method.")
  doublets <- results_df[results_df$spot_class == "doublet_certain",]
  first_DGE <- Matrix(0, nrow = dim(doublets)[1], ncol = length(gene_list))
  second_DGE <- Matrix(0, nrow = dim(doublets)[1], ncol = length(gene_list))
  rownames(first_DGE) = rownames(doublets); rownames(second_DGE) = rownames(doublets)
  colnames(first_DGE) = gene_list; colnames(second_DGE) = gene_list
  for(ind in 1:dim(doublets)[1]) {
    print(ind)
    barcode = rownames(doublets)[ind]
    doub_res <- decompose_doublet_fast(puck@counts[gene_list,barcode], weights_doublet[barcode,], gene_list, cell_type_info, results_df[barcode,"first_type"],results_df[barcode,"second_type"])
    first_DGE[barcode,] <- doub_res$expect_1; second_DGE[barcode,] <- doub_res$expect_2
  }
  singlet_id <- results_df$spot_class == "singlet"
  all_DGE <- rbind(first_DGE, second_DGE, t(puck@counts[gene_list, singlet_id]))
  cell_type_labels <- unlist(list(doublets$first_type, doublets$second_type, results_df[singlet_id, "first_type"]))
  coords <- rbind(puck@coords[rownames(doublets),c('x','y')], puck@coords[rownames(doublets),c('x','y')], puck@coords[singlet_id,c('x','y')])
  nUMI <- c(weights_doublet[rownames(doublets),"first_type"] *puck@nUMI[rownames(first_DGE)], weights_doublet[rownames(doublets),"second_type"]*puck@nUMI[rownames(second_DGE)], puck@nUMI[singlet_id])
  rownames(coords) = 1:dim(coords)[1]; names(nUMI) = 1:dim(coords)[1]
  rownames(all_DGE) = 1:dim(coords)[1]
  puck_d <- SpatialRNA(coords, t(all_DGE), nUMI, require_int = F)
  return(puck_d)
}

#' Assigns a cell type `weights` matrix to an \code{\linkS4class{RCTD}} object
#'
#' @param myRCTD a \code{\linkS4class{RCTD}} object to be assigned weights.
#' @param weights a matrix of weights (pixels by cell types). weights must be normalized
#' to have rows sum to 1. Furthermore, rownames and colnames must be assigned as
#' pixel names and cell types respectively.
#' @return the \code{\linkS4class{RCTD}} object with weights assigned.
#' @export
import_weights <- function(myRCTD, weights) {
  if(is.null(rownames(weights)) || is.null(colnames(weights)))
    stop('import_weights: weights must contain rownames and colnames. rownames and colnames must be assigned as pixel names and cell types respectively.')
  eps <- 1e-6
  if(any(abs(rowSums(weights) - 1) > eps))
    stop('import_weights: weights must be normalized to have rows sum to 1.')
  myRCTD@results$weights <- weights
  myRCTD@internal_vars$proportions <- colMeans(weights)/sum(colMeans(weights))
  myRCTD@internal_vars$cell_types_assigned <- TRUE
  myRCTD@internal_vars$sigma <- 1
  set_global_Q_all()
  myRCTD@internal_vars$Q_mat <- Q_mat_all[['100']]
  myRCTD@internal_vars$X_vals <- X_vals
  return(myRCTD)
}

#' Normalizes the `weights` matrix from the RCTD results object
#'
#' @param weights a matrix of weights to be normalized
#' @return norm.weights a normalized matrix of weights where rows sum to one.
#' @export
normalize_weights <- function(weights) {
  return(sweep(weights, 1, rowSums(weights), '/'))
}
