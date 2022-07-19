#' Calculates weight matrix from doublet mode results
#'
#' Used during \code{\link{initialize.subtypes}} when unsupervised was run on doublet mode
#'
#' @param RCTD a \code{\linkS4class{RCTD}} object
#' @return a matrix of doublet mode weights
#' @export
weights_from_results <- function(RCTD) {
  df <- RCTD@results$results_df
  wd <- RCTD@results$weights_doublet
  singlets <- rownames(df[df$spot_class == 'singlet', ])
  doublets <- rownames(df[df$spot_class == 'doublet_certain', ])
  barcodes <- c(singlets, doublets)
  cell_types <- levels(df$first_type)
  weights <- matrix(0, nrow = length(barcodes), ncol = length(cell_types))
  rownames(weights) <- barcodes; colnames(weights) <- cell_types
  for (b in singlets) {
    weights[b, as.character(df[b, 'first_type'])] <- 1
  }
  for (b in doublets) {
    weights[b, as.character(df[b, 'first_type'])] <- wd[b, 'first_type']
    weights[b, as.character(df[b, 'second_type'])] <- wd[b, 'second_type']
  }
  return(weights)
}

#' Calculates the weight change between two \code{\linkS4class{RCTD}} objects
#'
#' Used during \code{\link{iter.optim}} as the termination criterion
#'
#' @param RCTD1 a \code{\linkS4class{RCTD}} object
#' @param RCTD2 a \code{\linkS4class{RCTD}} object
#' @return \code{cell_type_info}, which is ready to run the \code{\link{create.RCTD.unsupervised}} function
#' @export
weights_change <- function(RCTD1, RCTD2) {
  weights1 <- RCTD1@results$weights
  weights2 <- RCTD2@results$weights
  return(norm(weights1 - weights2) / dim(weights1)[1])
}

#' Generates \code{cell_type_info} from cluster assignments on a \code{\linkS4class{SpatialRNA}} object
#'
#' @param puck a \code{\linkS4class{SpatialRNA}} object
#' @param clusters a \code{data.frame} with barcodes as rownames and cluster assignments under a \code{cell_types} column
#' @return \code{cell_type_info}, which is ready to run the \code{\link{create.RCTD.unsupervised}} function
#' @export
cell_type_info_from_clusters <- function(puck, clusters) {
  counts <- puck@counts[, rownames(clusters)]
  nUMI <- puck@nUMI[rownames(clusters)]
  cell_types <- clusters$cell_types
  names(cell_types) <- rownames(clusters)
  cell_type_info <- get_cell_type_info(counts, cell_types, nUMI)
  return(list(info = cell_type_info, renorm = cell_type_info))
}

#' Generates \code{cell_type_info} from CSIDE on a \code{\linkS4class{RCTD}} object
#'
#' @param RCTD a \code{\linkS4class{RCTD}} object to run CSIDE
#' @param cell_types the cell types used for CSIDE. If null, all cell types in \code{cell_type_info} will be chosen.
#' @param gene_list the genes to compute expression for. If null, all genes will be included.
#' @return \code{cell_type_info}, which is ready to run the \code{\link{create.RCTD.unsupervised}} function
#' @export
cell_type_info_from_de <- function(RCTD,
                                   cell_types = NULL,
                                   gene_list = NULL) {
  if (is.null(cell_types))
    cell_types <- RCTD@cell_type_info$info[[2]]
  barcodes <- intersect(names(RCTD@spatialRNA@nUMI), colnames(RCTD@spatialRNA@counts))
  X <- as.matrix(rep(1, length(barcodes)))
  rownames(X) <- barcodes
  RCTD <- run.CSIDE(RCTD, X, barcodes, cell_types, gene_list_tot = gene_list, doublet_mode = (RCTD@config$RCTDmode == "doublet"), 
                    cell_type_threshold = 0, gene_threshold = -1, sigma_gene = F, test_genes_sig = F, params_to_test = 1, logs = T)
  cell_type_info <- list(as.data.frame(exp(RCTD@de_results$gene_fits$mean_val)), cell_types, length(cell_types))
  return(list(info = cell_type_info, renorm = cell_type_info))
}
