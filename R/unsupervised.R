#' Creates an \code{\linkS4class{RCTD}} object from a scRNA-seq reference \code{Reference} object and a \code{\linkS4class{SpatialRNA}} object
#'
#' @param spatialRNA a \code{\linkS4class{SpatialRNA}} object to run RCTD on
#' @param cell_type_info a \code{cell_type_info} initialization
#' @param gene_list (optional) complete list of genes to be included.
#' @param gene_cutoff_reg minimum normalized gene expression for genes to be included.
#' @param fc_cutoff_reg minimum log-fold-change (across cell types) for genes to be included.
#' @param UMI_min minimum UMI per pixel included in the analysis
#' @param UMI_max maximum UMI per pixel included in the analysis
#' @param UMI_min_sigma minimum UMI per pixel for the \link{choose_sigma_c} function
#' @param max_cores for parallel processing, the number of cores used. If set to 1, parallel processing is not used. The system will additionally be checked for
#' number of available cores.
#' @param class_df (optional) if not NULL, then a dataframe mapping each cell type to a cell class, so that RCTD will report confidence on the class level.
#' @param CONFIDENCE_THRESHOLD (Default 10) the minimum change in likelihood (compared to other cell types) necessary to determine a cell type identity with confidence
#' @param DOUBLET_THRESHOLD (Default 25) the penalty weight of predicting a doublet instead of a singlet for a pixel
#' @return an \code{\linkS4class{RCTD}} object, which is ready to run the \code{\link{run.RCTD}} function
#' @export
create.RCTD.unsupervised <- function(spatialRNA, cell_type_info, max_cores = 4, gene_list = NULL, gene_cutoff_reg = 0.0002, fc_cutoff_reg = 0.75, UMI_min = 100, UMI_max = 20000000, 
                                     UMI_min_sigma = 300, class_df = NULL, CONFIDENCE_THRESHOLD = 10, DOUBLET_THRESHOLD = 25) {

  config <- list(gene_cutoff_reg = gene_cutoff_reg, fc_cutoff_reg = fc_cutoff_reg, UMI_min = UMI_min, UMI_min_sigma = UMI_min_sigma, max_cores = max_cores,
                 N_epoch = 8, N_X = 50000, K_val = 100, N_fit = 1000, N_epoch_bulk = 30,
                 MIN_CHANGE_BULK = 0.0001, MIN_CHANGE_REG = 0.001, MIN_CHANGE_DE = 0.001, UMI_max = UMI_max,
                 MIN_OBS = 3, CONFIDENCE_THRESHOLD = CONFIDENCE_THRESHOLD, DOUBLET_THRESHOLD = DOUBLET_THRESHOLD)
  reference <- new("Reference", cell_types = factor(), counts = as(matrix(), 'dgCMatrix'))
  puck.original = restrict_counts(spatialRNA, rownames(spatialRNA@counts), UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
  message('create.RCTD.unsupervised: getting regression differentially expressed genes: ')
  if(is.null(gene_list))
    gene_list = get_de_genes(cell_type_info$info, puck.original, fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg, MIN_OBS = config$MIN_OBS)
  if(length(gene_list) == 0)
    stop("create.RCTD.unsupervised: Error: 0 regression differentially expressed genes found")
  puck = restrict_counts(puck.original, gene_list, UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
  puck = restrict_puck(puck, colnames(puck@counts))
  if(is.null(class_df))
    class_df <- data.frame(cell_type_info$info[[2]], row.names = cell_type_info$info[[2]]); colnames(class_df)[1] = "class"
  cell_type_info$info[[1]] <- cell_type_info$info[[1]][gene_list, ]
  cell_type_info$renorm[[1]] <- cell_type_info$renorm[[1]][gene_list, ]
  internal_vars <- list(gene_list_reg = gene_list, proportions = NULL, class_df = class_df, cell_types_assigned = F)
  new("RCTD", spatialRNA = puck, originalSpatialRNA = puck.original, reference = reference, config = config, cell_type_info = cell_type_info, internal_vars = internal_vars)
}


#' Runs the unsupervised pipeline on a \code{\linkS4class{RCTD}} object
#'
#' Equivalent to sequentially running the functions \code{\link{choose_sigma_c}} and \code{\link{iter.optim}}
#'
#' If in doublet mode, fits at most two cell types per pixel. It classifies each pixel as 'singlet' or 'doublet' and searches for the cell types
#' on the pixel. If in full mode, can fit any number of cell types on each pixel. If in subtype mode, fits only pixels belonging to given cell type(s).
#'
#' @param RCTD an \code{\linkS4class{RCTD}} object created using the \code{\link{create.RCTD}} function.
#' @param doublet_mode \code{character string}, either "doublet", "subtype", or "full" on which mode to run iter.optim. Please see above description.
#' @param cell_types the cell types used for CSIDE. If null, all cell types in \code{cell_type_info} will be chosen.
#' @param n_iter maximum number of optimization iterations
#' @param MIN_CHANGE minimum change required to terminate optimization
#' @return an \code{\linkS4class{RCTD}} object containing the results of the unsupervised algorithm. Please see \code{\linkS4class{RCTD}}
#' documentation for more information on interpreting the content of the RCTD object.
#' @export
run.unsupervised <- function(RCTD, doublet_mode = 'doublet', cell_types = NULL, n_iter = 50, MIN_CHANGE = 0.001) {
  if (!(doublet_mode %in% c('doublet', 'full', 'subtype')))
    stop(paste0("run.unsupervised: doublet_mode=", doublet_mode, " is not a valid choice. Please set doublet_mode=doublet, full, or subtype."))
  if (doublet_mode == 'subtype')
    cell_types = RCTD@internal_vars$subtypes
  RCTD@config$RCTDmode <- doublet_mode
  RCTD <- choose_sigma_c(RCTD)
  RCTD <- iter.optim(RCTD, doublet_mode = doublet_mode, cell_types = cell_types, n_iter = n_iter, MIN_CHANGE = MIN_CHANGE)
}


#' Runs the unsupervised iterative optimization algorithm
#'
#' If in doublet mode, fits at most two cell types per pixel. It classifies each pixel as 'singlet' or 'doublet' and searches for the cell types
#' on the pixel. If in full mode, can fit any number of cell types on each pixel. If in subtype mode, fits only pixels belonging to given cell type(s).
#'
#' @param RCTD an \code{\linkS4class{RCTD}} object after running the \code{\link{choose_sigma_c}} function.
#' @param doublet_mode \code{character string}, either "doublet", "subtype", or "full" on which mode to run iter.optim. Please see above description.
#' @param cell_types the cell types used for CSIDE. If null, all cell types in \code{cell_type_info} will be chosen.
#' @param n_iter maximum number of optimization iterations
#' @param MIN_CHANGE minimum change required to terminate optimization
#' @return an \code{\linkS4class{RCTD}} object containing the results of the unsupervised algorithm.
#' @export
iter.optim <- function(RCTD, doublet_mode = 'doublet', cell_types = NULL, n_iter = 50, MIN_CHANGE = 0.001) {
  if (is.null(cell_types)) 
    cell_types <- RCTD@cell_type_info$info[[2]]
  barcodes <- intersect(names(RCTD@spatialRNA@nUMI), colnames(RCTD@spatialRNA@counts))
  X <- as.matrix(rep(1, length(barcodes)))
  rownames(X) <- barcodes
  RCTD@de_results$gene_fits$mean_val <- as.matrix(log(RCTD@cell_type_info$info[[1]]))
  originalSpatialRNA <- RCTD@originalSpatialRNA
  DOUBLET_THRESHOLD <- RCTD@config$DOUBLET_THRESHOLD
  RCTD@originalSpatialRNA <- RCTD@spatialRNA
  RCTD@config$DOUBLET_THRESHOLD <- 10
  RCTD@config$MIN_CHANGE_REG <- 1e-2
  RCTD@config$MIN_CHANGE_DE <- 1e-2

  message('iter.optim: assigning initial cell types')
  if (doublet_mode == 'subtype') {
    initialSol <- list(weights = RCTD@results$weights)
  } else {
    initialSol <- NULL
  }
  RCTD <- fitPixels(RCTD, doublet_mode = doublet_mode, initialSol = initialSol)
  RCTD_list <- list(RCTD)
  for (i in 1:n_iter) {
    message(paste('iter.optim: running iteration', i))
    message('fitting gene profiles')
    initialSol <- RCTD@de_results$gene_fits$mean_val
    RCTD <- run.CSIDE(RCTD, X, barcodes, cell_types, doublet_mode = (doublet_mode == 'doublet'), cell_type_threshold = 0, gene_threshold = -1,
                      sigma_gene = F, test_genes_sig = F, params_to_test = 1, logs = T, initialSol = initialSol[, cell_types])
    info <- as.data.frame(exp(RCTD@de_results$gene_fits$mean_val))
    RCTD@cell_type_info$info[[1]][, cell_types] <- info 
    RCTD@cell_type_info$renorm[[1]][, cell_types] <- info 

    message('fitting cell types')
    initialSol <- list(weights = RCTD@results$weights, doublet_mat = RCTD@results$doublet_weights)
    RCTD <- fitPixels(RCTD, doublet_mode = doublet_mode, initialSol = initialSol)

    change <- weights_change(RCTD_list[[i]], RCTD)
    message(paste('change:', change))
    RCTD@config$MIN_CHANGE_REG <- max(min(1e-2, change ** 2), 1e-3)
    RCTD@config$MIN_CHANGE_DE <- max(min(1e-2, change ** 2), 1e-3)
    if (change < MIN_CHANGE) break
    RCTD_list[[i+1]] <- RCTD
  }
  if (doublet_mode == 'doublet') {
    results <- RCTD@results$results_df
    reassign <- rownames(results[results$spot_class == 'doublet_certain' & results$singlet_score - results$min_score < DOUBLET_THRESHOLD, ])
    if (length(reassign) > 0)
      RCTD@results$results_df[reassign, ]$spot_class <- 'singlet'
  }
  RCTD@originalSpatialRNA <- originalSpatialRNA
  RCTD_list[[i+1]] <- RCTD
  RCTD_list
}
