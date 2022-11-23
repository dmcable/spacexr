

process_cell_type_info <- function(reference, cell_type_names, CELL_MIN = 25) {
   message("Begin: process_cell_type_info")
   message(paste("process_cell_type_info: number of cells in reference:", dim(reference@counts)[2]))
   message(paste("process_cell_type_info: number of genes in reference:", dim(reference@counts)[1]))
   cell_counts = table(reference@cell_types)
   print(cell_counts)

   if(min(cell_counts) < CELL_MIN)
      stop(paste0("process_cell_type_info error: need a minimum of ",CELL_MIN, " cells for each cell type in the reference"))
   cell_type_info <- get_cell_type_info(reference@counts, reference@cell_types, reference@nUMI
                                        , cell_type_names = cell_type_names)
   message("End: process_cell_type_info")
   return(cell_type_info)
}

#' Creates an \code{\linkS4class{RCTD}} object from a scRNA-seq reference \code{Reference} object and a \code{\linkS4class{SpatialRNA}} object
#'
#' @param spatialRNA a \code{\linkS4class{SpatialRNA}} object to run RCTD on
#' @param reference a \code{\linkS4class{Reference}} object scRNA-seq reference used for RCTD
#' @param gene_cutoff minimum normalized gene expression for genes to be included in the platform effect normalization step.
#' @param fc_cutoff minimum log-fold-change (across cell types) for genes to be included in the platform effect normalization step.
#' @param gene_cutoff_reg minimum normalized gene expression for genes to be included in the RCTD step.
#' @param fc_cutoff_reg minimum log-fold-change (across cell types) for genes to be included in the RCTD step.
#' @param UMI_min minimum UMI per pixel included in the analysis
#' @param UMI_max maximum UMI per pixel included in the analysis
#' @param counts_MIN (default 10) minimum total counts per pixel of genes used in the analysis.
#' @param UMI_min_sigma minimum UMI per pixel for the \link{choose_sigma_c} function
#' @param max_cores for parallel processing, the number of cores used. If set to 1, parallel processing is not used. The system will additionally be checked for
#' number of available cores.
#' @param class_df (optional) if not NULL, then a dataframe mapping each cell type to a cell class, so that RCTD will report confidence on the class level.
#' @param CELL_MIN_INSTANCE minimum number of cells required per cell type. Default 25, can be lowered if desired.
#' @param cell_type_names A list of cell types to be included from the reference. If NULL, uses all cell types
#' @param MAX_MULTI_TYPES (multi-mode only) Default 4, max number of cell types per pixel
#' @param cell_type_profiles Default NULL, option to pass in cell type profiles in directly as a genes by cell type matrix, including gene names and cell type names.
#' If this option is used, reference will be ignored.
#' @param keep_reference (Default FALSE) if true, keeps the \code{reference} object stored within the \code{\linkS4class{RCTD}} object
#' @param CONFIDENCE_THRESHOLD (Default 10) the minimum change in likelihood (compared to other cell types) necessary to determine a cell type identity with confidence
#' @param DOUBLET_THRESHOLD (Default 25) the penalty weight of predicting a doublet instead of a singlet for a pixel
#' @return an \code{\linkS4class{RCTD}} object, which is ready to run the \code{\link{run.RCTD}} function
#' @export
create.RCTD <- function(spatialRNA, reference, max_cores = 4, test_mode = FALSE, gene_cutoff = 0.000125, fc_cutoff = 0.5, gene_cutoff_reg = 0.0002, fc_cutoff_reg = 0.75, UMI_min = 100, UMI_max = 20000000, counts_MIN = 10, UMI_min_sigma = 300,
                         class_df = NULL, CELL_MIN_INSTANCE = 25, cell_type_names = NULL, MAX_MULTI_TYPES = 4, keep_reference = F, cell_type_profiles = NULL, CONFIDENCE_THRESHOLD = 10, DOUBLET_THRESHOLD = 25) {

   config <- list(gene_cutoff = gene_cutoff, fc_cutoff = fc_cutoff, gene_cutoff_reg = gene_cutoff_reg, fc_cutoff_reg = fc_cutoff_reg, UMI_min = UMI_min, UMI_min_sigma = UMI_min_sigma, max_cores = max_cores,
                 N_epoch = 8, N_X = 50000, K_val = 100, N_fit = 1000, N_epoch_bulk = 30,
                 MIN_CHANGE_BULK = 0.0001, MIN_CHANGE_REG = 0.001, UMI_max = UMI_max, counts_MIN = counts_MIN,
                 MIN_OBS = 3, MAX_MULTI_TYPES = MAX_MULTI_TYPES, CONFIDENCE_THRESHOLD = CONFIDENCE_THRESHOLD, DOUBLET_THRESHOLD = DOUBLET_THRESHOLD)
   if(test_mode)
     config <- list(gene_cutoff = .00125, fc_cutoff = 0.5, gene_cutoff_reg = 0.002, fc_cutoff_reg = 0.75, UMI_min = 1000,
          N_epoch = 1, N_X = 50000, K_val = 100, N_fit = 50, N_epoch_bulk = 4, MIN_CHANGE_BULK = 1,
          MIN_CHANGE_REG = 0.001, UMI_max = 200000, MIN_OBS = 3, max_cores = 1, counts_MIN = 5,
          UMI_min_sigma = 300, MAX_MULTI_TYPES = MAX_MULTI_TYPES, CONFIDENCE_THRESHOLD = CONFIDENCE_THRESHOLD, DOUBLET_THRESHOLD = DOUBLET_THRESHOLD)
   if(is.null(cell_type_profiles)) {
      if(is.null(cell_type_names))
         cell_type_names <- levels(reference@cell_types)
      cell_type_info <- list(info = process_cell_type_info(reference, cell_type_names = cell_type_names, CELL_MIN = CELL_MIN_INSTANCE), renorm = NULL)
   } else {
      cell_type_names <- colnames(cell_type_profiles)
      cell_type_info <- list(info = list(cell_type_profiles, cell_type_names, length(cell_type_names)),
                             renorm = NULL)
   }
   if(!keep_reference)
      reference <- create_downsampled_data(reference, n_samples = 5)
   puck.original = restrict_counts(spatialRNA, rownames(spatialRNA@counts),
                                   UMI_thresh = config$UMI_min, UMI_max = config$UMI_max,
                                   counts_thresh = config$counts_MIN)
   message('create.RCTD: getting regression differentially expressed genes: ')
   #puckMeans <- rowMeans(sweep(puck@counts, 2 , puck@nUMI, '/'))
   gene_list_reg = get_de_genes(cell_type_info$info, puck.original, fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg, MIN_OBS = config$MIN_OBS)
   if(length(gene_list_reg) < 10)
      stop("create.RCTD: Error: fewer than 10 regression differentially expressed genes found")
   message('create.RCTD: getting platform effect normalization differentially expressed genes: ')
   gene_list_bulk = get_de_genes(cell_type_info$info, puck.original, fc_thresh = config$fc_cutoff, expr_thresh = config$gene_cutoff, MIN_OBS = config$MIN_OBS)
   if(length(gene_list_bulk) < 10)
      stop("create.RCTD: Error: fewer than 10 bulk differentially expressed genes found")
   puck = restrict_counts(puck.original, gene_list_bulk, UMI_thresh = config$UMI_min,
                          UMI_max = config$UMI_max, counts_thresh = config$counts_MIN)
   puck = restrict_puck(puck, colnames(puck@counts))
   if(is.null(class_df))
      class_df <- data.frame(cell_type_info$info[[2]], row.names = cell_type_info$info[[2]]); colnames(class_df)[1] = "class"
   internal_vars <- list(gene_list_reg = gene_list_reg, gene_list_bulk = gene_list_bulk, proportions = NULL, class_df = class_df, cell_types_assigned = F)
   new("RCTD", spatialRNA = puck, originalSpatialRNA = puck.original, reference = reference, config = config, cell_type_info = cell_type_info, internal_vars = internal_vars)
}

#' Runs the RCTD pipeline on a \code{\linkS4class{RCTD}} object
#'
#' Equivalent to sequentially running the functions \code{\link{fitBulk}}, \code{\link{choose_sigma_c}}, and \code{\link{fitPixels}}
#'
#' If in doublet mode, fits at most two cell types per pixel. It classifies each pixel as 'singlet' or 'doublet' and searches for the cell types
#' on the pixel. If in full mode, can fit any number of cell types on each pixel. In multi mode, cell types are added using a greedy algorithm,
#' up to a fixed number.
#'
#' @param RCTD an \code{\linkS4class{RCTD}} object created using the \code{\link{create.RCTD}} function.
#' @return an \code{\linkS4class{RCTD}} object containing the results of the RCTD algorithm. Please see \code{\linkS4class{RCTD}}
#' documentation for more information on interpreting the content of the RCTD object.
#' @export
run.RCTD <- function(RCTD, doublet_mode = "doublet") {
   if(!(doublet_mode %in% c('doublet','multi','full')))
      stop(paste0("run.RCTD: doublet_mode=",doublet_mode, " is not a valid choice. Please set doublet_mode=doublet, multi, or full."))
   RCTD@config$RCTDmode <- doublet_mode
   RCTD <- fitBulk(RCTD)
   RCTD <- choose_sigma_c(RCTD)
   RCTD <- fitPixels(RCTD, doublet_mode = doublet_mode)
}

#exports RCTD results as csv files
export.RCTD <- function(RCTD, datadir) {
   doublet_mode <- myRCTD@config$doublet_mode
   if(is.null(doublet_mode))
      stop("RCTD@config$doublet_mode is NULL. Please set to one of 'doublet', 'multi', 'full'.")
   if(!(doublet_mode %in% c('doublet','multi','full')))
      stop(paste0("export.RCTD: RCTD@config$doublet_mode=",doublet_mode, " is not a valid choice. Please set doublet_mode=doublet, multi, or full."))
   if(doublet_mode == 'multi')
      stop("export.RCTD not implemented for RCTD@config$doublet_mode = 'multi'. Please contact the developers for assistance.")
   write.csv(RCTD@results$weights, file.path(datadir, 'weights.csv'))
   if(doublet_mode == 'doublet') {
      write.csv(RCTD@results$weights_doublet, file.path(datadir, 'weights_doublet.csv'))
      write.csv(RCTD@results$results_df, file.path(datadir, 'results_df.csv'))
   }
}

check_vector <- function(variable, var_name, f_name, require_int = FALSE) {
   if(!is.atomic(variable))
      stop(paste0(f_name,': ',var_name,' is not an atomic vector. Please format ',var_name,' as an atomic vector.'))
   if(!is.numeric(variable))
      stop(paste0(f_name,': ',var_name,' is not numeric'))
   if(is.null(names(variable)))
      stop(paste0(f_name,': names(',var_name,') is null. Please enter names'))
   if(length(variable) == 1)
      stop(paste0(f_name,': the length of ',var_name,' is 1, indicating only one element present. Please format ',var_name,' so that
         the length is greater than 1.'))
   if(require_int) {
      if(max(abs(variable %% 1)) > 1e-6)
         stop(paste0(f_name,': variable does not contain integers'))
   }
}


#' Updates an old \code{\linkS4class{RCTD}} object to be compatible with the current
#' version of \code{spacexr}.
#'
#' @param RCTD an \code{\linkS4class{RCTD}} object (potentially from an older version.
#' @return an \code{\linkS4class{RCTD}} object updated to be compatible with the current version
#' of \code{spacexr}.
#' @export
convert.old.RCTD <- function(myRCTD) {
   if(class(myRCTD@reference )[1] == 'Seurat') {
      ref <- convert_old_reference(myRCTD@reference)
   } else {
      ref <- myRCTD@reference
   }
   if(attr(class(myRCTD@spatialRNA),'package') != 'spacexr')
      myRCTD@spatialRNA <- coerce_old(myRCTD@spatialRNA)
   if(attr(class(myRCTD@originalSpatialRNA),'package') != 'spacexr')
      myRCTD@originalSpatialRNA <- coerce_old(myRCTD@originalSpatialRNA)
   if(attr(class(ref),'package') != 'spacexr')
      ref <- coerce_deglam_reference(ref)
   new("RCTD", spatialRNA = myRCTD@spatialRNA, originalSpatialRNA = myRCTD@spatialRNA, reference = ref,
       config = myRCTD@config, cell_type_info = myRCTD@cell_type_info, internal_vars = myRCTD@internal_vars, results = myRCTD@results)
}
