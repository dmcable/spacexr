

process_cell_type_info <- function(reference, cell_type_names, CELL_MIN = 25) {
   print("Begin: process_cell_type_info")
   print(paste("process_cell_type_info: number of cells in reference:", dim(reference@counts)[2]))
   print(paste("process_cell_type_info: number of genes in reference:", dim(reference@counts)[1]))
   cell_counts = table(reference@cell_types)
   print(cell_counts)

   if(min(cell_counts) < CELL_MIN)
      stop(paste0("process_cell_type_info error: need a minimum of ",CELL_MIN, " cells for each cell type in the reference"))
   cell_type_info <- get_cell_type_info(reference@counts, reference@cell_types, reference@nUMI
                                        , cell_type_names = cell_type_names)
   print("End: process_cell_type_info")
   return(cell_type_info)
}

#' Creates an  \code{\linkS4class{RCTD}} object from a scRNA-seq reference \code{Seurat} object and a \code{\linkS4class{SpatialRNA}} object
#'
#' @param spatialRNA a \code{\linkS4class{SpatialRNA}} object to run RCTD on
#' @param reference a \code{\linkS4class{Seurat}} object scRNA-seq reference used for RCTD
#' @param gene_cutoff minimum normalized gene expression for genes to be included in the platform effect normalization step.
#' @param fc_cutoff minimum log-fold-change (across cell types) for genes to be included in the platform effect normalization step.
#' @param gene_cutoff_reg minimum normalized gene expression for genes to be included in the RCTD step.
#' @param fc_cutoff_reg minimum log-fold-change (across cell types) for genes to be included in the RCTD step.
#' @param UMI_min minimum UMI per pixel included in the analysis
#' @param UMI_max maximum UMI per pixel included in the analysis
#' @param max_cores for parallel processing, the number of cores used. If set to 1, parallel processing is not used. The system will additionally be checked for
#' number of available cores.
#' @param class_df (optional) if not NULL, then a dataframe mapping each cell type to a cell class, so that RCTD will report confidence on the class level.
#' @param CELL_MIN_INSTANCE minimum number of cells required per cell type. Default 25, can be lowered if desired.
#' @param cell_type_names A list of cell types to be included from the reference. If NULL, uses all cell types
#' @return an \code{\linkS4class{RCTD}} object, which is ready to run the \code{\link{run.RCTD}} function
#' @export
create.RCTD <- function(spatialRNA, reference, max_cores = 8, test_mode = FALSE, gene_cutoff = 0.000125, fc_cutoff = 0.5, gene_cutoff_reg = 0.0002, fc_cutoff_reg = 0.75, UMI_min = 100, UMI_max = 200000,
                         class_df = NULL, CELL_MIN_INSTANCE = 25, cell_type_names = NULL) {

   config <- list(gene_cutoff = gene_cutoff, fc_cutoff = fc_cutoff, gene_cutoff_reg = gene_cutoff_reg, fc_cutoff_reg = fc_cutoff_reg, UMI_min = UMI_min, max_cores = max_cores,
                 N_epoch = 8, N_X = 50000, K_val = 100, N_fit = 1000, N_epoch_bulk = 30, MIN_CHANGE_BULK = 0.0001, MIN_CHANGE_REG = 0.001, UMI_max = UMI_max, MIN_OBS = 3)
   if(test_mode)
     config <- list(gene_cutoff = .00125, fc_cutoff = 0.5, gene_curoff_reg = 0.002, fc_cutoff_reg = 0.75, UMI_min = 1000,
          N_epoch = 1, N_X = 50000, K_val = 100, N_fit = 50, N_epoch_bulk = 4, MIN_CHANGE_BULK = 1, MIN_CHANGE_REG = 0.001, UMI_max = 200000, MIN_OBS = 3, max_cores = 1)
   if(is.null(cell_type_names))
      cell_type_names <- levels(reference@cell_types)
   cell_type_info <- list(info = process_cell_type_info(reference, CELL_MIN = CELL_MIN_INSTANCE), renorm = NULL)

   puck = restrict_counts(spatialRNA, rownames(spatialRNA@counts), UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
   print('create.RCTD: getting regression differentially expressed genes: ')
   #puckMeans <- rowMeans(sweep(puck@counts, 2 , puck@nUMI, '/'))
   gene_list_reg = get_de_genes(cell_type_info$info, puck, fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg, MIN_OBS = config$MIN_OBS)
   if(length(gene_list_reg) == 0)
      stop("create.RCTD: Error: 0 regression differentially expressed genes found")
   print('create.RCTD: getting platform effect normalization differentially expressed genes: ')
   gene_list_bulk = get_de_genes(cell_type_info$info, puck, fc_thresh = config$fc_cutoff, expr_thresh = config$gene_cutoff, MIN_OBS = config$MIN_OBS)
   if(length(gene_list_bulk) == 0)
      stop("create.RCTD: Error: 0 bulk differentially expressed genes found")
   puck = restrict_counts(puck, gene_list_bulk, UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
   puck = restrict_puck(puck, colnames(puck@counts))
   if(is.null(class_df))
      class_df <- data.frame(cell_type_info$info[[2]], row.names = cell_type_info$info[[2]]); colnames(class_df)[1] = "class"
   internal_vars <- list(gene_list_reg = gene_list_reg, gene_list_bulk = gene_list_bulk, proportions = NULL, class_df = class_df)
   new("RCTD", spatialRNA = puck, reference = reference, config = config, cell_type_info = cell_type_info, internal_vars = internal_vars)
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
#' @param doublet_mode \code{character string}, either "doublet", "multi", or "full" on which mode to run RCTD. Please see above description.
#' @return an \code{\linkS4class{RCTD}} object containing the results of the RCTD algorithm.
#' @export
run.RCTD <- function(RCTD, doublet_mode = "doublet") {
   if(!(doublet_mode %in% c('doublet','multi','full')))
      stop(paste0("run.RCTD: doublet_mode=",doublet_mode, " is not a valid choice. Please set doublet_mode=doublet, multi, or full."))
   RCTD <- fitBulk(RCTD)
   RCTD <- choose_sigma_c(RCTD)
   RCTD <- fitPixels(RCTD, doublet_mode = doublet_mode)
}

