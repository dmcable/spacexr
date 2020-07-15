
process_cell_type_info <- function(reference) {
   print("Begin: process_cell_type_info")
   print(paste("init_RCTD: number of cells in reference:", dim(reference@assays$RNA@counts)[2]))
   print(paste("init_RCTD: number of genes in reference:", dim(reference@assays$RNA@counts)[1]))
   cell_counts = table(reference@meta.data$liger_ident_coarse)
   print(cell_counts)
   CELL_MIN = 25 # need at least this for each cell type
   if(min(cell_counts) < CELL_MIN)
      stop(paste0("init_RCTD error: need a minimum of ",CELL_MIN, " cells for each cell type in the reference"))
   cell_type_info <- get_cell_type_info(reference@assays$RNA@counts, reference@meta.data$liger_ident_coarse, reference@meta.data$nUMI)
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
#' @return an \code{\linkS4class{RCTD}} object, which is ready to run the \code{\link{run.RCTD}} function
#' @export
create.RCTD <- function(spatialRNA, reference, max_cores = 8, test_mode = FALSE, gene_cutoff = 0.000125, fc_cutoff = 0.5, gene_cutoff_reg = 0.0002, fc_cutoff_reg = 0.75, UMI_min = 100, UMI_max = 200000,
                         class_df = NULL) {
   config <- list(gene_cutoff = gene_cutoff, fc_cutoff = fc_cutoff, gene_cutoff_reg = gene_cutoff_reg, fc_cutoff_reg = fc_cutoff_reg, UMI_min = UMI_min, max_cores = max_cores,
                 N_epoch = 8, N_X = 50000, K_val = 100, N_fit = 1000, N_epoch_bulk = 30, MIN_CHANGE_BULK = 0.0001, UMI_max = UMI_max, MIN_OBS = 3)
   if(test_mode)
     config <- list(gene_cutoff = .00125, fc_cutoff = 0.5, gene_curoff_reg = 0.002, fc_cutoff_reg = 0.75, UMI_min = 1000,
          N_epoch = 1, N_X = 50000, K_val = 100, N_fit = 50, N_epoch_bulk = 4, MIN_CHANGE_BULK = 1, UMI_max = 200000, MIN_OBS = 3, max_cores = 1)
   cell_type_info <- list(info = process_cell_type_info(reference), renorm = NULL)
   puck = restrict_counts(spatialRNA, rownames(spatialRNA@counts), UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
   print('create.RCTD: getting regression differentially expressed genes: ')
   gene_list_reg = get_de_genes(cell_type_info$info, puck, fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg, MIN_OBS = config$MIN_OBS)
   print('create.RCTD: getting platform effect normalization differentially expressed genes: ')
   gene_list_bulk = get_de_genes(cell_type_info$info, puck, fc_thresh = config$fc_cutoff, expr_thresh = config$gene_cutoff, MIN_OBS = config$MIN_OBS)
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
#' @param RCTD an \code{\linkS4class{RCTD}} object created using the \code{\link{create.RCTD}} function.
#' @param doublet_mode \code{logical} of whether RCTD should be run in doublet_mode.
#' @return an \code{\linkS4class{RCTD}} object containing the results of the RCTD algorithm.
#' @export
run.RCTD <- function(RCTD, doublet_mode = T) {
   RCTD <- fitBulk(RCTD)
   RCTD <- choose_sigma_c(RCTD)
   RCTD <- fitPixels(RCTD, doublet_mode = doublet_mode)
}

