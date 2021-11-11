#' RCTD: an R package for assigning cell types to spatial transcriptomics data.
#'
#' @section Running RCTD:
#'
#' To get started, create a \code{\linkS4class{SpatialRNA}} object (called \code{puck} here) for the
#' spatial transcriptomics data and a \code{\linkS4class{Reference}} object (called \code{reference} here) for the
#' scRNA-seq data. Then simply run RCTD as:
#'
#' \code{myRCTD <- create.RCTD(puck, reference)}
#'
#' \code{myRCTD <- run.RCTD(myRCTD)}
#'
#' @docType package
#' @name RCTD
NULL

#' An S4 class to represent Spatial Transcriptomic data
#'
#' @slot coords a dataframe with x and y coordinates of each pixel
#' @slot counts a sparse matrix of raw counts for each gene (rowname)
#' and each pixel (colnames or barcodes)
#' @slot n_cell_type the number of cell types
#' @slot cell_type_names a list of cell type names
#' @slot nUMI a named list (by barcode) of total UMIs per pixel
#' @slot cell_labels a factor of cell type labels for each pixel
#' @export
#' @import Matrix
#' @import doParallel
#' @import foreach
#' @importClassesFrom Matrix Matrix dgCMatrix
setClass("SpatialRNA",
         slots = c(
           coords = "data.frame",
           counts = "dgCMatrix",
           nUMI = "numeric"
         ),
         prototype = list(
           coords = data.frame(NULL),
           counts = NULL,
           nUMI = NA_integer_
         )
)

#' An S4 class to represent Single-Cell RNA-seq reference
#'
#' @slot cell_types a factor of cell type identities for each cell
#' @slot counts a sparse matrix of raw counts for each gene (rowname)
#' and each cell (colnames or barcodes)
#' @slot nUMI an atomic vector of numeric UMI counts per cell
#' @importClassesFrom Matrix Matrix dgCMatrix
#' @export
setClass("Reference",
  slots = c(
   cell_types = "factor",
   counts = "dgCMatrix",
   nUMI = "numeric"
  ),
  prototype = list(
   cell_types = NULL,
   counts = NULL,
   nUMI = NA_integer_
  )
)

#' An S4 class used to run the RCTD algorithm
#'
#' Created using the \code{\link{create.RCTD}} function, a user can run RCTD using the \code{\link{run.RCTD}} function.
#'
#' @slot spatialRNA a \code{\linkS4class{SpatialRNA}} object containing the Spatial RNA dataset to be used for RCTD
#' @slot originalSpatialRNA a \code{\linkS4class{SpatialRNA}} object containing the Spatial RNA dataset with all genes
#' @slot reference a \code{\linkS4class{Reference}} object containing the cell type-labeled single cell reference
#' @slot config a list of configuration options, set using the \code{\link{create.RCTD}} function
#' @slot cell_type_info a named list of cell type profiles (means), containing two elements: \code{info}, directly calculated from the scRNA-seq reference, and
#' \code{renorm}, renormalized the match the SpatialRNA dataset.
#' @slot internal_vars a list of internal variables used by RCTD's computation
#' @slot results (created after running RCTD) a list of results_df (a dataframe of RCTD results in doublet mode),
#' weights (a dataframe of RCTD predicted weights in full mode), and weights_doublet (a
#' dataframe of predicted weights in doublet mode, with cell type information in results_df).
#'
#' In doublet-mode, The results of 'doublet_mode' are stored in `@results$results_df` and `@results$weights_doublet`, the weights of each cell type.
#' More specifically, the `results_df` object contains one column per pixel (barcodes as rownames). Important columns are:
#' * `spot_class`, a factor variable representing RCTD's classification in doublet mode: "singlet" (1 cell type on pixel), "doublet_certain" (2 cell types on pixel), "doublet_uncertain" (2 cell types on pixel, but only confident of 1), "reject" (no prediction given for pixel).
#' * Next, the `first_type` column gives the first cell type predicted on the bead (for all spot_class conditions except "reject").
#' * The `second_type column` gives the second cell type predicted on the bead for doublet spot_class conditions (not a confident prediction for "doublet_uncertain").
#'
#' Note that in multi-mode, results consists of a list of results for each pixel, which contains all_weights (weights from full mode),
#' cell_type_list (cell types on multi mode), conf_list (which cell types are confident on multi mode) and
#' sub_weights (proportions of cell types on multi mode).
#' @slot de_results results of the RCTDE algorithm. Contains `gene_fits`, which contains the results of fits on individual genes,
#' whereas `res_gene_list` is a list, for each cell type, of significant genes detected by RCTDE.
#' @slot internal_vars_de a list of variables that are used internally by RCTDE
#' @export
#' @import Matrix
#' @import doParallel
#' @import foreach
#' @importClassesFrom Matrix Matrix dgCMatrix
setClass("RCTD",
  slots = c(
   spatialRNA = "SpatialRNA",
   originalSpatialRNA = "SpatialRNA",
   reference = "Reference",
   config = "list",
   cell_type_info = 'list',
   internal_vars = 'list',
   results = 'list',
   de_results = 'list',
   internal_vars_de = 'list'
  ),
  prototype = list(
   spatialRNA = NULL,
   originalSpatialRNA = NULL,
   reference = NULL,
   config = list(),
   cell_type_info = list(info = NULL, renorm = NULL),
   internal_vars = list(),
   results = list(),
   de_results = list(),
   internal_vars_de = list()
  )
)
