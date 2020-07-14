#' An S4 class to represent Spatial Transcriptomic data
#'
#' @slot coords a dataframe with x and y coordinates of each pixel
#' @slot counts a dataframe of raw counts for each gene (rowname)
#' and each pixel (colnames or barcodes)
#' @slot n_cell_type the number of cell types
#' @slot cell_type_names a list of cell type names
#' @slot nUMI a named list (by barcode) of total UMIs per pixel
#' @slot cell_labels a factor of cell type labels for each pixel
#' @export
#' @import Matrix
#' @import doParallel
#' @import foreach
#' @importClassesFrom Seurat Seurat
#' @importClassesFrom Matrix Matrix dgCMatrix
setClass("SpatialRNA",
         slots = c(
           coords = "data.frame",
           counts = "dgCMatrix",
           n_cell_type = "integer",
           cell_type_names = "character",
           nUMI = "numeric",
           cell_labels = "factor"
         ),
         prototype = list(
           cell_type_names = NA_character_,
           n_cell_type = NA_integer_,
           coords = data.frame(NULL),
           counts = NULL,
           nUMI = NA_integer_,
           cell_labels = factor(NULL)
         )
)

#' An S4 class used to run the RCTD algorithm
#'
#' Created using the \code{\link{create.RCTD}} function, a user can run RCTD using the \code{\link{run.RCTD}} function.
#'
#' @slot spatialRNA a \code{\linkS4class{SpatialRNA}} object containing the Spatial RNA dataset to be processed
#' @slot reference a Seurat object containing the annotated single cell reference
#' @slot config a list of configuration options, set using the \code{\link{create.RCTD}} function
#' @slot cell_type_info a named list of cell type profiles (means), containing two elements: \code{info}, directly calculated from the scRNA-seq reference, and
#' \code{renorm}, renormalized the match the SpatialRNA dataset.
#' @slot internal_vars a list of internal variables used by RCTD's computation
#' @slot results (created after running RCTD) a list of results_df (a dataframe of RCTD results in doublet mode),
#' weights (a dataframe of RCTD predicted weights in full mode), and weights_doublet (a
#' dataframe of predicted weights in doublet mode, with cell type information in results_df)
#' @export
#' @import Matrix
#' @import doParallel
#' @import foreach
#' @importClassesFrom Seurat Seurat
#' @importClassesFrom Matrix Matrix dgCMatrix
setClass("RCTD",
         slots = c(
           spatialRNA = "SpatialRNA",
           reference = "Seurat",
           config = "list",
           cell_type_info = 'list',
           internal_vars = 'list',
           results = 'list'
         ),
         prototype = list(
           spatialRNA = NULL,
           reference = NULL,
           config = list(),
           cell_type_info = list(info = NULL, renorm = NULL),
           internal_vars = list(),
           results = list()
         )
)
