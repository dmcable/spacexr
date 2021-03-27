
#' constructor of \code{\linkS4class{Reference}} object
#' @param counts A matrix (or dgCmatrix) representing Digital Gene Expression (DGE). Rownames should be genes
#' and colnames represent barcodes/cell names.
#' @param cell_types A named (by cell barcode) factor of cell type for each cell. The 'levels' of the factor would be the possible
#' cell type identities.
#' @param nUMI Optional, a named (by cell barcode) list of total counts or UMI's appearing at each pixel. If not provided,
#' nUMI will be assumed to be the total counts appearing on each pixel.
#'
#' Counts should be untransformed count-level data
#'
#' @return Returns a \code{\linkS4class{Reference}} object containing the counts matrix, cell type labels, and UMI vector
#' from the input files
#' @export
Reference <- function(counts, cell_types, nUMI = NULL) {
  counts <- check_counts(counts, 'Reference', require_2d = T)
  if(is.null(nUMI)) {
    nUMI = colSums(counts)
    names(nUMI) <- colnames(counts)
  } else {
    check_UMI(nUMI, 'Reference', require_2d = T)
  }
  check_cell_types(cell_types)
  barcodes <- intersect(intersect(names(nUMI), names(cell_types)), colnames(counts))
  if(length(barcodes) == 0)
    stop('Reference: cell_types, counts, and nUMI do not share any barcode names. Please ensure that names(cell_types)
         matches colnames(counts) and names(nUMI)')
  if(length(barcodes) < max(length(nUMI),length(cell_types),dim(counts)[2]))
    warning('Reference: some barcodes in nUMI, cell_types, or counts were not mutually shared. Such barcodes were removed.')
  if(sum(nUMI[barcodes] != colSums(counts[,barcodes])) > 0)
    warning('Reference: nUMI does not match colSums of counts. If this is unintended, please correct this discrepancy. If this
            is intended, there is no problem.')
  missing_cell_types <- names(which(table(cell_types[barcodes]) == 0))
  if(length(missing_cell_types) > 0)
    warning(paste('Reference: missing cell types with no occurences: ',paste(missing_cell_types,collapse=', ')))
  new("Reference", cell_types = cell_types[barcodes], counts = counts[,barcodes], nUMI = nUMI[barcodes])
}

check_cell_types <- function(cell_types) {
  if(class(cell_types) != 'factor')
    stop('Reference: cell_types is not a factor. Please format cell_types as a factor.')
  if(length(cell_types) < 2)
    stop('Reference: length(cell_types) < 2. cell_types needs to be a factor with length equal to the number of cells.')
  if(length(levels(cell_types)) < 2)
    stop('Reference: length(levels(cell_types)) < 2. cell_types needs to be a factor with multiple levels for each cell type.')
  if(is.null(names(cell_types)))
    stop('Reference: names(cell_types) is null. Please enter cell barcodes as names')
}
