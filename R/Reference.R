
#' constructor of \code{\linkS4class{Reference}} object
#' @param counts A matrix (or dgCmatrix) representing Digital Gene Expression (DGE). Rownames should be genes
#' and colnames represent barcodes/cell names.
#' @param cell_types A named (by cell barcode) factor of cell type for each cell. The 'levels' of the factor would be the possible
#' cell type identities.
#' @param nUMI Optional, a named (by cell barcode) list of total counts or UMI's appearing at each pixel. If not provided,
#' nUMI will be assumed to be the total counts appearing on each pixel.
#' @param min_UMI (default 100) minimum UMI count for cells to be included in the reference.
#' @param n_max_cells (default 10,000) the maximum number of cells per cell type. Will downsample if this number is exceeded.
#'
#' Counts should be untransformed count-level data
#'
#' @return Returns a \code{\linkS4class{Reference}} object containing the counts matrix, cell type labels, and UMI vector
#' from the input files
#' @export
Reference <- function(counts, cell_types, nUMI = NULL, require_int = TRUE, n_max_cells = 10000, min_UMI = 100) {
  counts <- check_counts(counts, 'Reference', require_2d = T, require_int = require_int)
  if(is.null(nUMI)) {
    nUMI = colSums(counts)
    names(nUMI) <- colnames(counts)
  } else {
    check_UMI(nUMI, 'Reference', require_2d = T, require_int = require_int, min_UMI = min_UMI)
  }
  check_cell_types(cell_types)
  barcodes <- intersect(intersect(names(nUMI), names(cell_types)), colnames(counts))
  if(length(barcodes) == 0)
    stop('Reference: cell_types, counts, and nUMI do not share any barcode names. Please ensure that names(cell_types) matches colnames(counts) and names(nUMI)')
  if(length(barcodes) < max(length(nUMI),length(cell_types),dim(counts)[2]))
    warning('Reference: some barcodes in nUMI, cell_types, or counts were not mutually shared. Such barcodes were removed.')
  barcodes <- names(which(nUMI[barcodes] >= min_UMI))
  if(length(barcodes) < 1)
    stop('Reference: no barcodes were included with nUMI at least min_UMI. Please lower the parameter min_UMI or ensure that cells have sufficient UMI counts.')
  if(sum(nUMI[barcodes] != colSums(counts[,barcodes])) > 0)
    warning('Reference: nUMI does not match colSums of counts. If this is unintended, please correct this discrepancy. If this is intended, there is no problem.')
  missing_cell_types <- names(which(table(cell_types[barcodes]) == 0))
  if(length(missing_cell_types) > 0)
    warning(paste('Reference: missing cell types with no occurences: ',paste(missing_cell_types,collapse=', ')))
  reference <- new("Reference", cell_types = cell_types[barcodes], counts = counts[,barcodes], nUMI = nUMI[barcodes])
  cur_count <- max(table(reference@cell_types))
  if(cur_count > n_max_cells) {
    warning(paste0('Reference: number of cells per cell type is ', cur_count, ', larger than maximum allowable of ', n_max_cells,
                   '. Downsampling number of cells to: ', n_max_cells))
    reference <- create_downsampled_data(reference, n_samples = n_max_cells)
  }
  return(reference)
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
  if(min(unlist(lapply(levels(cell_types),nchar))) == 0)
    stop('Reference: levels(cell_types) contains a cell type with an empty name "". Please ensure all cell type names are nonempty strings.')
  cell_type_names <- levels(cell_types)
  prohibited_character = '/'
  if(any(grepl(prohibited_character, cell_type_names)))
    stop(paste0('Reference: levels(cell_types) contains a cell type with name containing prohibited character ', prohibited_character,'. Please rename this cell type.'))
}

convert_old_reference <- function(old_reference, n_max_cells = 10000) {
  cell_types <- old_reference@meta.data$liger_ident_coarse
  nUMI <- old_reference@meta.data$nUMI
  names(cell_types) <- rownames(old_reference@meta.data);
  names(nUMI) <- rownames(old_reference@meta.data);
  Reference(old_reference@assays$RNA@counts, cell_types, nUMI, n_max_cells = n_max_cells)
}

coerce_deglam_reference <- function(old_reference) {
  return(Reference(old_reference@counts, old_reference@cell_types,
                   old_reference@nUMI, n_max_cells = max(table(old_reference@cell_types)) + 1,
                   min_UMI = 1))
}

restrict_reference <- function(reference, barcodes) {
  reference@counts <- reference@counts[,barcodes]
  reference@nUMI <- reference@nUMI[barcodes]
  reference@cell_types <- reference@cell_types[barcodes]
  return(reference)
}

restrict_reference_cell_types <- function(reference, cell_type_list) {
  new_ref <- (restrict_reference(reference, names(reference@cell_types)[reference@cell_types %in% cell_type_list]))
  new_ref@cell_types <- droplevels(new_ref@cell_types)
  return(new_ref)
}

create_downsampled_data <- function(reference, cell_types_keep = NULL, n_samples = 10000) {
  if(is.null(cell_types_keep))
    cell_types_keep = levels(reference@cell_types)
  cell_types_keep = cell_types_keep[unlist(lapply(cell_types_keep, function(x) nchar(x) > 0))]
  index_keep = c(); i = 1
  repeat{
    new_index = which(reference@cell_types == cell_types_keep[i])
    new_samples = min(n_samples, length(new_index))
    index_keep = c(index_keep, sample(new_index,new_samples,replace=FALSE))
    if((i = i + 1) > length(cell_types_keep))
      break
  }
  reference@counts = reference@counts[,index_keep]
  reference@cell_types = reference@cell_types[index_keep]
  reference@cell_types = droplevels(reference@cell_types)
  reference@nUMI = reference@nUMI[index_keep]
  return(reference)
}

save.Reference <- function(reference, save.folder) {
  dir.create(save.folder)
  write.csv(reference@cell_types, file.path(save.folder,'cell_types.csv'))
  write.csv(reference@nUMI, file.path(save.folder,'nUMI.csv'))
  write.csv(as.matrix(reference@counts), file.path(save.folder,'counts.csv'))
}

