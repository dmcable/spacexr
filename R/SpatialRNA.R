
fake_coords <- function(counts) {
  coords <- data.frame(Matrix(0,nrow=dim(counts)[2],ncol=2))
  colnames(coords) <- c('x','y')
  rownames(coords) <- colnames(counts)
  return(coords)
}

#' constructor of SpatialRNA object
#' @param counts A matrix (or dgCmatrix) representing Digital Gene Expression (DGE). Rownames should be genes
#' and colnames represent barcodes/pixel names.
#' @param coords A data.frame (or matrix) representing the spatial pixel locations. rownames are barcodes/pixel names,
#' and there should be two columns for 'x' and for 'y'.
#' @param nUMI Optional, a named (by pixel barcode) list of total counts or UMI's appearing at each pixel. If not provided,
#' nUMI will be assumed to be the total counts appearing on each pixel.
#' @param use_fake_coords logical, FALSE by default. If true, the 'coords' parameter will be ignored, and replaced with a placeholder
#' coords matrix.
#'
#' Counts should be untransformed count-level data
#'
#' @return Returns a \code{\linkS4class{SpatialRNA}} object containing the coordinates and counts
#' from the input files
#' @export


#' Creates a SpatialRNA object from a 10x Genomics Visium `outs` directory
#'
#' Given a SpatialRNA directory 10x Genomics Visium `outs` directory and returns a SpatialRNA object.
#'
#' @param datadir (string) full path to the 10x Genomics Visium `outs` directory
#' @return Returns a \code{\linkS4class{SpatialRNA}} object containing the coordinates and counts
#' from the input files
#' @export
read.VisiumSpatialRNA <- function (datadir) 
{
  coords <- readr::read_csv(file = paste(datadir, "spatial/tissue_positions_list.csv", 
                                         sep = "/"), 
                            col_names = c("barcodes", "in_tissue", "x", "y", "pxl_col_in_fullres", "pxl_row_in_fullres"))
  coords = tibble::column_to_rownames(coords, var = "barcodes")
  counts <- Seurat::Read10X_h5(paste0(datadir, "/filtered_feature_bc_matrix.h5"))
  puck = SpatialRNA(coords, counts)
  restrict_puck(puck, colnames(puck@counts))
}


save.SpatialRNA <- function(puck, save.folder) {
  dir.create(save.folder)
  write.csv(puck@coords, file.path(save.folder,'coords.csv'))
  white.csv(puck@nUMI, file.path(save.folder,'nUMI.csv'))
}

SpatialRNA <- function(coords, counts, nUMI = NULL, use_fake_coords = FALSE) {
  counts <- check_counts(counts, 'SpatialRNA')
  if(use_fake_coords)
    coords <- use_fake_coords(counts)
  else
    coords <- check_coords(coords)
  if(is.null(nUMI)) {
    nUMI = colSums(counts)
    names(nUMI) <- colnames(counts)
  } else {
    check_UMI(nUMI, 'SpatialRNA')
  }
  barcodes <- intersect(intersect(names(nUMI), rownames(coords)), colnames(counts))
  if(length(barcodes) == 0)
    stop('SpatialRNA: coords, counts, and nUMI do not share any barcode names. Please ensure that rownames(coords)
         matches colnames(counts) and names(nUMI)')
  if(length(barcodes) < max(length(nUMI),dim(coords)[1],dim(counts)[2]))
    warning('SpatialRNA: some barcodes in nUMI, coords, or counts were not mutually shared. Such barcodes were removed.')
  if(sum(nUMI[barcodes] != colSums(counts[,barcodes])) > 0)
    warning('SpatialRNA: nUMI does not match colSums of counts. If this is unintended, please correct this discrepancy. If this
            is intended, there is no problem.')
  new("SpatialRNA", coords = coords[barcodes,], counts = counts[,barcodes], nUMI = nUMI[barcodes])
}

check_UMI <- function(nUMI, f_name, require_2d = F) {
  if(!is.atomic(nUMI))
    stop(paste0(f_name,': nUMI is not an atomic vector. Please format nUMI as an atomic vector.'))
  if(!is.numeric(nUMI))
    stop(paste0(f_name,': nUMI is not numeric'))
  if(max(abs(nUMI %% 1)) > 1e-6)
    stop(paste0(f_name,': nUMI does not contain integers'))
  if(is.null(names(nUMI)))
    stop(paste0(f_name,': names(nUMI) is null. Please enter barcodes as names'))
  if(length(nUMI) == 1)
    if(require_2d)
      stop(paste0(f_name,': the length of nUMI is 1, indicating only one cell present. Please format nUMI so that
           the length is greater than 1.'))
    else
      warning(paste0(f_name,': the length of nUMI is 1, indicating only one cell present. If this is unintended,
        please format nUMI so that the length is greater than 1.'))
}

check_counts <- function(counts, f_name, require_2d = F) {
  if(class(counts) != 'dgCMatrix') {
    if(class(counts) != 'matrix')
      tryCatch({
        counts <- as(counts,'matrix')
      }, error = function(e) {
        stop(paste0(f_name,': could not convert counts to matrix using as(counts,\'matrix\'). Please check that
             counts is coercible to matrix, such as a matrix, dgCmatrix, or data.frame.'))
      })
    counts <- as(counts,"dgCMatrix")
  }
  if(dim(counts)[1] == 1) #check more than one gene
    stop(paste0(f_name,': the first dimension of counts is 1, indicating only one gene present. Please format counts so that
           the first dimension is greater than 1.'))
  if(dim(counts)[2] == 1)
    if(require_2d)
      stop(paste0(f_name,': the second dimension of counts is 1, indicating only one cell present. Please format counts so that
           the second dimension is greater than 1.'))
    else
      warning(paste0(f_name,': the second dimension of counts is 1, indicating only one cell/pixel present. If this is unintended,
        please format counts so that the second dimension is greater than 1.'))
  if(!is.numeric(counts[1,1]))
    stop(paste0(f_name,': elements of counts are not numeric'))
  if(max(abs(counts %% 1)) > 1e-6)
    stop(paste0(f_name,': counts does not contain integers'))
  if(is.null(rownames(counts)))
    stop(paste0(f_name,': rownames(counts) is null. Please enter gene names as rownames'))
  if(is.null(colnames(counts)))
    stop(paste0(f_name,': colnames(counts) is null. Please enter barcodes as colnames'))
  return(counts)
}

check_coords <- function(coords) {
  if(class(coords) != 'data.frame') {
    tryCatch({
      coords <- as(coords,'data.frame')
    }, error = function(e) {
      stop('SpatialRNA: could not convert coords to data.frame using as(coords,\'data.frame\'). Please check that
           coords is coercible to data.frame, such as a matrix or data.frame .')
    })
  }
  if(dim(coords)[2] != 2) #check more than one gene
    stop('SpatialRNA: the second dimension of coords is not 2. Please enforce that dim(coords)[2] == 2 (x and y coordinates).')
  colnames(coords) <- c('x','y')
  if(!(is.numeric(coords$x) & is.numeric(coords$y)))
    stop('SpatialRNA: coords is not numeric')
  if(is.null(rownames(coords)))
    stop('SpatialRNA: rownames(coords) is null. Please enter barcodes as rownames')
  return(coords)
}



#' Restricts a SpatialRNA object to a subset of genes (and applies a UMI threshold)
#'
#' @param puck a \code{\linkS4class{SpatialRNA}} object
#' @param gene_list a list of gene names
#' @param UMI_thresh minimum UMI per pixel
#' @param UMI_max maximum UMI per pixel
#' @return Returns a \code{\linkS4class{SpatialRNA}} with counts filtered based on UMI threshold and gene list
#' @export
restrict_counts <- function(puck, gene_list, UMI_thresh = 1, UMI_max = 20000) {
  keep_loc = (puck@nUMI >= UMI_thresh) & (puck@nUMI <= UMI_max)
  puck@counts = puck@counts[gene_list,keep_loc]
  puck@nUMI = puck@nUMI[keep_loc]
  return(puck)
}

#' Restricts a SpatialRNA object to a subset of pixels
#'
#' Given a \code{\linkS4class{SpatialRNA}} object and a list of barcodes (pixels), will return a
#' \code{\linkS4class{SpatialRNA}} object restricted to the barcodes.
#'
#' @param puck a \code{\linkS4class{SpatialRNA}} object
#' @param barcodes a list of barcode names, a subset of \code{rownames(puck@coords)}
#' @return Returns a \code{\linkS4class{SpatialRNA}} object subsampled to the barcodes
#' @export
restrict_puck <- function(puck, barcodes) {
  barcodes = intersect(colnames(puck@counts), barcodes)
  puck@counts = puck@counts[,barcodes]
  puck@nUMI = puck@nUMI[barcodes]
  puck@coords = puck@coords[barcodes,]
  return(puck)
}


#coerces an old SpatialRNA object
coerce_old <- function(puck) {
  new("SpatialRNA", coords = puck@coords, counts = puck@counts, nUMI = puck@nUMI)
}


