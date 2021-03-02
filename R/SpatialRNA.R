

#' Creates a SpatialRNA object from a directory
#'
#' Given a SpatialRNA directory folder with 2 files: \code{BeadLocationsForR.csv} and \code{MappedDGEForR.csv}.
#' and returns a SpatialRNA object.
#'
#' @section Input file format (contained in datadir):
#' \enumerate{
#' \item \code{BeadLocationsForR.csv} # a CSV file (with 3 columns, with headers "barcodes", "xcoord", and "ycoord") containing the spatial locations
#' of the pixels.
#' \item \code{MappedDGEForR.csv} # a DGE (gene counts by barcodes) CSV file. Represents raw counts at each pixel.
#' }
#'
#' @param datadir (string) the directory of the SpatialRNA dataset
#' @param count_file (optional, string) the file location for the DGE
#' @return Returns a \code{\linkS4class{SpatialRNA}} object containing the coordinates and counts
#' from the input files
#' @export
read.SpatialRNA <- function(datadir, count_file = "MappedDGEForR.csv") {
  coords <- readr::read_csv(file = paste(datadir,"BeadLocationsForR.csv",sep="/"))
  counts <- readr::read_csv(file = paste(datadir,count_file,sep="/"))
  colnames(coords)[2] = 'x' #renaming xcoord -> x
  colnames(coords)[3] = 'y' #renaming ycoord -> y
  counts = tibble::column_to_rownames(counts, var = colnames(counts)[1])
  #rownames(counts) = counts[,1]
  coords = tibble::column_to_rownames(coords, var = "barcodes")
  #rownames(coords) <- coords$barcodes
  coords$barcodes <- NULL
  puck = SpatialRNA(coords, as(as(counts,"matrix"),"dgCMatrix"))
  restrict_puck(puck, colnames(puck@counts))
}


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

fake_coords <- function(counts) {
  coords <- data.frame(Matrix(0,nrow=dim(counts)[2],ncol=2))
  colnames(coords) <- c('x','y')
  rownames(coords) <- colnames(counts)
  return(coords)
}

#constructor of SpatialRNA object
SpatialRNA <- function(coords = NULL, counts, nUMI = NULL) {
  if(is.null(coords)) {
    coords <- fake_coords(counts)
  }
  if(is.null(nUMI)) {
    nUMI = colSums(counts)
  }
  names(nUMI) <- colnames(counts)
  new("SpatialRNA", coords = coords, counts = counts, nUMI = nUMI)
}

#Use the seurat object to create a 'fake' SpatialRNA object
seurat.to.SpatialRNA <- function(reference, cell_type_info) {
  cell_labels = reference@meta.data$liger_ident_coarse
  nUMI = reference@meta.data$nUMI
  counts = reference@assays$RNA@counts
  names(nUMI) = colnames(counts)
  names(cell_labels) = colnames(counts)
  cell_type_names = cell_type_info[[2]]
  n_cell_type = cell_type_info[[3]]
  coords <- fake_coords(counts)
  new("SpatialRNA", coords = coords, counts = counts, cell_labels = cell_labels,
      cell_type_names = cell_type_names, nUMI = nUMI, n_cell_type = n_cell_type)
}


#given a puck object, returns a puck with counts filtered based on UMI threshold and gene list
restrict_counts <- function(puck, gene_list, UMI_thresh = 1, UMI_max = 20000) {
  keep_loc = (puck@nUMI >= UMI_thresh) & (puck@nUMI <= UMI_max)
  puck@counts = puck@counts[gene_list,keep_loc]
  if(length(puck@cell_labels) > 0) #check cell_labels non null
    puck@cell_labels = puck@cell_labels[keep_loc]
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
  if(length(puck@cell_labels) > 0) #check cell_labels non null
    puck@cell_labels = puck@cell_labels[barcodes]
  puck@nUMI = puck@nUMI[barcodes]
  puck@coords = puck@coords[barcodes,]
  return(puck)
}

#given a dataframe with points consecutively, crops the points left (from perspective of first point) those lines
crop_puck_line <- function(puck, line_dat) {
  keep_loc = puck@coords$x == puck@coords$x
  for(i in 1:(dim(line_dat)[1]/2)) {
    y1 = line_dat$Y[2*i - 1]
    x1 = line_dat$X[2*i - 1]
    y2 = line_dat$Y[2*i]
    x2 = line_dat$X[2*i]
    keep_loc = keep_loc & ((puck@coords$y - y1)*(x2-x1) > (puck@coords$x - x1) * (y2-y1))
  }
  return (restrict_puck(puck, rownames(puck@coords[keep_loc,])))
}


split_puck <- function(puck, SpatialRNAdir, n_folds) {
  splitdir <- file.path(SpatialRNAdir, "SplitPuck")
  if(!dir.exists(splitdir))
    dir.create(splitdir)
  splitresdir <- file.path(SpatialRNAdir, "SplitPuckResults")
  if(!dir.exists(splitresdir))
    dir.create(splitresdir)
  do.call(file.remove, list(list.files(splitdir, full.names = TRUE)))
  do.call(file.remove, list(list.files(splitresdir, full.names = TRUE)))
  barcodes <- colnames(puck@counts)
  N <- length(barcodes)
  for (i in 1:n_folds) {
    puck_fold <- restrict_puck(puck, barcodes[(round((i-1)*N/n_folds) + 1):round(i*N/n_folds)])
    saveRDS(puck_fold, file.path(splitdir, paste0("puck",i,".RDS")))
  }
}

#coerces an old SpatialRNA object
coerce_old <- function(puck) {
  new("SpatialRNA", coords = puck@coords, counts = puck@counts, cell_labels = puck@cell_labels,
      cell_type_names = puck@cell_type_names, nUMI = puck@nUMI, n_cell_type = puck@n_cell_type)
}


