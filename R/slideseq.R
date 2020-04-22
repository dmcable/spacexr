
#SpatialRNA data object to contain the counts and coords information.
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

#reads counts and coords from a data director and returns a SpatialRNA object
read.SpatialRNA <- function(datadir, count_file = NULL) {
  coords <- readr::read_csv(file = paste(datadir,"BeadLocationsForR.csv",sep="/"))
  if(is.null(count_file))
    counts <- readr::read_csv(file = paste(datadir,"MappedDGEForR.csv",sep="/"))
  else
    counts <- readr::read_csv(file = paste(datadir,count_file,sep="/"))
  colnames(coords)[2] = 'x' #renaming xcoord -> x
  colnames(coords)[3] = 'y' #renaming ycoord -> y
  counts = tibble::column_to_rownames(counts, var = colnames(counts)[1])
  #rownames(counts) = counts[,1]
  coords = tibble::column_to_rownames(coords, var = "barcodes")
  #rownames(coords) <- coords$barcodes
  coords$barcodes <- NULL
  counts = counts[,2:dim(counts)[2]]
  SpatialRNA(coords, as(as(counts,"matrix"),"dgCMatrix"))
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

#restricts a puck by barcodes
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
