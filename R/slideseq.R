library(readr) #
library(tibble) #

#Slideseq data object to contain the counts and coords information.
setClass("Slideseq", 
   slots = c(
     coords = "data.frame", 
     counts = "dgCMatrix",
     n_cell_type = "integer",
     cell_type_names = "character",
     center = "numeric",
     radius = "numeric",
     nUMI = "numeric",
     cell_labels = "factor"
   ), 
   prototype = list(
     cell_type_names = NA_character_,
     n_cell_type = NA_integer_,
     coords = data.frame(NULL),
     counts = NULL,
     center = NaN,
     radius = NaN,
     nUMI = NA_integer_,
     cell_labels = factor(NULL)
   )
)

center_dist <- function(center, a, b) pointDistance(center,c(a,b),lonlat=FALSE)

#reads counts and coords from a data director and returns a slideseq object
read.slideseq <- function(datadir, count_file = NULL) {
  coords <- read_csv(file = paste(datadir,"BeadLocationsForR.csv",sep="/"))
  if(is.null(count_file))
    counts <- read_csv(file = paste(datadir,"MappedDGEForR.csv",sep="/"))
  else
    counts <- read_csv(file = paste(datadir,count_file,sep="/"))
  colnames(coords)[2] = 'x' #renaming xcoord -> x
  colnames(coords)[3] = 'y' #renaming ycoord -> y
  counts = column_to_rownames(counts, var = colnames(counts)[1])
  #rownames(counts) = counts[,1]
  coords = column_to_rownames(coords, var = "barcodes")
  #rownames(coords) <- coords$barcodes
  coords$barcodes <- NULL
  counts = counts[,2:dim(counts)[2]]
  Slideseq(coords, as(as(counts,"matrix"),"dgCMatrix"))
}

#constructor of slideseq object
Slideseq <- function(coords, counts) {
  nUMI = colSums(counts)
  center = c(mean(coords$x), mean(coords$y))
  radius = .98 * max(apply(coords[,c('x','y')], 1, function(z) center_dist(center, z['x'],z['y']))) #arbitrary shrinking factor
  new("Slideseq", coords = coords, counts = counts, center = center, radius = radius, nUMI = nUMI)
}

#Use the seurat object to create a 'fake' slideseq object
seurat.to.slideseq <- function(reference, cell_type_info) {
  cell_labels = reference@meta.data$liger_ident_coarse
  nUMI = reference@meta.data$nUMI
  counts = reference@assays$RNA@counts
  names(nUMI) = colnames(counts)
  names(cell_labels) = colnames(counts)
  cell_type_names = cell_type_info[[2]]
  n_cell_type = cell_type_info[[3]]
  new("Slideseq", counts = counts, cell_labels = cell_labels, cell_type_names = cell_type_names,
      nUMI = nUMI, n_cell_type = n_cell_type)
}

#get a uniformly spaced set of points within the range of the puck
get_uniform_points <- function(puck, delta = 50) {
  x_vals = seq(puck@center[1] - puck@radius,puck@center[1] + puck@radius,by=delta)
  y_vals = seq(puck@center[2] - puck@radius,puck@center[2] + puck@radius,by=delta)
  point_list <- expand.grid(x = x_vals, y = y_vals)
  #filter out the ones outside circle
  point_list = point_list[apply(point_list[,c('x','y')], 1, function(x) center_dist(puck@center,x[1],x[2])) < puck@radius,]
}

#get data for prediction on uniform points
get_uniform_data <- function(puck, delta = 50) {
  point_list = get_uniform_points(puck, delta)
  prediction_data = augment_prediction_data(puck,point_list)
}

#given a puck object, returns a puck with counts filtered based on UMI threshold and gene list
restrict_counts <- function(puck, gene_list, UMI_thresh = 0, UMI_max = 10000) {
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

crop_puck_circle <- function(puck, radius_mult) {
  keep_loc = (puck@coords$x - puck@center[1])^2 + (puck@coords$y - puck@center[2])^2 < (radius_mult*puck@radius)^2
  return (restrict_puck(puck, rownames(puck@coords[keep_loc,])))
}
