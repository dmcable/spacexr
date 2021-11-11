read.Slideseq.mtx <- function(datadir, puck_no) {
  mtx_file <- file.path(datadir, paste0(puck_no, '.matched.digital_expression_matrix.mtx'))
  barcodes_file <- file.path(datadir,  paste0(puck_no, '.matched.digital_expression_barcodes.tsv'))
  coords_file <- file.path(datadir,  paste0(puck_no, '_barcode_matching.txt' ))
  features_file <- file.path(datadir,  paste0(puck_no, '.matched.digital_expression_features.tsv'))
  counts <- readMM(mtx_file)
  barcodes <- read.csv(barcodes_file, header = F)
  pre_coords <- read.csv(coords_file, header = F, sep = '\t')
  pre_coords <- pre_coords[!duplicated(pre_coords[,2]), ]
  coords <- data.frame(pre_coords[,3:4])
  rownames(coords) <- pre_coords[,2]
  features <- read.csv(features_file, header = F, sep = '\t')
  rownames(counts) <- features[,2]
  colnames(counts) <- barcodes[,1]
  puck <- SpatialRNA(coords,counts)
  return(puck)
}
library(fields)
crop.puck.convex <- function(puck, puck_old, puckCropped_old, DIST.THRESH = 50, DIST.THRESH.BAD = 50) {
  chull <- convhulln(puckCropped_old@coords)
  in_region <- inhulln(chull, as.matrix(puck@coords))
  names(in_region) <- rownames(puck@coords)
  d_mat <- rdist(puck@coords, puckCropped_old@coords)
  d_mat_bad <- rdist(puck@coords, puck_old@coords[! (rownames(puck_old@coords) %in% rownames(puckCropped_old@coords)),])
  min_dist <- apply(d_mat,1,min)
  max_dist <- apply(d_mat_bad, 1, min)
  Y <- as.integer(in_region) #Oligo
  names(Y) <- names(in_region)
  Y[names(Y) %in% rownames(puckCropped_old@coords)] <- Y[names(Y) %in% rownames(puckCropped_old@coords)] + 1
  Y[(min_dist > DIST.THRESH) & (max_dist < DIST.THRESH.BAD)] <- 0
  p <- plot_puck_continuous(puck,rownames(puck@coords), Y, ylimit = c(0,2)) #Oligo
  print(p)
  puckCropped = restrict_puck(puck,names(which(Y > 0)))
  return(puckCropped)
}
