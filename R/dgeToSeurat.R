#' Creates a scRNA-seq Seurat object
#'
#' Given a reference directory folder with 3 files: \code{dge.csv}, \code{meta_data.csv},
#' and \code{cell_type_dict.csv}. Saves the reference as refdir/'SCRef.RDS' and also
#' saves a reference downsampled to 1000 cells (per cell type) as refdir/'scRefSubsampled.RDS'
#'
#' @section Input file format (contained in refdir):
#' \enumerate{
#' \item \code{meta_data.csv} # a CSV file (with 3 columns, with headers "barcode", "cluster", and "nUMI") containing the numeric cluster assignment for each cell.
#' \item \code{cell_type_dict.csv} # a CSV file (with 2 columns, with headers "Cluster" and "Name") containing the mapping between numeric cluster ID and cluster name. If you want a cluster to be filtered out of the single cell reference, you can leave the cluster name blank. The cell types must not contain the character '/' or '-'.
#' \item \code{dge.csv} # a DGE (barcodes by gene counts) CSV file in the standard 10x format.
#' }
#'
#' @param refdir (string) the directory of the reference
#' @return Returns a \code{\link{Seurat}} object containing the scRNA-seq reference, downsampled to 1000 cells
#' (per cell type).
#' @export
dgeToSeurat <- function(refdir) {
  #check for , at end of header of DGE
  dge_file = file.path(refdir,"dge.csv")
  conn <- file(dge_file,open="r")
  lines <-readLines(conn)
  text_len = nchar(lines[1])
  if(substr(lines[1],text_len,text_len) == ",")
   lines[1] <- substr(lines[1],1,text_len-1)
  writeLines(lines,dge_file)
  close(conn)
  dge <- readr::read_csv(file.path(refdir,"dge.csv"))
  gene_list <- dge$X1
  raw.data <- as(as(dge[,-1],"matrix"),"dgCMatrix")
  rownames(raw.data) <- gene_list
  meta_data = read.csv(file.path(refdir,"meta_data.csv"))
  rownames(meta_data) = meta_data$barcode
  meta_data$barcode = NULL
  common_barcodes = intersect(colnames(raw.data), rownames(meta_data))
  raw.data = raw.data[,common_barcodes]
  meta_data = meta_data[common_barcodes, ]
  cell_dict_file <- paste(refdir,"cell_type_dict.csv",sep="/")
  true_type_names <- remap_celltypes(cell_dict_file, meta_data$cluster)
  meta_data$liger_ident_coarse = true_type_names
  reference = Seurat::CreateSeuratObject(raw.data, meta.data = meta_data)
  saveRDS(reference, paste(refdir,"SCRef.RDS",sep="/"))
  dref <- create_downsampled_data(reference, refdir)
  return(dref)
}

#create down sampled data
#if each cell_type, then it takes n_samples from each cell type. Otherwise, it takes randomly from whole dataset
create_downsampled_data <- function(reference, refdir, cell_type_import = NULL,n_samples = 1000, each_cell_type = T,save.file = T) {
  if(!each_cell_type)
    index_keep = sample(which(reference@meta.data$liger_ident_coarse != ""),n_samples,replace=FALSE)
  else {
    if(is.null(cell_type_import))
      cell_type_import = levels(reference@meta.data$liger_ident_coarse)
    cell_type_import = cell_type_import[unlist(lapply(cell_type_import, function(x) nchar(x) > 0))]
    index_keep = c(); i = 1
    repeat{
      new_index = which(reference@meta.data$liger_ident_coarse == cell_type_import[i])
      new_samples = min(n_samples, length(new_index))
      index_keep = c(index_keep, sample(new_index,new_samples,replace=FALSE))
      if((i = i + 1) > length(cell_type_import))
        break
    }
  }
  reference@assays$RNA@counts = reference@assays$RNA@counts[,index_keep]
  reference@meta.data = reference@meta.data[index_keep,]
  reference@meta.data$liger_ident_coarse = droplevels(reference@meta.data$liger_ident_coarse)
  reference@assays$RNA@data <-matrix(0,nrow=2,ncol=2)
  reference@assays$RNA@scale.data <- matrix(0,nrow=2,ncol=2)
  if(save.file)
    saveRDS(reference, paste(refdir,"/scRefSubsampled", n_samples, ".RDS",sep=""))
  return(reference)
}
