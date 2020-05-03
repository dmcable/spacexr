remap_celltypes <- function(cell_dict_file, cell_ident) {
  cell_type_dict <- read.csv(file=cell_dict_file, header=TRUE, sep=",")
  rownames(cell_type_dict) = cell_type_dict[,'Cluster']
  true_type_names = lapply(cell_ident, function(x) cell_type_dict[as.character(x),"Name"])
  true_type_names = unlist(true_type_names)
}

#finds DE genes
#Genes must be observed a minimum of MIN_OBS times to mitigate sampling noise in the
#Platform effect estimation

#' Returns a list of differentially expressed genes
#'
#' For each cell type, chooses genes that have a minimum average normalized expression in that cell
#' type, and whose expression is larger in that cell type than the average of all cell types.
#' Filters out mitochondrial genes.
#'
#' @param puck an object of type \linkS4class{SpatialRNA}
#' @param cell_type_info cell type information and profiles of each cell, calculated from the scRNA-seq
#' reference (see \code{\link{get_cell_type_info}})
#' @param MIN_OBS the minimum number of occurances of each gene in the SpatialRNA object.
#' @param fc_thresh minimum \code{log_e} fold change required for a gene.
#' @param expr_thresh minimum expression threshold, as normalized expression (proportion out of 1, or counts per 1).
#' @return a list of differntially expressed gene names
#' @export
get_de_genes <- function(cell_type_info, puck, fc_thresh = 1.25, expr_thresh = .00015, MIN_OBS = 3) {
  total_gene_list = c()
  epsilon = 1e-9
  bulk_vec = rowSums(puck@counts)
  gene_list = rownames(cell_type_info[[1]])
  if(length(grep("mt-",gene_list)) > 0)
    gene_list = gene_list[-grep("mt-",gene_list)]
  gene_list = intersect(gene_list,names(bulk_vec))
  gene_list = gene_list[bulk_vec[gene_list] >= MIN_OBS]
  for(cell_type in cell_type_info[[2]]) {
    other_mean = rowMeans(cell_type_info[[1]][gene_list,cell_type_info[[2]] != cell_type])
    logFC = log(cell_type_info[[1]][gene_list,cell_type] + epsilon) - log(other_mean + epsilon)
    type_gene_list = which((logFC > fc_thresh) & (cell_type_info[[1]][gene_list,cell_type] > expr_thresh))
    print(paste0("get_de_genes: ", cell_type, " found DE genes: ",length(type_gene_list)))
    total_gene_list = union(total_gene_list, type_gene_list)
  }
  total_gene_list = gene_list[total_gene_list]
  print(paste0("get_de_genes: total DE genes: ",length(total_gene_list)))
  return(total_gene_list)
}

#min weight to be considered a singlet as a function of nUMI
UMI_cutoff <- function(nUMI) {
  return (pmax(0.25, 2 - log(nUMI,2) / 5))
}

prepareBulkData <- function(bulkdir, cell_type_means, puck, gene_list) {
  bulk_vec = rowSums(puck@counts)
  nUMI = sum(puck@nUMI)
  X = cell_type_means[gene_list,] * nUMI
  b = bulk_vec[gene_list]
  write.csv(as.matrix(X),file.path(bulkdir,"X_bulk.csv"))
  write.csv(as.matrix(b),file.path(bulkdir,"b_bulk.csv"))
  return(list(X=X, b=b))
}

#if test_reference is not null, it will convert this to a SpatialRNA object rather than read in one
#if puck_file is not null, then reads in puck from this file
#if load_info_renorm, loads cell type info from MetaData/cell_type_info_renorm.RDS. Takes gene_list to be rownames of cell_type_info (renorm)
#get_proportions -> calculates cell type info renorm

#' Initializes RCTD, returns meta-data
#'
#' Returns a named list of data used for running RCTD. Reads in config file
#' \code{conf/dataset.yml}, and then reads in the conf/config_mode.yml file
#' using the config_mode from the first config file. These config files must be configured
#' to point RCTD to the correct datasets, and to provide parameters for RCTD.
#' Will search for scRNA-seq reference in 'reffolder/reffile'.
#' Will search for SpatialRNA data in 'SpatialRNAfolder/puckrds'
#'
#' @param gene_list_reg logical, if T uses the 'reg' cutoffs for gene list, rather
#' than the bulk cutoffs.
#' @param get_proportions logical, whether to load proportions from 'results/Bulk/weights.RDS'
#' @param load_info_renorm logical, if true, reads \code{cell_type_info} from \code{MetaData/cell_type_info_renorm.RDS}
#' @return A named list of meta data, including
#' \itemize{
#'   \item refdir, the directory of the scRNA-seq reference
#'   \item SpatialRNAdir, the directory of the SpatialRNA data
#'   \item bulkdir, the directory of results for the bulk analysis
#'   \item reference, a \code{\link{Seurat}} object of the scRNA-seq reference
#'   \item proportions, a a named list (for each cell type) of proportion of the cell type on the bulk dataset
#'   \item gene_list, a list of genes to be used for the normalization
#'   \item puck, an object of type \linkS4class{SpatialRNA}
#'   \item cell_type_info, cell type information and profiles of each cell, calculated from the scRNA-seq
#' reference (see \code{\link{get_cell_type_info}})
#'   \item config, additional configuration parameters
#'   \item n_puck_folds, the number of folds to split the puck for RCTD
#'   \item puckrds, the filename of the SpatialRNA object
#' }
#' @export
init_RCTD <- function(gene_list_reg = T, get_proportions = F, test_reference = NULL, puck_file = NULL, MIN_OBS = 1, load_info_renorm = F, load_info = F) {
  print("init_RCTD: begin")
  config_data <- config::get(file = "conf/dataset.yml", use_parent = FALSE)
  print(paste("init_RCTD: using config mode:",config_data$config_mode))
  config <- config::get(file = paste0("conf/",config_data$config_mode,".yml"), use_parent = FALSE)
  SpatialRNAdir <- config_data$SpatialRNAfolder
  resultsdir = file.path(SpatialRNAdir,"results")
  if(!dir.exists(resultsdir))
    dir.create(resultsdir)
  bulkdir <- paste(SpatialRNAdir,"results/Bulk",sep="/")
  if(!dir.exists(bulkdir))
    dir.create(bulkdir)
  if(!dir.exists('logs'))
    dir.create('logs')
  refdir <- config_data$reffolder
  if(load_info_renorm) {
    cell_type_info <- readRDS(file.path(SpatialRNAdir, "MetaData/cell_type_info_renorm.RDS"))
    reference <- NULL;
  } else if(load_info) {
    cell_type_info <- readRDS(file.path(refdir, "MetaData/cell_type_info.RDS"))
    reference <- NULL;
  } else {
    reference <- readRDS(paste(refdir,config_data$reffile,sep="/"))
    print(paste("init_RCTD: number of cells in reference:", dim(reference@assays$RNA@counts)[2]))
    print(paste("init_RCTD: number of genes in reference:", dim(reference@assays$RNA@counts)[1]))
    cell_counts = table(reference@meta.data$liger_ident_coarse)
    print(cell_counts)
    CELL_MIN = 25 # need at least this for each cell type
    if(min(cell_counts) < CELL_MIN)
      stop(paste0("init_RCTD error: need a minimum of ",CELL_MIN, " cells for each cell type in the reference"))
    cell_type_info <- get_cell_type_info(reference@assays$RNA@counts, reference@meta.data$liger_ident_coarse, reference@meta.data$nUMI)
  }
  print(paste("init_RCTD: number of cell types used:", cell_type_info[[3]]))
  proportions <- NULL
  if(get_proportions) {
    proportions <- readRDS(file.path(bulkdir,"weights.RDS"))
    names(proportions) = cell_type_info[[2]]
    proportions <- proportions # / sum(proportions)
    print("init_RCTD: estimated bulk composition: ")
    print(proportions)
  }
  if(is.null(test_reference)) {
    if(is.null(puck_file))
      puck = readRDS(file.path(SpatialRNAdir, config_data$puckrds))
    else
      puck = readRDS(file.path(SpatialRNAdir, puck_file))
  }
  else
    puck <- seurat.to.SpatialRNA(test_reference, cell_type_info)
  puck = restrict_counts(puck, rownames(puck@counts), UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
  if(load_info_renorm)
    gene_list = rownames(cell_type_info[[1]])
  else {
    if(gene_list_reg)
      gene_list = get_de_genes(cell_type_info, puck, fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg, MIN_OBS = MIN_OBS)
    else
      gene_list = get_de_genes(cell_type_info, puck, fc_thresh = config$fc_cutoff, expr_thresh = config$gene_cutoff, MIN_OBS = MIN_OBS)
  }
  print(paste("init_RCTD: number of genes used:", length(gene_list)))
  puck = restrict_counts(puck, gene_list, UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
  puck = restrict_puck(puck, colnames(puck@counts))
  print(paste("init_RCTD: number of pixels used in test data passing UMI threshold:", dim(puck@counts)[2]))
  print("init_RCTD: end")
  return(list(refdir = refdir, SpatialRNAdir = SpatialRNAdir, bulkdir = bulkdir, reference = reference,
              proportions = proportions, gene_list = gene_list, puck = puck, cell_type_info = cell_type_info,
              config = config, n_puck_folds = config_data$n_puck_folds, puckrds = config_data$puckrds))
}


get_class_df <- function(cell_type_names, use_classes = F) {
  class_df = data.frame(cell_type_names, row.names = cell_type_names)
  colnames(class_df)[1] = "class"
  if(use_classes) {
    class_df["Bergmann","class"] = "Astrocytes"
    class_df["Fibroblast","class"] = "Endothelial"
    class_df["MLI2","class"] = "MLI1"
    class_df["Macrophages","class"] = "Microglia"
    class_df["Polydendrocytes","class"] = "Oligodendrocytes"
  }
  return(class_df)
}
