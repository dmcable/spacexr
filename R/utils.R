remap_celltypes <- function(cell_dict_file, cell_ident) {
  cell_type_dict <- read.csv(file=cell_dict_file, header=TRUE, sep=",")
  cell_type_dict$Name <- factor(cell_type_dict$Name)
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
  if(length(gene_list) == 0)
    stop("get_de_genes: Error: 0 common genes between SpatialRNA and Reference objects. Please check for gene list nonempty intersection.")
  gene_list = gene_list[bulk_vec[gene_list] >= MIN_OBS]
  for(cell_type in cell_type_info[[2]]) {
    other_mean = rowMeans(cell_type_info[[1]][gene_list,cell_type_info[[2]] != cell_type])
    logFC = log(cell_type_info[[1]][gene_list,cell_type] + epsilon) - log(other_mean + epsilon)
    type_gene_list = which((logFC > fc_thresh) & (cell_type_info[[1]][gene_list,cell_type] > expr_thresh)) #| puck_means[gene_list] > expr_thresh)
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

prepareBulkData <- function(cell_type_means, puck, gene_list, MIN_OBS = 10) {
  bulk_vec = rowSums(puck@counts)
  gene_list <- intersect(names(which(bulk_vec >= MIN_OBS)),gene_list)
  nUMI = sum(puck@nUMI)
  X = as.matrix(cell_type_means[gene_list,] * nUMI)
  b = bulk_vec[gene_list]
  return(list(X=X, b=b))
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
