get_cell_type_info <- function(raw.data, cell_types, nUMI, cell_type_names = NULL) {
  if(is.null(cell_type_names))
    cell_type_names = levels(cell_types)

  n_cell_types = length(cell_type_names)

  get_cell_mean <- function(cell_type) {
    cell_type_data = raw.data[,cell_types == cell_type]
    cell_type_umi = nUMI[cell_types == cell_type]
    normData = sweep(cell_type_data,2,cell_type_umi,`/`)
    return(rowSums(normData) / dim(normData)[2])
  }

  cell_type = cell_type_names[1]
  cell_type_means <- data.frame(get_cell_mean(cell_type))
  colnames(cell_type_means)[1] = cell_type
  for (cell_type in cell_type_names[2:length(cell_type_names)]) {
    cell_type_means[cell_type] = get_cell_mean(cell_type)
  }
  return(list(cell_type_means, cell_type_names, n_cell_types))
}

#renormalizes cell_type_means to have average the same as the puck
#proportions is the estimated proportion of each cell type
get_norm_ref <- function(puck, cell_type_means, gene_list, proportions) {
  bulk_vec = rowSums(puck@counts)
  weight_avg = rowSums(sweep(cell_type_means[gene_list,],2,proportions / sum(proportions),'*'))
  target_means = bulk_vec[gene_list]/sum(puck@nUMI)
  cell_type_means_renorm = sweep(cell_type_means[gene_list,],1,weight_avg / target_means,'/')
}

#restrict a reference to certain cell types
restrict_ref <- function(reference, cell_type_names) {
  keep = reference@meta.data$liger_ident_coarse %in% cell_type_names
  reference@meta.data = reference@meta.data[keep,]
  reference@assays$RNA@counts = reference@assays$RNA@counts[,keep]
  return(reference)
}
