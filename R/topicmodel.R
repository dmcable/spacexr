test_single_beads <- function(puck, gene_list, cell_type_info, trust_model = FALSE, constrain = T) {
  cell_type_names = cell_type_info[[2]]; n_cell_types = cell_type_info[[3]]
  beads = t(as.matrix(puck@counts[gene_list,]))
  weights = decompose_batch(puck@nUMI, cell_type_info[[1]], beads, gene_list, constrain = constrain)
  pred_labels = unlist(lapply(weights,function(x) which.max(x)))
  cell_type_lev = factor(1:n_cell_types)
  cell_type_map = data.frame(cindex = 1:n_cell_types, row.names = cell_type_names)
  if(trust_model)
    true_labels = pred_labels
  else
    true_labels = lapply(puck@cell_labels, function(x) cell_type_map[as.character(x),"cindex"])
  conf_mat = caret::confusionMatrix(factor(pred_labels,cell_type_lev),factor(true_labels,cell_type_lev))
  rownames(conf_mat$table) = cell_type_names; colnames(conf_mat$table) = cell_type_names
  return(list(conf_mat,weights,pred_labels))
}


#main function for assigning cell type labels and decompositions to a dataset.
process_data <- function(puck, gene_list, cell_type_info, proportions, trust_model = FALSE, constrain = T) {
  proportions = smooth_proportions(proportions, constrain = F)
  cell_type_info_renorm = cell_type_info
  cell_type_info_renorm[[1]] = get_norm_ref(puck, cell_type_info[[1]], gene_list, proportions)
  test_results <- test_single_beads(puck, gene_list, cell_type_info_renorm, trust_model = trust_model, constrain = constrain)
  return(test_results)
}
