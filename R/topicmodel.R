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

test_doublet_beads <- function(puck, gene_list, cell_type_info, SINGLET_THRESH = 0.8) {
  weights <- process_beads_batch(cell_type_info, gene_list, puck)
  first_label = names(unlist(lapply(weights,function(x) which.max(x))))
  sec_label = names(unlist(lapply(weights,function(x) which.min(x))))
  singlet = unlist(lapply(weights,function(x) as.numeric(max(x) > SINGLET_THRESH)))
  max_weight = unlist(lapply(weights,function(x) max(x)))
  results_df <- data.frame(first_label, sec_label, singlet, max_weight)
}

#main function for assigning cell type labels and decompositions to a dataset.
#if proportions is null, does not renormalize cell type means
process_data <- function(puck, gene_list, cell_type_info, proportions = NULL, trust_model = FALSE, constrain = T) {
  cell_type_info_renorm = cell_type_info
  if(!is.null(proportions)) {
    proportions = smooth_proportions(proportions, smoothing_par = 0, constrain = F)
    cell_type_info_renorm[[1]] = get_norm_ref(puck, cell_type_info[[1]], gene_list, proportions)
  }
  test_results <- test_single_beads(puck, gene_list, cell_type_info_renorm, trust_model = trust_model, constrain = constrain)
  return(test_results)
}
