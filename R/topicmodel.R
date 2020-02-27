test_single_beads <- function(puck, gene_list, cell_type_info, trust_model = FALSE, constrain = T, OLS = F) {
  cell_type_names = cell_type_info[[2]]; n_cell_types = cell_type_info[[3]]
  beads = t(as.matrix(puck@counts[gene_list,]))
  weights = decompose_batch(puck@nUMI, cell_type_info[[1]], beads, gene_list, constrain = constrain, OLS = OLS)
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
process_data <- function(puck, gene_list, cell_type_info, proportions = NULL, trust_model = FALSE, constrain = T, OLS = F) {
  cell_type_info_renorm = cell_type_info
  if(!is.null(proportions)) {
    proportions = smooth_proportions(proportions, smoothing_par = 0, constrain = F)
    cell_type_info_renorm[[1]] = get_norm_ref(puck, cell_type_info[[1]], gene_list, proportions)
  }
  test_results <- test_single_beads(puck, gene_list, cell_type_info_renorm, trust_model = trust_model, constrain = constrain, OLS = OLS)
  return(test_results)
}


#calculates the quasilikelihood function as a matrix
get_Q_mat <- function(delta = 1e-5, L = 100000, mult = 100, Y_max = 100) {
  Q_mat = matrix(0, nrow = Y_max, ncol = 2*L)
  for (Y in 0:(Y_max-1)) {
    x = (1:(L*mult)) * delta
    V <- get_V(x)
    sigma <- sqrt(V)
    sigma_bar <- pmax(sigma, 1)
    scaled_residual <- (Y - x)/sigma_bar
    results <- sigma_bar*phi(scaled_residual)/(sigma^2)
    precise_J <- cumsum(results)*delta
    precise_J <- precise_J - precise_J[round(Y/delta+1)]
    Q_mat[Y+1,1:L] = precise_J[1:L]
    Q_mat[Y+1,(L+1):(2*L)] = precise_J[(1:L)*mult]
  }
  return(Q_mat)
}

#assumes Q_mat was initialized as a global variable using get_Q_mat with same params for delta, L, mult, and Y_max
get_QL <- function(Y, x, delta = 1e-5, L = 100000, mult = 100, Y_max = 100) {
  if(Y > Y_max - 1)
    Y = Y_max - 1
  my_index = x/delta
  if(my_index < 1)
    my_index = 1
  if(my_index <= L) {
    first_ind = floor(my_index)
    other_p = my_index - first_ind
    return ((1-other_p)*Q_mat[Y+1,first_ind] + other_p*Q_mat[Y+1,first_ind+1])
  } else {
    my_index = my_index / mult
    if(my_index >= L)
      my_index = L - 1e-9
    first_ind = floor(my_index)
    other_p = my_index - first_ind
    return ((1-other_p)*Q_mat[Y+1,L + first_ind] + other_p*Q_mat[Y+1,L + first_ind+1])
  }
}

