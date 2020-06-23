#decompose with just two cell types
#if score_mode, then returns the objective function score
#if denoise, then it fits a "noise" dimension as the mean of all the data
decompose_sparse <- function(cell_type_profiles, nUMI, bead, type1=NULL, type2=NULL, score_mode = FALSE, plot = F, custom_list = NULL, verbose=F, constrain = T) {
  if(is.null(custom_list))
    cell_types = c(type1,type2)
  else
    cell_types = custom_list
  reg_data = data.matrix(cell_type_profiles[,cell_types])
  if(score_mode)
    n.iter = 25
  else
    n.iter = 50
  results = solveIRWLS.weights(reg_data,bead,nUMI,OLS = FALSE, constrain = constrain, verbose = verbose, n.iter = n.iter)
  if(! score_mode) {
    results$weights = results$weights / sum(results$weights)
    return(results)
  } else {
    prediction = reg_data %*% results$weights
    total_score = calc_log_l_par(gene_list, prediction, bead)
    return (total_score)
  }
}

#decompose with all cell types
decompose_full <- function(cell_type_profiles, nUMI, bead, constrain = TRUE, OLS = FALSE, verbose = F, n.iter = 50, MIN_CHANGE = 0.001) {
  results = solveIRWLS.weights(cell_type_profiles,bead,nUMI,OLS = OLS, constrain = constrain, verbose = verbose, n.iter = n.iter, MIN_CHANGE = MIN_CHANGE)
  return(results)
}

check_pairs_type <- function(cell_type_profiles, bead, UMI_tot, score_mat, min_score, my_type, class_df, QL_score_cutoff, constrain) {
  candidates = rownames(score_mat)
  singlet_score = get_singlet_score(cell_type_profiles, bead, UMI_tot, my_type, constrain)
  all_pairs = T; all_pairs_class = !is.null(class_df)
  other_class = my_type #other types present from this class
  for(i in 1:(length(candidates)-1)) {
    type1 = candidates[i]
    for(j in (i+1):length(candidates)) {
      type2 = candidates[j]
      if(score_mat[i,j] < min_score + QL_score_cutoff) {
        if(type1 != my_type && type2 != my_type)
          all_pairs = F
        if(!is.null(class_df)) {
          first_class = class_df[my_type,"class"] == class_df[type1,"class"]
          second_class = class_df[my_type,"class"] == class_df[type2,"class"]
          if(!first_class && !second_class)
            all_pairs_class = F
          if(first_class && ! (type1 %in% other_class))
            other_class = c(other_class, type1)
          if(second_class && ! (type2 %in% other_class))
            other_class = c(other_class, type2)
        }
      }
    }
  }
  if(is.null(class_df))
    all_pairs_class = all_class
  if(all_pairs_class && !all_pairs && length(other_class) > 1) {
    for (type in other_class[2:length(other_class)])
      singlet_score = min(singlet_score, get_singlet_score(cell_type_profiles, bead, UMI_tot, type, constrain))
  }
  return(list(all_pairs = all_pairs, all_pairs_class = all_pairs_class, singlet_score = singlet_score))
}

#Decomposing a single bead via doublet search
process_bead_doublet <- function(cell_type_info, gene_list, UMI_tot, bead, class_df = NULL, constrain = T, verbose = F) {
  cell_type_profiles <- cell_type_info[[1]][gene_list,]
  cell_type_profiles = cell_type_profiles * UMI_tot
  cell_type_profiles = data.matrix(cell_type_profiles)
  QL_score_cutoff = 10; doublet_like_cutoff = 25
  results_all = decompose_full(cell_type_profiles, UMI_tot, bead, constrain = constrain, verbose = verbose)
  all_weights <- results_all$weights
  conv_all <- results_all$converged
  initial_weight_thresh = 0.01; cell_type_names = cell_type_info[[2]]
  candidates <- names(which(all_weights > initial_weight_thresh))
  if(length(candidates) == 1)
    if(candidates[1] == cell_type_info[[2]][1])
      candidates = c(candidates, cell_type_info[[2]][2])
    else
      candidates = c(candidates, cell_type_info[[2]][1])
  score_mat = Matrix(0, nrow = length(candidates), ncol = length(candidates))
  rownames(score_mat) = candidates; colnames(score_mat) = candidates
  min_score = 0
  first_type = NULL; second_type = NULL
  first_class = F; second_class = F #indicates whether the first (resp second) refers to a class rather than a type
  for(i in 1:(length(candidates)-1)) {
    type1 = candidates[i]
    for(j in (i+1):length(candidates)) {
      type2 = candidates[j]
      score = decompose_sparse(cell_type_profiles, UMI_tot, bead, type1, type2, score_mode = T, constrain = constrain, verbose = verbose)
      score_mat[i,j] = score; score_mat[j,i] = score
      if(is.null(second_type) || score < min_score) {
        first_type <- type1; second_type <- type2
        min_score = score
      }
    }
  }
  type1_pres = check_pairs_type(cell_type_profiles, bead, UMI_tot, score_mat, min_score, first_type, class_df, QL_score_cutoff, constrain)
  type2_pres = check_pairs_type(cell_type_profiles, bead, UMI_tot, score_mat, min_score, second_type, class_df, QL_score_cutoff, constrain)
  if(!type1_pres$all_pairs_class && !type2_pres$all_pairs_class) {
    spot_class <- "reject"
    singlet_score = min_score + 2 * doublet_like_cutoff #arbitrary
  }
  else if(type1_pres$all_pairs_class && !type2_pres$all_pairs_class) {
    first_class <- !type1_pres$all_pairs
    singlet_score = type1_pres$singlet_score
    spot_class = "doublet_uncertain"
  } else if(!type1_pres$all_pairs_class && type2_pres$all_pairs_class) {
    first_class <- !type2_pres$all_pairs
    singlet_score = type2_pres$singlet_score
    temp = first_type; first_type = second_type; second_type = temp
    spot_class = "doublet_uncertain"
  } else {
    spot_class = "doublet_certain"
    singlet_score = min(type1_pres$singlet_score, type2_pres$singlet_score)
    first_class <- !type1_pres$all_pairs; second_class <- !type2_pres$all_pairs
    if(type2_pres$singlet_score < type1_pres$singlet_score) {
      temp = first_type; first_type = second_type; second_type = temp
      first_class <- !type2_pres$all_pairs; second_class <- !type1_pres$all_pairs
    }
  }
  if(singlet_score - min_score < doublet_like_cutoff)
    spot_class = "singlet"
  doublet_results = decompose_sparse(cell_type_profiles, UMI_tot, bead, first_type, second_type, constrain = constrain)
  doublet_weights = doublet_results$weights; conv_doublet = doublet_results$converged
  spot_class <- factor(spot_class, c("reject", "singlet", "doublet_certain", "doublet_uncertain"))
  return(list(all_weights = all_weights, spot_class = spot_class, first_type = first_type, second_type = second_type,
              doublet_weights = doublet_weights, min_score = min_score, singlet_score = singlet_score,
              conv_all = conv_all, conv_doublet = conv_doublet, score_mat = score_mat,
              first_class = first_class, second_class = second_class))
}

get_prediction_sparse <- function(cell_type_profiles, UMI_tot, p, type1, type2) {
  cell_types = c(type1,type2)
  reg_data = cell_type_profiles[,cell_types]
  prediction = reg_data %*% c(p,1-p)
  return(prediction)
}

get_singlet_score <- function(cell_type_profiles, bead, UMI_tot, type, constrain) {
  if(!constrain)
    return(decompose_sparse(cell_type_profiles, UMI_tot, bead, type1=type, score_mode = T, constrain = constrain))
  dummy_type = cell_type_info[[2]][1]
  if(dummy_type == type)
    dummy_type = cell_type_info[[2]][2]
  prediction <- get_prediction_sparse(cell_type_profiles, UMI_tot, 1, type, dummy_type)
  log_l <- calc_log_l_par(gene_list, prediction, bead)
  return(log_l)
}

#decompose a doublet into two cells
decompose_doublet_fast <- function(bead, weights, gene_list, cell_type_info, type1, type2) {
  N_genes = length(gene_list)
  expect_1 = vector(mode="numeric",length = N_genes)
  expect_2 = vector(mode="numeric",length = N_genes)
  variance = vector(mode="numeric",length = N_genes)
  names(expect_1) = gene_list; names(expect_2) = gene_list; names(variance) = gene_list
  epsilon = 1e-10
  for(ind in which(bead > 0)) {
    gene = gene_list[ind]
    denom = weights[1] * cell_type_info[[1]][gene,type1] + weights[2] * cell_type_info[[1]][gene,type2] + epsilon
    posterior_1 = (weights[1] * cell_type_info[[1]][gene,type1] + epsilon / 2) / denom
    expect_1[[ind]] = posterior_1 * bead[gene]
    expect_2[[ind]] = bead[gene] - posterior_1 * bead[gene]
    variance[[ind]] = posterior_1 * bead[gene] * (1 - posterior_1)
  }
  return(list(expect_1 = expect_1, expect_2 = expect_2, variance = variance))
}
