test_single_beads <- function(puck, gene_list, cell_type_info, trust_model = FALSE, constrain = T, OLS = F) {
  cell_type_names = cell_type_info[[2]]; n_cell_types = cell_type_info[[3]]
  beads = t(as.matrix(puck@counts[gene_list,]))
  weights = decompose_batch(puck@nUMI, cell_type_info[[1]], beads, gene_list, constrain = constrain, OLS = OLS)
  pred_labels = unlist(lapply(weights,function(x) which.max(x$weights)))
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
#if proportions is null, does not renormalize cell type means
process_data <- function(puck, gene_list, cell_type_info, proportions = NULL, trust_model = FALSE, constrain = T, OLS = F) {
  cell_type_info_renorm = cell_type_info
  if(!is.null(proportions)) {
    cell_type_info_renorm[[1]] = get_norm_ref(puck, cell_type_info[[1]], gene_list, proportions)
  }
  test_results <- test_single_beads(puck, gene_list, cell_type_info_renorm, trust_model = trust_model, constrain = constrain, OLS = OLS)
  return(test_results)
}


#doublet decomposition for the whole puck.
process_beads_batch <- function(cell_type_info, gene_list, puck, class_df = NULL, constrain = T, doublet_mode = F) {
  beads = t(as.matrix(puck@counts[gene_list,]))
  out_file = "logs/process_beads_log.txt"
  if (file.exists(out_file))
    file.remove(out_file)
  numCores = parallel::detectCores(); MAX_CORES = 8
  if(parallel::detectCores() > MAX_CORES)
    numCores <- MAX_CORES
  cl <- parallel::makeCluster(numCores,outfile="") #makeForkCluster
  doParallel::registerDoParallel(cl)
  environ = c('process_bead','decompose_full','decompose_sparse','solveIRWLS.weights','solveOLS','solveWLS','Q_mat','X_vals','K_val', 'use_Q')
  results <- foreach::foreach(i = 1:(dim(beads)[1]), .packages = c("quadprog"), .export = environ) %dopar% {
    if(i %% 100 == 0)
      cat(paste0("Finished sample: ",i,"\n"), file=out_file, append=TRUE)
    assign("Q_mat",Q_mat, envir = globalenv()); assign("X_vals",X_vals, envir = globalenv())
    assign("K_val",K_val, envir = globalenv()); assign("use_Q",use_Q, envir = globalenv())
    if(doublet_mode)
      result = process_bead_doublet(cell_type_info, gene_list, puck@nUMI[i], beads[i,], class_df = class_df, constrain = constrain)
    else
      result = process_bead(cell_type_info, gene_list, puck@nUMI[i], beads[i,], class_df = class_df, constrain = constrain)
    result
  }
  parallel::stopCluster(cl)
  return(results)
}

#in parallel, does the (all cell type) decomposition of a batch of beads
decompose_batch <- function(nUMI, cell_type_means, beads, gene_list, constrain = T, OLS = F) {
  out_file = "logs/decompose_batch_log.txt"
  if (file.exists(out_file))
    file.remove(out_file)
  numCores = parallel::detectCores()
  if(parallel::detectCores() > 8)
    numCores <- 8
  cl <- parallel::makeCluster(numCores,outfile="") #makeForkCluster
  doParallel::registerDoParallel(cl)
  environ = c('decompose_full','solveIRWLS.weights',
              'solveOLS','solveWLS', 'Q_mat', 'K_val','X_vals','use_Q')
  #for(i in 1:100) {
  weights <- foreach::foreach(i = 1:(dim(beads)[1]), .packages = c("quadprog"), .export = environ) %dopar% {
    if(i %% 100 == 0)
      cat(paste0("Finished sample: ",i,"\n"), file=out_file, append=TRUE)
    assign("Q_mat",Q_mat, envir = globalenv()); assign("X_vals",X_vals, envir = globalenv())
    assign("K_val",K_val, envir = globalenv()); assign("use_Q",use_Q, envir = globalenv())
    decompose_full(cell_type_means, gene_list, nUMI[i], beads[i,], constrain = constrain, OLS = OLS)
  }
  parallel::stopCluster(cl)
  return(weights)
}


