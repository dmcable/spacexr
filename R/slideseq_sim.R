#given a cell and number of UMI_sample, samples without replacement from the cell
sub_sample_cell <- function(gene_list, raw.data, cell_index, UMI_sample) {
  sub_sample = sample(rep(rownames(raw.data),raw.data[,cell_index]),UMI_sample,replace=FALSE)
  sub_sample = table(sub_sample)[as.character(gene_list)]
  names(sub_sample) = gene_list
  sub_sample[is.na(sub_sample)] <- 0
  return(sub_sample)
}

#creates a simulated bead as a mixture of two downsampled cells.
#nUMI is vector of total UMI in the reference
#UMI1 is how much to sample from first cell, and UMI2 is how much to sample from second cell
bead_mix <- function(test_ref, gene_list, UMI1, UMI2, type1, type2) {
  nUMI = test_ref@meta.data$nUMI
  firstInd = sample(intersect(which(test_ref@meta.data$liger_ident_coarse == type1) , which(nUMI > 1000)),1)
  secondInd = sample(intersect(which(test_ref@meta.data$liger_ident_coarse == type2) , which(nUMI > 1000)),1)
  bead = sub_sample_cell(gene_list, test_ref@assays$RNA@counts, firstInd, UMI1)
  bead = bead + sub_sample_cell(gene_list, test_ref@assays$RNA@counts, secondInd, UMI2)
  return(data.matrix((as.matrix(bead))))
}

bead_singlet <- function(test_ref, gene_list, UMI, cell_type) {
  nUMI = test_ref@meta.data$nUMI
  firstInd = sample(intersect(which(test_ref@meta.data$liger_ident_coarse == cell_type) , which(nUMI > 1000)),1)
  bead = sub_sample_cell(gene_list, test_ref@assays$RNA@counts, firstInd, UMI)
  return(data.matrix((as.matrix(bead))))
}

#decompose with just two cell types
#if score_mode, then returns the objective function score
#if denoise, then it fits a "noise" dimension as the mean of all the data
decompose_sparse <- function(cell_type_means, gene_list, nUMI, bead, type1=NULL, type2=NULL, score_mode = FALSE, plot = F, custom_list = NULL, verbose=F, constrain = T) {
  if(is.null(custom_list))
    cell_types = c(type1,type2)
  else
    cell_types = custom_list
  reg_data = cell_type_means[gene_list,] * nUMI
  reg_data = data.matrix(reg_data)
  reg_data = data.matrix(reg_data[,cell_types])
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

#get weights and doublet
decompose_sparse_report <- function(cell_type_info, gene_list, UMI_tot, bead, type1, type2, constrain = F, verbose = F) {
  results <- decompose_sparse(cell_type_info[[1]], gene_list, UMI_tot, bead, type1, type2, score_mode = F, constrain = constrain, verbose = verbose)
  decompose_results <- decompose_doublet(bead, results$weights, gene_list, cell_type_info, type1, type2)
  return(list(weights = results$weights, decompose_results = decompose_results))
}

get_likelihood <- function(gene_list, prediction, bead) {
  delta = 1e-5
  total_score=0
  for(gene in gene_list) {
    Y = bead[gene]
    #x = round(prediction[gene,]/delta) + 1
    total_score = total_score + get_QL(Y, prediction[gene,]) #J_mat[Y+1,x]
  }
  return(-total_score)
}

decompose_sparse_score <- function(cell_type_means, gene_list, nUMI, bead, type1=NULL, type2=NULL) {
  cell_types = c(type1,type2)
  reg_data = cell_type_means[gene_list,] * nUMI
  reg_data = data.matrix(reg_data)
  reg_data = reg_data[,cell_types]
  M = 100
  best_score = Inf
  for(p in (0:M)/M) {
    weights=c(p, 1-p)
    prediction = reg_data %*% weights
    scaled_residual = (prediction - bead)/(sqrt(get_V(prediction)))
    my_score = sum(abs(scaled_residual))
    if(my_score < best_score)
      best_score = my_score
  }
  return(best_score)
}

#decompose with all cell types
decompose_full <- function(cell_type_means, gene_list, nUMI, bead, constrain = TRUE, OLS = FALSE, verbose = F, n.iter = 50, MIN_CHANGE = 0.001) {
  reg_data = cell_type_means[gene_list,] * nUMI
  reg_data = data.matrix(reg_data)
  results = solveIRWLS.weights(reg_data,bead,nUMI,OLS = OLS, constrain = constrain, verbose = verbose, n.iter = n.iter, MIN_CHANGE = MIN_CHANGE)
  results$weights <- results$weights
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

check_pairs_type <- function(cell_type_info, gene_list, bead, UMI_tot, score_mat, min_score, my_type, class_df, QL_score_cutoff, constrain) {
  candidates = rownames(score_mat)
  singlet_score = get_singlet_score(cell_type_info, gene_list, bead, UMI_tot, my_type, constrain)
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
      singlet_score = min(singlet_score, get_singlet_score(cell_type_info, gene_list, bead, UMI_tot, type, constrain))
  }
  return(list(all_pairs = all_pairs, all_pairs_class = all_pairs_class, singlet_score = singlet_score))
}

#Decomposing a single bead via doublet search
process_bead_doublet <- function(cell_type_info, gene_list, UMI_tot, bead, class_df = NULL, constrain = T, verbose = F) {
  QL_score_cutoff = 10; doublet_like_cutoff = 25
  results_all = decompose_full(cell_type_info[[1]], gene_list, UMI_tot, bead, constrain = constrain, verbose = verbose)
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
      score = decompose_sparse(cell_type_info[[1]], gene_list, UMI_tot, bead, type1, type2, score_mode = T, constrain = constrain, verbose = verbose)
      score_mat[i,j] = score; score_mat[j,i] = score
      if(is.null(second_type) || score < min_score) {
        first_type <- type1; second_type <- type2
        min_score = score
      }
    }
  }
  type1_pres = check_pairs_type(cell_type_info, gene_list, bead, UMI_tot, score_mat, min_score, first_type, class_df, QL_score_cutoff, constrain)
  type2_pres = check_pairs_type(cell_type_info, gene_list, bead, UMI_tot, score_mat, min_score, second_type, class_df, QL_score_cutoff, constrain)
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
  doublet_results = decompose_sparse(cell_type_info[[1]], gene_list, UMI_tot, bead, first_type, second_type, constrain = constrain)
  doublet_weights = doublet_results$weights; conv_doublet = doublet_results$converged
  spot_class <- factor(spot_class, c("reject", "singlet", "doublet_certain", "doublet_uncertain"))
  return(list(all_weights = all_weights, spot_class = spot_class, first_type = first_type, second_type = second_type,
              doublet_weights = doublet_weights, min_score = min_score, singlet_score = singlet_score,
              conv_all = conv_all, conv_doublet = conv_doublet, score_mat = score_mat,
              first_class = first_class, second_class = second_class))
}

#main function for decomposing a single bead
process_bead <- function(cell_type_info, gene_list, UMI_tot, bead, class_df = NULL, constrain = T) {
  QL_score_cutoff = 10; doublet_like_cutoff = 25
  results_all = decompose_full(cell_type_info[[1]], gene_list, UMI_tot, bead, constrain = constrain)
  all_weights <- results_all$weights
  conv_all <- results_all$converged
  first_type = names(which.max(all_weights))
  min_score = 0
  second_score = 0
  other_type = NULL
  second_type = NULL
  for (type in cell_type_info[[2]])
    if(type != first_type) {
      score = decompose_sparse(cell_type_info[[1]], gene_list, UMI_tot, bead, type, first_type, score_mode = T, constrain = constrain)
      if(is.null(second_type) || score < min_score) {
        other_type <- second_type
        second_score <- min_score
        min_score = score
        second_type = type
      } else if(is.null(other_type) || score < second_score) {
        second_score <- score
        other_type <- type
      }
    }
  doublet_results = decompose_sparse(cell_type_info[[1]], gene_list, UMI_tot, bead, first_type, second_type, constrain = constrain)
  doublet_weights = doublet_results$weights; conv_doublet = doublet_results$converged
  dlik <- delta_likelihood(first_type, second_type, doublet_weights[1], bead, cell_type_info[[1]], gene_list, UMI_tot)
  class_other = F
  if(!is.null(class_df))
    if(class_df[second_type,"class"] == class_df[other_type,"class"])
      class_other = T
  if(max(all_weights) < UMI_cutoff(UMI_tot))
    spot_class <- "reject"
  else if(dlik < doublet_like_cutoff)
    spot_class <- "singlet"
  else if(second_score - min_score > QL_score_cutoff)
    spot_class <- "doublet_certain"
  else if(class_other)
    spot_class <- "doublet_certain_class"
  else
    spot_class <- "doublet_uncertain"
  spot_class <- factor(spot_class, c("reject", "singlet", "doublet_certain", "doublet_certain_class", "doublet_uncertain"))
  return(list(all_weights = all_weights, spot_class = spot_class, first_type = first_type, second_type = second_type,
              doublet_weights = doublet_weights, min_score = min_score, second_score = second_score,
              other_type = other_type, conv_all = conv_all, conv_doublet = conv_doublet))
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

#doublet decomposition for the whole puck.
process_beads_sparse <- function(cell_type_info, gene_list, puck, meta_df, constrain = F) {
  beads = t(as.matrix(puck@counts[gene_list,]))
  out_file = "logs/process_beads_log.txt"
  if (file.exists(out_file))
    file.remove(out_file)
  numCores = parallel::detectCores(); MAX_CORES = 8
  if(parallel::detectCores() > MAX_CORES)
    numCores <- MAX_CORES
  cl <- parallel::makeCluster(numCores,outfile="") #makeForkCluster
  doParallel::registerDoParallel(cl)
  environ = c('process_bead','decompose_full','decompose_sparse','solveIRWLS.weights','solveOLS','solveWLS','Q_mat','X_vals','K_val', 'use_Q', 'decompose_sparse_report')
  results <- foreach::foreach(i = 1:(dim(beads)[1]), .packages = c("quadprog"), .export = environ) %dopar% {
    if(i %% 100 == 0)
      cat(paste0("Finished sample: ",i,"\n"), file=out_file, append=TRUE)
    assign("Q_mat",Q_mat, envir = globalenv()); assign("X_vals",X_vals, envir = globalenv())
    assign("K_val",K_val, envir = globalenv()); assign("use_Q",use_Q, envir = globalenv())
    result = decompose_sparse_report(cell_type_info, gene_list, puck@nUMI[i], beads[i,], meta_df$first_type[i], meta_df$second_type[i], constrain = constrain)
    result
  }
  parallel::stopCluster(cl)
  return(results)
}


#run the weight recovery test. creates random mixtures between two cell types and tests the model's ability to
#recover the weights
#conditions should divide into UMI_tot
#trials is number of trials per condition
#there is actually 1 more condition than given by the conditions variable (the 0 condition)
weight_recovery_test <- function(test_ref, gene_list, cell_type_means, type1, type2, UMI_tot = 1000, conditions = 10, trials = 30, save.file = F, mydir = "") {
  UMI_step = round(UMI_tot / conditions)
  UMI_tot = UMI_step * conditions
  type1_avg = vector(mode="numeric",length = conditions + 1)
  type1_moment2 = vector(mode="numeric",length = conditions + 1)
  UMI1_vec = 0:conditions*UMI_step
  for (prop_ind in 1:(conditions + 1)) {
    UMI1 = UMI1_vec[prop_ind]
    UMI2 = UMI_tot - UMI1
    for (t in 1:trials) {
      bead= bead_mix(test_ref, gene_list, UMI1, UMI2, type1, type2)
      weights = decompose_sparse(cell_type_means, gene_list, UMI_tot, bead, type1, type2, constrain = F)$weights
      type1_avg[prop_ind] = type1_avg[prop_ind] + weights[1]
      type1_moment2[prop_ind] = type1_moment2[prop_ind] + (weights[1]^2)
    }
  }
  type1_avg = type1_avg / trials
  type1_moment2 = type1_moment2 / trials
  st_dev = sqrt(type1_moment2 - type1_avg^2 + 1e-12) #avoid sqrt of negative
  proportion = UMI1_vec / UMI_tot
  const_stderr = 1.96
  plot_df <- data.frame(proportion,type1_avg,const_stderr*st_dev)
  p<- ggplot2::ggplot(plot_df, ggplot2::aes(x=proportion, y=type1_avg, colour = "type1_avg")) +
    ggplot2::geom_line() +
    ggplot2::geom_point()+
    ggplot2::geom_line(ggplot2::aes(y=proportion,colour = "proportion")) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=type1_avg-st_dev, ymax=type1_avg+st_dev), width=.05,
                  position=ggplot2::position_dodge(0.05))
  print(p)
  if(save.file) {
    out_dir = paste(mydir,"Output/CompositionRecovery",sep="/")
    saveRDS(plot_df,file = paste(out_dir,"weight_recovery_df.RDS",sep="/"))
  }
}

doublet_accuracy_test <- function(test_ref, gene_list, cell_type_info, type1, type2, UMI_tot = 1000, conditions = 10, trials = 30, save.file = F, mydir = "") {
  UMI_step = round(UMI_tot / conditions)
  UMI_tot = UMI_step * conditions
  UMI1_vec = 0:conditions*UMI_step
  success = vector(mode="numeric",length = conditions + 1)
  first_failed_vec = vector(mode="numeric",length = conditions + 1)
  for (prop_ind in 1:(conditions+1)) {
    UMI1 = UMI1_vec[prop_ind]
    UMI2 = UMI_tot - UMI1
    for (t in 1:trials) {
      bead = bead_mix(test_ref, gene_list, UMI1, UMI2, type1, type2)
      weights = process_bead(cell_type_info, gene_list, UMI_tot, bead, constrain = F)
      max_name = names(which.max(weights))
      if(max_name != type1 && max_name != type2)
        first_failed_vec[prop_ind] = first_failed_vec[prop_ind] + 1
      else {
        min_name = names(which.min(weights))
        if(min_name == type1 || min_name == type2)
          success[prop_ind] = success[prop_ind] + 1
      }
    }
  }
  plot(UMI1_vec/UMI_tot,success/trials)
  if(save.file) {
    out_dir = paste(mydir,"Output/CompositionRecovery",sep="/")
    saveRDS(success,file = paste(out_dir,"success.RDS",sep="/"))
    saveRDS(first_failed_vec,file = paste(out_dir,"first_failed.RDS",sep="/"))
  }
}

#creates a plot (for the simulated data) of marker_score vs recovered weight
create_marker_plot <- function(test_ref,type1,type2, cell_type_means, marker_data, gene_list, UMI_tot = 1000, trials = 100) {
  weight_vec = vector(mode="numeric",length = trials)
  marker_vec = vector(mode="numeric",length = trials)
  for (ind in 1:trials) {
    UMI_first = sample(0:UMI_tot,1)
    bead = bead_mix(test_ref, union(rownames(marker_data),gene_list), UMI_first, UMI_tot - UMI_first, type1, type2)
    weight_vec[ind] = decompose_sparse(cell_type_means_renorm[gene_list,], gene_list, UMI_tot, bead[gene_list,], type1, type2)[1]
    marker_vec[ind] = get_marker_score(marker_data,cell_type_means, bead[marker_list,],type1,type2,UMI_tot)
  }
  plot(marker_vec,weight_vec)
}

#gets the score of what percent marker character is the bead of type1 vs type2 marker genes
get_marker_score <- function(marker_data, cell_type_means, bead, type1, type2, UMI_tot) {
  type1_score = get_marker_score_type(marker_data, bead, UMI_tot, type1, score_threshold = 10)
  type2_score = get_marker_score_type(marker_data, bead, UMI_tot, type2, score_threshold = 10)
  norm_score = type1_score / (type1_score + type2_score)
  return(norm_score)
}

#returns the marker score for cell_type
#gene_means is the mean of the gene in the test dataset or the reference dataset
get_marker_score_type <- function(marker_data, bead, UMI_tot, cell_type, gene_means,score_threshold = 10, score_mode = T) {
  gene_list = rownames(marker_data)[marker_data$cell_type==cell_type]
  mark_genes = unlist(lapply(bead[gene_list]/(gene_means[gene_list]*UMI_tot), function(x) min(x,score_threshold)))
  if(score_mode)
    return(mean(mark_genes))
  else
    return(tail(mark_genes[order(mark_genes)],10))
}

#get marker scores for all cell types
get_marker_scores <- function(marker_data, puck, cell_type_info, score_threshold = 10) {
  gene_means = rowMeans(cell_type_info[[1]][rownames(marker_data),])
  score_df <- Matrix(0, nrow = dim(puck@counts)[2], ncol = cell_type_info[[3]])
  rownames(score_df) = colnames(puck@counts); colnames(score_df) = cell_type_info[[2]]
  for(cell_type in cell_type_info[[2]]) {
    scores = numeric(dim(puck@counts)[2])
    for (i in 1:dim(puck@counts)[2]) {
      scores[i] <- get_marker_score_type(marker_data, puck@counts[,i], puck@nUMI[i], cell_type, gene_means, score_threshold)
    }
    score_df[,cell_type] = scores
  }
  return(score_df)
}


get_prediction_sparse <- function(cell_type_means, gene_list, UMI_tot, p, type1, type2) {
  cell_types = c(type1,type2)
  reg_data = cell_type_means[gene_list,] * 1000
  reg_data = data.matrix(reg_data)
  reg_data = reg_data[,cell_types]
  prediction = reg_data %*% c(p,1-p)
  return(prediction)
}

delta_likelihood <- function(type1, type2, p, bead, cell_type_means, gene_list, UMI_tot) {
  prediction <- get_prediction_sparse(cell_type_means, gene_list, UMI_tot, 1, type1, type2)
  log_l <- calc_log_l_par(gene_list, prediction, bead)
  prediction_doub <- get_prediction_sparse(cell_type_means, gene_list, UMI_tot, p, type1, type2)
  log_l_doub <- calc_log_l_par(gene_list, prediction_doub, bead)
  return(-log_l_doub + log_l)
}

get_singlet_score <- function(cell_type_info, gene_list, bead, UMI_tot, type, constrain) {
  if(!constrain)
    return(decompose_sparse(cell_type_info[[1]], gene_list, UMI_tot, bead, type1=type, score_mode = T, constrain = constrain))
  dummy_type = cell_type_info[[2]][1]
  if(dummy_type == type)
    dummy_type = cell_type_info[[2]][2]
  prediction <- get_prediction_sparse(cell_type_info[[1]], gene_list, UMI_tot, 1, type, dummy_type)
  log_l <- calc_log_l_par(gene_list, prediction, bead)
  return(log_l)
}

#decompose a doublet into two cells
decompose_doublet <- function(bead, weights, gene_list, cell_type_info, type1, type2) {
  N_genes = length(gene_list)
  expect_1 = vector(mode="numeric",length = N_genes)
  expect_2 = vector(mode="numeric",length = N_genes)
  variance = vector(mode="numeric",length = N_genes)
  names(expect_1) = gene_list; names(expect_2) = gene_list; names(variance) = gene_list
  epsilon = 1e-10
  for(ind in 1:N_genes) {
    gene = gene_list[ind]
    denom = weights[1] * cell_type_info[[1]][gene,type1] + weights[2] * cell_type_info[[1]][gene,type2] + epsilon
    posterior_1 = (weights[1] * cell_type_info[[1]][gene,type1] + epsilon / 2) / denom
    expect_1[[ind]] = posterior_1 * bead[gene]
    expect_2[[ind]] = bead[gene] - posterior_1 * bead[gene]
    variance[[ind]] = posterior_1 * bead[gene] * (1 - posterior_1)
  }
  return(list(expect_1 = expect_1, expect_2 = expect_2, variance = variance))
}
