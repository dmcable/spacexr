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
decompose_sparse <- function(cell_type_means, gene_list, nUMI, bead, type1=NULL, type2=NULL, score_mode = FALSE, plot = F, denoise = F, custom_list = NULL, verbose=F, constrain = T) {
  if(is.null(custom_list))
    cell_types = c(type1,type2)
  else
    cell_types = custom_list
  reg_data = cell_type_means[gene_list,] * nUMI
  reg_data = data.matrix(reg_data)
  if(denoise) {
    reg_data = cbind(reg_data[,cell_types],rowMeans(reg_data))
    colnames(reg_data)[length(cell_types) + 1] = "Noise"
  } else
    reg_data = reg_data[,cell_types]
  weights= solveIRWLS.weights(reg_data,bead,nUMI,OLS = FALSE, constrain = constrain, verbose = verbose)
  weights=weights / sum(weights)
  if(! score_mode)
    return(weights)
  else {
    prediction = reg_data %*% weights
    j<-19; threshold<-nUMI/2^(j-1)
    prediction[which(prediction < threshold)] = threshold
    scaled_residual = (reg_data %*% weights - bead)/(sqrt(prediction))
    my_score = sum(abs(scaled_residual))
    if(plot)
      hist(abs(scaled_residual)[abs(scaled_residual) > 5], breaks=60)
    return(my_score)
  }
}

#decompose
decompose <- function(cell_type_means, gene_list, nUMI, bead, constrain = TRUE, OLS = FALSE, verbose = F) {
  reg_data = cell_type_means[gene_list,] * nUMI
  reg_data = data.matrix(reg_data)
  weights= solveIRWLS.weights(reg_data,bead,nUMI,OLS = OLS, constrain = constrain, verbose = verbose)
  return(weights)
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
  environ = c('decompose','solveIRWLS.weights',
              'solveOLS','solveWLS')
  weights <- foreach::foreach(i = 1:(dim(beads)[1]), .packages = c("quadprog"), .export = environ) %dopar% {
    if(i %% 100 == 0)
      cat(paste0("Finished sample: ",i,"\n"), file=out_file, append=TRUE)
    decompose(cell_type_means, gene_list, nUMI[i], beads[i,], constrain = constrain, OLS = OLS)
  }
  parallel::stopCluster(cl)
  return(weights)
}

#main function for decomposing a single bead
process_bead <- function(cell_type_info, gene_list, UMI_tot, bead) {
  weights = decompose(cell_type_info[[1]], gene_list, UMI_tot, bead)
  first_type = names(which.max(weights))
  min_score = 0
  second_type = NULL
  for (type in cell_type_info[[2]])
    if(type != first_type) {
      score = decompose_sparse(cell_type_info[[1]], gene_list, UMI_tot, bead, type, first_type, score_mode = T)
      if(is.null(second_type) || score < min_score) {
        min_score = score
        second_type = type
      }
    }
  weights = decompose_sparse(cell_type_info[[1]], gene_list, UMI_tot, bead, first_type, second_type)
  return(weights)
}

#doublet decomposition for the whole puck.
process_beads_batch <- function(cell_type_info, gene_list, puck, constrain = T) {
  beads = t(as.matrix(puck@counts[gene_list,]))
  out_file = "logs/decompose_batch_log.txt"
  if (file.exists(out_file))
    file.remove(out_file)
  numCores = parallel::detectCores()
  if(parallel::detectCores() > 8)
    numCores <- 8
  cl <- parallel::makeCluster(numCores,outfile="") #makeForkCluster
  doParallel::registerDoParallel(cl)
  environ = c('process_bead','decompose','decompose_sparse','solveIRWLS.weights',
              'solveOLS','solveWLS')
  weights <- foreach::foreach(i = 1:(dim(beads)[1]), .packages = c("quadprog"), .export = environ) %dopar% {
    if(i %% 100 == 0)
      cat(paste0("Finished sample: ",i,"\n"), file=out_file, append=TRUE)
    process_bead(cell_type_info, gene_list, puck@nUMI[i], beads[i,])
  }
  parallel::stopCluster(cl)
  return(weights)
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
      weights = decompose_sparse(cell_type_means, gene_list, UMI_tot, bead, type1, type2, constrain = F)
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
      weights = process_bead(cell_type_info, gene_list, UMI_tot, bead)
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

#get marker scores for some cell type
get_marker_scores <- function(marker_data, puck, cell_type, cell_type_means, score_threshold = 10) {
  gene_list = rownames(marker_data)[marker_data$cell_type==cell_type]
  gene_means = rowMeans(cell_type_means[rownames(marker_data),])
  scores = numeric(dim(puck@counts)[2])
  for (i in 1:dim(puck@counts)[2]) {
    scores[i] <- get_marker_score_type(marker_data[gene_list,], puck@counts[gene_list,i], puck@nUMI[i], cell_type, gene_means, score_threshold)
  }
  names(scores) = colnames(puck@counts)
  return(scores)
}

#marker scores for all cell types
get_norm_marker_scores <- function(marker_data, puck, cell_type, cell_type_means, score_threshold = 10) {
  gene_list = rownames(marker_data)[marker_data$cell_type==cell_type]
  gene_means = rowMeans(cell_type_means[rownames(marker_data),])
  scores = numeric(dim(puck@counts)[2])
  for (i in 1:dim(puck@counts)[2]) {
    scores[i] <- get_marker_score_type(marker_data[gene_list,], puck@counts[gene_list,i], puck@nUMI[i], cell_type, gene_means, score_threshold)
  }
  names(scores) = colnames(puck@counts)
  return(scores)
}

#TODO: remember to take into account the mean of the new dataset
#simple classifier based on the marker_genes
#cell_type_means is from the reference
classify_cell_marker <- function(marker_data, bead, UMI_tot, cell_type_means, score_threshold = 10) {
  max_score = 0
  max_type = NULL
  gene_means = rowMeans(cell_type_means[rownames(marker_data),])
  for(cell_type in cell_type_names) {
    cell_score = get_marker_score_type(marker_data, bead, UMI_tot, cell_type, gene_means, score_threshold = score_threshold)
    if(is.null(max_type) || cell_score > max_score) {
      max_score = cell_score
      max_type = cell_type
    }
  }
  return(max_type)
}

#does the (all cell type) classification of a batch of beads
classify_marker_batch <- function(nUMI, cell_type_means, beads, marker_data) {
  types <- foreach(i = 1:(dim(beads)[1])) %do% {
    classify_cell_marker(marker_data, beads[i,], nUMI[i], cell_type_means)
  }
  return(types)
}
