downsample <- function(raw.data, cell.names, UMI = 500, replicates = 1) {
  datalist = list()
  labellist = list()
  for (i in 1:dim(raw.data)[2]) { #dim(reference@assays$RNA@counts)[2]) {
    datalist[[i]] <- Matrix::Matrix(rmultinom(replicates,UMI,raw.data[,i]),sparse=TRUE)
    labellist[[i]] <- rep(cell.names[i],replicates)
  }
  sub_beads = do.call(cbind, datalist)
  sub_labels = unlist(labellist)
  return (list(sub_beads,sub_labels))
}

#create down sampled data
#if each cell_type, then it takes n_samples from each cell type. Otherwise, it takes randomly from whole dataset
create_downsampled_data <- function(reference, refdir, cell_type_import = NULL,n_samples = 1000, each_cell_type = T,save.file = T) {
  if(!each_cell_type)
    index_keep = sample(which(reference@meta.data$liger_ident_coarse != ""),n_samples,replace=FALSE)
  else {
    if(is.null(cell_type_import))
      cell_type_import = levels(reference@meta.data$liger_ident_coarse)
    cell_type_import = cell_type_import[unlist(lapply(cell_type_import, function(x) nchar(x) > 0))]
    index_keep = c(); i = 1
    repeat{
      new_index = which(reference@meta.data$liger_ident_coarse == cell_type_import[i])
      new_samples = min(n_samples, length(new_index))
      index_keep = c(index_keep, sample(new_index,new_samples,replace=FALSE))
      if((i = i + 1) > length(cell_type_import))
        break
    }
  }
  reference@assays$RNA@counts = reference@assays$RNA@counts[,index_keep]
  reference@meta.data = reference@meta.data[index_keep,]
  reference@meta.data$liger_ident_coarse = droplevels(reference@meta.data$liger_ident_coarse)
  reference@assays$RNA@data <-matrix(0,nrow=2,ncol=2)
  reference@assays$RNA@scale.data <- matrix(0,nrow=2,ncol=2)
  if(save.file)
    saveRDS(reference, paste(refdir,"/scRefSubsampled", n_samples, ".RDS",sep=""))
  return(reference)
}

get_marker_genes <- function(reference, cell_type_means, gene_list, deve_thresh = 0.6, verbose = FALSE) {
  #compute percent deviance explained
  raw.data = reference@assays$RNA@counts[gene_list,]#[gene_list[1:100],] #c("Pcp2", "Malat1")
  nUMI = colSums(reference@assays$RNA@counts)
  mat = do.call("cbind",lapply(reference@meta.data$liger_ident_coarse, function (x) cell_type_means[gene_list,][[as.character(x)]]))
  nmat = rowMeans(mat)
  Tmat = t(nUMI)
  R_dev <- vector("list", length(gene_list))
  names(R_dev) <- gene_list
  for (i in 1:length(gene_list)) {
    null_D = sum((raw.data[i,] - nmat[i] * nUMI + raw.data[i,] * (log(nmat[i] * nUMI + 1e-9) - log(raw.data[i,] + 1e-9)))/nUMI)
    model_D = sum((raw.data[i,] - mat[i,] * nUMI + raw.data[i,] * (log(mat[i,] * nUMI + 1e-9) - log(raw.data[i,] + 1e-9)))/nUMI)
    #null_D = rowSums(sweep((raw.data - nmat*Tmat) + raw.data * (log(nmat*Tmat + 1e-9) - log(raw.data + 1e-9)), 2, nUMI, '/'))
    #model_D = rowSums(sweep((raw.data - mat*Tmat) + raw.data * (log(mat*Tmat + 1e-9) - log(raw.data + 1e-9)), 2, nUMI, '/'))
    R_dev[[i]] = 1 - model_D / null_D ### SOME NAN VALS
    if(verbose && (i %% 100 == 0))
      print(paste("Progress:",i))
  }
  keep_index = R_dev > deve_thresh
  return(gene_list[keep_index])  #R_dev[keep_index]
}

#if marker_provided is TRUE, then sets marker_list equal to gene_list
get_marker_data <- function(cell_type_names, reference, cell_type_means, gene_list, deve_thresh = 0.6, marker_plot = FALSE, marker_provided = FALSE) {
  if(marker_provided) {
    marker_means = cell_type_means[gene_list,]
    marker_norm = marker_means / rowSums(marker_means)
    marker_data = as.data.frame(cell_type_names[max.col(marker_means)])
    marker_data$entropy = rowSums(marker_norm * log(marker_norm + 1e-12))
    colnames(marker_data) = c("cell_type",'entropy')
    rownames(marker_data) = gene_list
    return(marker_data)
  } else {
    marker_list_info = get_marker_genes(reference, cell_type_means, gene_list, deve_thresh) # 0.4 is the best so far
    marker_list = marker_list_info[[1]]; marker_R_vals = marker_list_info[[2]]
    marker_data = as.data.frame(marker_R_vals)
    colnames(marker_data) = "R_val"
    marker_means = cell_type_means[marker_list,]
    marker_data$mean = log(rowMeans(marker_means)) ## symmetric for now
    marker_norm = marker_means / rowSums(marker_means)
    #marker_data$entropy = rowSums(marker_norm * log(marker_norm + 1e-12))
    #marker_data$across = marker_list %in% marker_pass
    marker_data$cell_type = cell_type_names[max.col(marker_means)]
    if(marker_plot)
      ggplot2::ggplot(marker_data, ggplot2::aes(x=mean, y=entropy, color=across)) + ggplot2::geom_point();
  }
}

remap_celltypes <- function(cell_dict_file, cell_ident) {
  cell_type_dict <- read.csv(file=cell_dict_file, header=TRUE, sep=",")
  rownames(cell_type_dict) = cell_type_dict[,'Cluster']
  true_type_names = lapply(cell_ident, function(x) cell_type_dict[as.character(x),"Name"])
  true_type_names = unlist(true_type_names)
}

#get the worst genes responsible for the misclassification
#raw data would be like cb.obj@assays$RNA@counts
get_worst_gene <- function(raw_data, cell_type_means, gene_list, desired_cell_type, alt_cell_type, cell_index) {
  smooth_cell_type_means = log(cell_type_means + 1e-9)
  diff_scores = smooth_cell_type_means[gene_list, desired_cell_type] - smooth_cell_type_means[gene_list,alt_cell_type]
  this_cell = raw.data[gene_list,cell_index]
  total_score = sum(diff_scores * this_cell)
  worst_gene = gene_list[which.min(diff_scores * this_cell)]
  return(worst_gene)
}



#plots the deviance of genes for predicting a bead as a particular cell type.
#example: plot_deviance_hist(bead_vec, "Purkinje", cell_type_means, marker_list, nUMI)
#note: bead_vec must be same length as gene list
plot_deviance_hist <- function(bead_vec,cell_type, cell_type_means, gene_list, nUMI, to_disp = 10, verbose = F) {
  if(length(bead_vec) != length(gene_list))
    stop("plot_deviance_hist: bead_vec and gene_list need to be the same length")
  reg_data = cell_type_means[gene_list,] * nUMI
  EPSILON = 1e-9
  deviance = (bead_vec - reg_data[[cell_type]]) + bead_vec * (log(reg_data[[cell_type]] + EPSILON) - log(bead_vec + EPSILON))
  dev_df = as.data.frame(deviance)
  rownames(dev_df) = gene_list
  dev_df = dev_df[order(deviance), , drop=FALSE]
  if(verbose)
    print(sum(dev_df))
  return(head(dev_df,to_disp))
}

#gets a gene list which is determined by genes in the cell_type_means matrix, filtered by
#having a mean expression large enough, being present in both dropviz and slideseq puck.
get_gene_list <- function(cell_type_means, puck, cutoff_val = 1/20000) {
  bulk_vec = rowSums(puck@counts)
  gene_list = names(which(rowMeans(cell_type_means) > cutoff_val))
  gene_list = intersect(gene_list,names(bulk_vec))
  gene_list = gene_list[bulk_vec[gene_list] > 0]
  return(gene_list)
}

#finds DE genes
#Genes must be observed a minimum of MIN_OBS times to mitigate sampling noise in the
#Platform effect estimation
get_de_genes <- function(cell_type_info, puck, fc_thresh = 1.25, expr_thresh = .00015, MIN_OBS = 3) {
  total_gene_list = c()
  epsilon = 1e-9
  bulk_vec = rowSums(puck@counts)
  gene_list = rownames(cell_type_info[[1]])
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

#given a list of cell labels, computes the (smoothed) proportion with each label.
#if not null nUMI, then it does it weighted by nUMI
#if constrain, sum to 1, otherwise sum to orginial sum
compute_proportions <- function(cell_labels, cell_type_names, nUMI = NULL, smoothing_par = 1/1000) {
  if(is.null(nUMI))
    true_counts = as.data.frame(table(cell_labels))
  else {
    true_counts = aggregate(puck@nUMI, by=list(puck@cell_labels),FUN=sum)
  }
  rownames(true_counts) = true_counts[,1]
  true_counts[,1] <- NULL
  true_counts = true_counts[cell_type_names,1]
  true_counts[is.na(true_counts)] <- 0
  names(true_counts) = cell_type_names
  true_counts = true_counts + sum(true_counts)*smoothing_par #smoothing
  proportions = true_counts / sum(true_counts)
}

#if constrain, sum to 1, otherwise sum to orginial sum
smooth_proportions <- function(proportions, smoothing_par = 1/300, constrain = T) {
  target_sum = 1
  if(!constrain)
    target_sum = sum(proportions)
  proportions = proportions + smoothing_par #smoothing
  proportions = target_sum * proportions / sum(proportions)
}

#min weight to be considered a singlet as a function of nUMI
UMI_cutoff <- function(nUMI) {
  return (pmax(0.25, 2 - log(nUMI,2) / 5))
}

#get correlation between weight and markers
get_corr <- function(cell_type_info, norm_marker_scores, weights) {
  correlations = numeric(cell_type_info[[3]])
  names(correlations) = cell_type_info[[2]]
  for (cell_type in cell_type_info[[2]])
    correlations[cell_type] = cor(norm_marker_scores[,cell_type], weights[,cell_type])
  write.csv(correlations, file.path(resultsdir,"correlation_markers.csv"))
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

#if test_reference is not null, it will convert this to a slideseq object rather than read in one
#if puck_file is not null, then reads in puck from this file
#if load_info_renorm, loads cell type info from MetaData/cell_type_info_renorm.RDS. Takes gene_list to be rownames of cell_type_info (renorm)
#get_proportions -> calculates cell type info renorm
init_RCTD <- function(gene_list_reg = T, get_proportions = F, test_reference = NULL, puck_file = NULL, MIN_OBS = 3, load_info_renorm = F) {
  print("init_RCTD: begin")
  config <- config::get(file = "conf/test.yml", use_parent = FALSE)
  config_data <- config::get(file = "conf/dataset.yml", use_parent = FALSE)
  slideseqdir <- file.path("Data/Slideseq",config_data$slideseqfolder)
  resultsdir = file.path(slideseqdir,"results")
  if(!dir.exists(resultsdir))
    dir.create(resultsdir)
  bulkdir <- paste(slideseqdir,"results/Bulk",sep="/")
  if(!dir.exists(bulkdir))
    dir.create(bulkdir)
  if(load_info_renorm) {
    cell_type_info <- readRDS(file.path(slideseqdir, "MetaData/cell_type_info_renorm.RDS"))
    reference <- NULL; refdir <- NULL
  } else {
    refdir <- file.path("Data/Reference",config_data$reffolder)
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
    proportions <- read.csv(file.path(bulkdir,"weights.csv"))$Weight
    names(proportions) = cell_type_info[[2]]
    proportions <- proportions / sum(proportions)
    print("init_RCTD: estimated bulk composition: ")
    print(proportions)
  }
  if(is.null(test_reference)) {
    if(is.null(puck_file))
      puck = readRDS(file.path(slideseqdir, config_data$puckrds))
    else
      puck = readRDS(file.path(slideseqdir, puck_file))
  }
  else
    puck <- seurat.to.slideseq(test_reference, cell_type_info)
  puck = restrict_counts(puck, rownames(puck@counts), UMI_thresh = config$UMI_min)
  if(load_info_renorm)
    gene_list = rownames(cell_type_info[[1]])
  else {
    if(gene_list_reg)
      gene_list = get_de_genes(cell_type_info, puck, fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg, MIN_OBS = MIN_OBS)
    else
      gene_list = get_de_genes(cell_type_info, puck, fc_thresh = config$fc_cutoff, expr_thresh = config$gene_cutoff, MIN_OBS = MIN_OBS)
  }
  print(paste("init_RCTD: number of genes used:", length(gene_list)))
  puck = restrict_counts(puck, gene_list, UMI_thresh = config$UMI_min)
  puck = restrict_puck(puck, colnames(puck@counts))
  print(paste("init_RCTD: number of spots used in test data passing UMI threshold:", dim(puck@counts)[2]))
  print("init_RCTD: end")
  return(list(refdir = refdir, slideseqdir = slideseqdir, bulkdir = bulkdir, reference = reference,
              proportions = proportions, gene_list = gene_list, puck = puck, cell_type_info = cell_type_info,
              config = config, n_puck_folds = config_data$n_puck_folds))
}

get_class_df <- function(cell_type_names) {
  class_df = data.frame(cell_type_names, row.names = cell_type_names)
  colnames(class_df)[1] = "class"
  class_df["Bergmann","class"] = "Astrocytes"
  class_df["Fibroblast","class"] = "Endothelial"
  class_df["MLI2","class"] = "MLI1"
  class_df["Macrophages","class"] = "Microglia"
  class_df["Polydendrocytes","class"] = "Oligodendrocytes"
  return(class_df)
}

chooseSigma <- function(prediction, counts, resultsdir, sigma_init = 1, N_epoch = 15, folder_id = "") {
  X = as.vector(prediction)
  X = pmax(X, 1e-4)
  Y = as.vector(counts)
  num_sample = min(500000, length(X)) #300000
  big_params = F
  use_ind = sample(1:length(X), num_sample)
  X = X[use_ind]
  Y = Y[use_ind]
  sigma = sigma_init
  alpha_init = 0.0001; batch = 25
  X_batch = list()
  Y_batch = list()
  for(k in as.numeric(names(table(Y)))) {
    X_vals <- X[Y==k]; N_X = length(X_vals)
    for(b in 1:ceiling(N_X/batch)) {
      X_ind = (batch*(b-1) + 1):min((batch*b),N_X)
      curr_X = X_vals[X_ind]
      X_batch[[length(X_batch) + 1]] <- curr_X
      Y_batch[[length(Y_batch) + 1]] <- k
    }
  }
  ordering <- sample(1:length(X_batch))
  sigma_vals = list(sigma)
  loss_vals = list()
  for(j in 1:N_epoch) {
    alpha = alpha_init / j
    total_loss = 0
    for(i in ordering) {
      Q = get_Q(X_batch[[i]], Y_batch[[i]], sigma, big_params = big_params)
      Q_d = get_Q_d(X_batch[[i]], Y_batch[[i]], sigma, big_params = big_params)
      sigma = sigma + sum(Q_d/Q)*alpha
      if(i%%100 == 0)
        print(sigma)
      total_loss = total_loss + sum(log(Q))
      sigma_vals[[length(sigma_vals) + 1]] <- sigma
    }
    print(total_loss)
    loss_vals[[length(loss_vals) + 1]] <- total_loss
  }
  sigresdir = file.path(resultsdir, paste0("sigma",folder_id))
  if(!dir.exists(sigresdir))
    dir.create(sigresdir)

  saveRDS(sigma_vals, file.path(sigresdir,"sigma_vals.RDS"))
  saveRDS(loss_vals, file.path(sigresdir,"loss_vals.RDS"))
  saveRDS(sigma, file.path(sigresdir,"sigma.RDS"))

  pdf(file.path(sigresdir,"sigma_trace.pdf"))
  plot(unlist(sigma_vals),type = 'n')
  lines(unlist(sigma_vals))
  dev.off()

  pdf(file.path(sigresdir,"loss_trace.pdf"))
  plot(unlist(loss_vals),type = 'n')
  lines(unlist(loss_vals))
  dev.off()
  return(sigma)
}
