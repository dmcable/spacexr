make_heat_map <- function(cell_type_means, gene_list) {
  heatmap(as.matrix(cell_type_means[gene_list,]))
}

#plotting where one cell type is bigger than others
plot_cell_types_spec <- function(puck, barcodes) {
  my_table = puck@coords[barcodes,]
  my_table$class = puck@cell_labels[barcodes]
  n_levels = length(levels(droplevels(puck@cell_labels[barcodes])))
  size_vec = c(3,rep(1, n_levels-1))
  gg <- ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) + ggplot2::geom_point(ggplot2::aes(shape=19,color=class,size=class)) +
    ggplot2::scale_color_manual(values = pals::kelly(n_levels+1)[2:(n_levels+1)])+ ggplot2::scale_shape_identity() + ggplot2::scale_size_manual(values=size_vec)+ggplot2::theme_bw()
  gg
}

plot_cell_types <- function(puck, barcodes, results_dir) {
  my_table = puck@coords[barcodes,]
  my_table$class = puck@cell_labels[barcodes]
  n_levels = length(levels(droplevels(puck@cell_labels[barcodes])))
  my_pal = pals::kelly(n_levels+1)[2:(n_levels+1)]
  if(n_levels > 21)
    my_pal = pals::polychrome(n_levels)
  if(n_levels > 36)
    stop("Plotting currently supports at most 36 cell types as colors")
  plot <- ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) + ggplot2::geom_point(ggplot2::aes(shape=19,color=class)) +
    ggplot2::scale_color_manual(values = my_pal)+ ggplot2::scale_shape_identity() + ggplot2::theme_bw()
  pdf(file.path(results_dir,"all_cell_types.pdf"))
  invisible(print(plot))
  dev.off()
}

#individually plots cell types
#if counter_stain = cell_type, then it also plots one cell type as a reference
plot_cell_types_ind <- function(puck, results_dir, counter_stain = NULL) {
  cell_types = levels(droplevels(puck@cell_labels))
  n_levels = length(cell_types)
  plots <- vector(mode = "list", length = n_levels)
  for(i in 1:n_levels) {
    cell_type = cell_types[i]
    curr_loc = puck@cell_labels==cell_type
    if(!is.null(counter_stain) && cell_type != counter_stain)
      curr_loc = curr_loc | puck@cell_labels==counter_stain
    my_table = puck@coords[curr_loc,]
    my_table$class = puck@cell_labels[curr_loc]
    plots[[i]] <- ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) + ggplot2::geom_point(ggplot2::aes(shape=19, color = class)) + ggplot2::scale_shape_identity() + ggplot2::theme_bw() + ggplot2::ggtitle(cell_type)
  }
  #l = mget(plots)
  if(!is.null(counter_stain))
    pdf(file.path(results_dir,paste0("cell_type_calls_counter_",counter_stain,".pdf")))
  else
    pdf(file.path(results_dir,"cell_type_calls.pdf"))
  invisible(lapply(plots, print))
  dev.off()
}


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
get_de_genes <- function(cell_type_means, puck, fc_thresh = 1.25, expr_thresh = .00015) {
  total_gene_list = c()
  epsilon = 1e-9
  bulk_vec = rowSums(puck@counts)
  gene_list = rownames(cell_type_means)
  gene_list = intersect(gene_list,names(bulk_vec))
  gene_list = gene_list[bulk_vec[gene_list] > 0]
  for(cell_type in cell_type_names) {
    other_mean = rowMeans(cell_type_means[gene_list,cell_type_names != cell_type])
    logFC = log(cell_type_means[gene_list,cell_type] + epsilon) - log(other_mean + epsilon)
    type_gene_list = which((logFC > fc_thresh) & (cell_type_means[gene_list,cell_type] > expr_thresh))
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
compute_proportions <- function(cell_labels, cell_type_names, nUMI = NULL, smoothing_par = 1/300) {
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
