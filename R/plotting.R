#simple moving average
ma <- function(x, n = 5){stats::filter(x, rep(1 / n, n), sides = 2)}

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

plot_doublets <- function(doublets, results_dir, cell_type_info) {
  barcodes = rownames(doublets)
  my_table = puck@coords[barcodes,]
  my_table$class = doublets$first_type
  my_table2 = puck@coords[barcodes,]
  jitter = 15
  my_table2$x = my_table2$x + jitter
  my_table2$class = doublets$second_type
  my_table = rbind(my_table, my_table2)
  n_levels = cell_type_info[[3]]
  my_pal = pals::kelly(n_levels+1)[2:(n_levels+1)]
  pres = unique(as.integer(my_table$class))
  pres = pres[order(pres)]
  if(n_levels > 21)
    my_pal = pals::polychrome(n_levels)
  if(n_levels > 36)
    stop("Plotting currently supports at most 36 cell types as colors")
  plot <- ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) + ggplot2::geom_point(ggplot2::aes(size = .15, shape=19,color=class)) +
    ggplot2::scale_color_manual(values = my_pal[pres])+ ggplot2::scale_shape_identity() + ggplot2::theme_bw() + ggplot2::scale_size_identity()
  pdf(file.path(results_dir,"all_doublets.pdf"))
  invisible(print(plot))
  dev.off()
}

plot_doublets_type <- function(doublets_base, results_dir, cell_type_info) {
  plots <- vector(mode = "list", length = cell_type_info[[3]])
  i = 1
  for (cell_type in cell_type_info[[2]]) {
    doublets = doublets_base[doublets_base$first_type == cell_type | doublets_base$second_type == cell_type,]
    barcodes = rownames(doublets)
    if(length(barcodes) > 0) {
      my_table = puck@coords[barcodes,]
      my_table$class = doublets$first_type
      my_table2 = puck@coords[barcodes,]
      jitter = 15
      my_table2$x = my_table2$x + jitter
      my_table2$class = doublets$second_type
      my_table = rbind(my_table, my_table2)
      n_levels = cell_type_info[[3]]
      my_pal = pals::kelly(n_levels+1)[2:(n_levels+1)]
      pres = unique(as.integer(my_table$class))
      pres = pres[order(pres)]
      if(n_levels > 21)
        my_pal = pals::polychrome(n_levels)
      if(n_levels > 36)
        stop("Plotting currently supports at most 36 cell types as colors")
      plots[[i]] <- ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) + ggplot2::geom_point(ggplot2::aes(size = .15, shape=19,color=class)) + ggplot2::scale_color_manual(values = my_pal[pres])+ ggplot2::scale_shape_identity() + ggplot2::theme_bw() + ggplot2::scale_size_identity() +ggplot2::ggtitle(cell_type)
    }
    i = i + 1
  }
  pdf(file.path(results_dir,"all_doublets_type.pdf"))
  invisible(lapply(plots, print))
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

#plots a continuous value over the puck
plot_puck_continuous <- function(puck, barcodes, plot_val, title = NULL, ylimit = NULL) {
  my_pal = pals::kovesi.rainbow(20)
  my_table = puck@coords[barcodes,]
  my_table$value = plot_val[barcodes]
  if(!is.null(ylimit))
    sc <- ggplot2::scale_colour_gradientn(colors = my_pal, limits = ylimit)
  else
    sc <- ggplot2::scale_colour_gradientn(colors = my_pal)
  plot <- ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) + ggplot2::geom_point(ggplot2::aes(size = 0.5, shape=19,color=value)) +
     sc + ggplot2::scale_shape_identity() + ggplot2::theme_bw() + ggplot2::scale_size_identity()
  if(!is.null(title))
    plot <- plot + ggplot2::ggtitle(title)
  plot
  #pdf(file.path(results_dir,"all_cell_types.pdf"))
  #invisible(print(plot))
  #dev.off()
}

#plots a continunous value over certain beads in the puck
plot_puck_wrapper <- function(puck, plot_val, cell_type = NULL, minUMI = 0, maxUMI = 20000, min_val = NULL, max_val = NULL, title = NULL, my_cond = NULL) {
  UMI_filter = (puck@nUMI > minUMI) & (puck@nUMI < maxUMI)
  ylimit = NULL
  if(!is.null(my_cond))
    my_cond = UMI_filter & my_cond
  else
    my_cond = UMI_filter
  if(!is.null(cell_type))
    my_cond = my_cond & (puck@cell_labels == cell_type)
  if(!is.null(min_val))
    my_cond = my_cond & (plot_val > min_val)
  if(!is.null(max_val)) {
    plot_val[plot_val > max_val] = max_val
    if(!is.null(min_val))
      ylimit = c(min_val, max_val)
  }
  plot_puck_continuous(puck, names(which(my_cond)), plot_val, title = title, ylimit = ylimit)
}

#plots a gene over the puck. Positive restricts plotting to positive.
plot_puck_gene <- function(puck, gene, cell_type = NULL, minUMI = 0, positive = F) {
  gene_vals = puck@counts[gene,]
  plot_puck_wrapper(puck, gene_vals, cell_type, minUMI, min_val = positive - 0.5)
}

#spatially plot the weights for each cell type
plot_weights <- function(cell_type_info, puck, resultsdir, weights) {
  plots <- vector(mode = "list", length = cell_type_info[[3]])
  for (i in 1:cell_type_info[[3]]) {
    cell_type = cell_type_info[[2]][i]
    my_cond = weights[,cell_type] > UMI_cutoff(puck@nUMI) # pmax(0.25, 0.7 - puck@nUMI / 1000)
    plot_var <- weights[,cell_type]; names(plot_var) = rownames(weights)
    if(sum(my_cond) > 0)
      plots[[i]] <- plot_puck_wrapper(puck, plot_var, NULL, minUMI = 100,min_val = 0, max_val = 1, title = cell_type, my_cond = my_cond)
  }
  pdf(file.path(resultsdir,"cell_type_weights.pdf"))
  invisible(lapply(plots, print))
  dev.off()
}

#spatially plot the weights for each cell type
plot_weights_nmf <- function(cell_type_info, puck, resultsdir, weights, thresh_nmf) {
  plots <- vector(mode = "list", length = cell_type_info[[3]])
  for (i in 1:cell_type_info[[3]]) {
    cell_type = cell_type_info[[2]][i]
    my_cond = weights[,cell_type] > thresh_nmf[cell_type,] # pmax(0.25, 0.7 - puck@nUMI / 1000)
    plot_var <- weights[,cell_type]; names(plot_var) = rownames(weights)
    if(sum(my_cond) > 0)
      plots[[i]] <- plot_puck_wrapper(puck, plot_var, NULL, minUMI = 100,min_val = 0, max_val = 1, title = cell_type, my_cond = my_cond)
  }
  pdf(file.path(resultsdir,"cell_type_weights.pdf"))
  invisible(lapply(plots, print))
  dev.off()
}

#plot the confident counts for each cell type
plot_cond_occur <- function(cell_type_info, resultsdir, weights) {
  occur <- numeric(cell_type_info[[3]])
  names(occur) = cell_type_info[[2]]
  for (i in 1:cell_type_info[[3]]) {
    cell_type = cell_type_info[[2]][i]
    my_cond = weights[,cell_type] > UMI_cutoff(puck@nUMI)
    occur[cell_type] = sum(my_cond)
  }
  df<- melt(as.list(occur)); colnames(df) = c('Count','Cell_Type')
  pdf(file.path(resultsdir,"cell_type_occur.pdf"))
  plot<-ggplot(df, aes(x=Cell_Type, y=Count, fill=Cell_Type)) +
    geom_bar(stat="identity")+theme_minimal()
  invisible(print(plot))
  dev.off()
}

#plot the confident counts for each cell type
plot_occur_unthreshold <- function(cell_type_info, resultsdir, weights) {
  occur <- table(apply(weights, 1, which.max))[as.character(1:cell_type_info[[3]])]
  names(occur) = cell_type_info[[2]]
  df<- melt(as.list(occur)); colnames(df) = c('Count','Cell_Type')
  pdf(file.path(resultsdir,"cell_type_occur_unthreshold.pdf"))
  plot<-ggplot(df, aes(x=Cell_Type, y=Count, fill=Cell_Type)) +
    geom_bar(stat="identity")+theme_minimal()
  invisible(print(plot))
  dev.off()
}

#plot the confident counts for each cell type
plot_cond_occur_nmf <- function(cell_type_info, resultsdir, weights, thresh_nmf) {
  occur <- numeric(cell_type_info[[3]])
  names(occur) = cell_type_info[[2]]
  for (i in 1:cell_type_info[[3]]) {
    cell_type = cell_type_info[[2]][i]
    my_cond = weights[,cell_type] > thresh_nmf[cell_type,]
    occur[cell_type] = sum(my_cond)
  }
  df<- melt(as.list(occur)); colnames(df) = c('Count','Cell_Type')
  pdf(file.path(resultsdir,"cell_type_occur.pdf"))
  plot<-ggplot(df, aes(x=Cell_Type, y=Count, fill=Cell_Type)) +
    geom_bar(stat="identity")+theme_minimal()
  invisible(print(plot))
  dev.off()
}

#spatially plot the weights_doublet for each cell type
plot_weights_doublet <- function(cell_type_info, puck, resultsdir, weights_doublet) {
  plots <- vector(mode = "list", length = cell_type_info[[3]])
  for (i in 1:cell_type_info[[3]]) {
    cell_type = cell_type_info[[2]][i]
    all_weights <- weights_doublet[results_df$spot_class == "doublet_certain" & results_df$second_type == cell_type,2]
    all_weights <- c(all_weights, weights_doublet[!(results_df$spot_class == "reject") & results_df$first_type == cell_type,1])
    plots[[i]] <- plot_puck_continuous(puck, names(all_weights), all_weights, title = cell_type)
  }
  pdf(file.path(resultsdir,"cell_type_weights_doublets.pdf"))
  invisible(lapply(plots, print))
  dev.off()
}

#spatially plot the weights for each cell type unthresholded
plot_weights_unthreshold <- function(cell_type_info, puck, resultsdir, weights) {
  plots <- vector(mode = "list", length = cell_type_info[[3]])
  for (i in 1:cell_type_info[[3]]) {
    cell_type = cell_type_info[[2]][i]
    plot_var <- weights[,cell_type]; names(plot_var) = rownames(weights)
    if(sum(weights[,cell_type]) > 0)
      plots[[i]] <- plot_puck_wrapper(puck, plot_var, NULL, minUMI = 100,min_val = 0, max_val = 1, title = cell_type)
  }
  pdf(file.path(resultsdir,"cell_type_weights_unthreshold.pdf"))
  invisible(lapply(plots, print))
  dev.off()
}


#plots the distribution of weights for each cell type across UMI
plot_weight_distribution <- function(cell_type_info, puck, resultsdir, weights) {
  pdf(file.path(resultsdir,"cell_type_weight_distribution.pdf"))
  for (i in 1:cell_type_info[[3]]) {
    cell_type = cell_type_info[[2]][i]
    x_vals = log(puck@nUMI,2)
    y_vals = weights[,cell_type]
    z_vals = UMI_cutoff(puck@nUMI)
    z_vals = z_vals[order(x_vals)]
    y_vals = y_vals[order(x_vals)]
    x_vals = x_vals[order(x_vals)]
    plot(x_vals, y_vals, xlab = "log2(N_UMI)", ylab = "weight", main = cell_type, pch = 20)
    lines(x_vals, z_vals, col = "red")
  }
  dev.off()
}

#confidence rate
plot_confidence_rate <- function(puck, resultsdir, weights) {
  pdf(file.path(resultsdir,"confidence_rate.pdf"))
  x_vals = log(puck@nUMI,2)
  y_vals = apply(weights, 1, max) > UMI_cutoff(puck@nUMI)
  y_vals = y_vals[order(x_vals)]
  labeled = ma(y_vals,500)
  x_vals = x_vals[order(x_vals)]
  plot(x_vals, labeled, type = 'n', xlab = "log2(N_UMI)", ylab = "Percent called as singlets", ylim = c(0, 1))
  lines(x_vals, labeled)
  lines(x_vals, (1:length(x_vals)) / length(x_vals), col ="red")
  dev.off()
}

#plot meta gene scores
plot_meta_genes <- function(cell_type_info, marker_scores_df, resultsdir) {
  plots <- vector(mode = "list", length = cell_type_info[[3]])
  for(i in 1:cell_type_info[[3]]) {
    cell_type = cell_type_info[[2]][[i]]
    plots[[i]] <- plot_puck_wrapper(puck, marker_scores_df[,cell_type], NULL, minUMI = 100, min_val = 0, max_val = 3, title = cell_type)
  }
  pdf(file.path(resultsdir,"meta_gene_scores.pdf"))
  invisible(lapply(plots, print))
  dev.off()
}

#plot norm marker gene scores
plot_norm_meta_genes <- function(cell_type_info, norm_marker_scores, resultsdir) {
  plots <- vector(mode = "list", length = cell_type_info[[3]])
  for(i in 1:cell_type_info[[3]]) {
    cell_type = cell_type_info[[2]][[i]]
    plots[[i]] <- plot_puck_wrapper(puck, norm_marker_scores[,cell_type], NULL, minUMI = 100, min_val = 0, max_val = 1, title = cell_type)
  }
  pdf(file.path(resultsdir,"normalized_meta_gene_scores.pdf"))
  invisible(lapply(plots, print))
  dev.off()
}

