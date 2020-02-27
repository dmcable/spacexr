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

plot_cell_types <- function(puck, barcodes, results_dir) {
  my_table = puck@coords[barcodes,]
  my_table$class = puck@cell_labels[barcodes]
  n_levels = length(levels(droplevels(puck@cell_labels[barcodes])))
  my_pal = pals::kelly(n_levels+1)[2:(n_levels+1)]
  if(n_levels > 21)
    my_pal = pals::polychrome(n_levels)
  if(n_levels > 36)
    stop("Plotting currently supports at most 36 cell types as colors")
  plot <- ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) + ggplot2::geom_point(ggplot2::aes(size = .15, shape=19,color=class)) +
    ggplot2::scale_color_manual(values = my_pal)+ ggplot2::scale_shape_identity() + ggplot2::theme_bw() + ggplot2::scale_size_identity()
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

#plots a continuous value over the puck
plot_puck_continuous <- function(puck, barcodes, plot_val, title = NULL) {
  my_pal = pals::kovesi.rainbow(20)
  my_table = puck@coords[barcodes,]
  my_table$value = plot_val[barcodes]
  plot <- ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) + ggplot2::geom_point(ggplot2::aes(size = 0.5, shape=19,color=value)) +
    ggplot2::scale_colour_gradientn(colors = my_pal) + ggplot2::scale_shape_identity() + ggplot2::theme_bw() + ggplot2::scale_size_identity()
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
  if(!is.null(my_cond))
    my_cond = UMI_filter & my_cond
  else
    my_cond = UMI_filter
  if(!is.null(cell_type))
    my_cond = my_cond & (puck@cell_labels == cell_type)
  if(!is.null(min_val))
    my_cond = my_cond & (plot_val > min_val)
  if(!is.null(max_val))
    plot_val[plot_val > max_val] = max_val
  plot_puck_continuous(puck, names(which(my_cond)), plot_val, title = title)
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
    plots[[i]] <- plot_puck_wrapper(puck, weights[,cell_type], NULL, minUMI = 100,min_val = 0, max_val = 1, title = cell_type, my_cond = my_cond)
  }
  pdf(file.path(resultsdir,"cell_type_weights.pdf"))
  invisible(lapply(plots, print))
  dev.off()
}

#spatially plot the weights for each cell type unthresholded
plot_weights_unthreshold <- function(cell_type_info, puck, resultsdir, weights) {
  plots <- vector(mode = "list", length = cell_type_info[[3]])
  for (i in 1:cell_type_info[[3]]) {
    cell_type = cell_type_info[[2]][i]
    plots[[i]] <- plot_puck_wrapper(puck, weights[,cell_type], NULL, minUMI = 100,min_val = 0, max_val = 1, title = cell_type)
  }
  pdf(file.path(resultsdir,"cell_type_weights_unthreshold.pdf"))
  invisible(lapply(plots, print))
  dev.off()
}


#plots the distribution of weights for each cell type across UMI
plot_weight_distribution <- function(cell_type_info, puck, resultsdir) {
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
  plot(x_vals, labeled, type = 'n', xlab = "log2(N_UMI)", ylab = "Percent called as singlets")
  lines(x_vals, labeled)
  lines(x_vals, (1:length(x_vals)) / length(x_vals), col ="red")
  dev.off()
}

#plot meta gene scores
plot_meta_genes <- function(cell_type_info, marker_scores_df, resultsdir) {
  plots <- vector(mode = "list", length = cell_type_info[[3]])
  for(i in 1:n_cell_types) {
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
  for(i in 1:n_cell_types) {
    cell_type = cell_type_info[[2]][[i]]
    plots[[i]] <- plot_puck_wrapper(puck, norm_marker_scores[,cell_type], NULL, minUMI = 100, min_val = 0, max_val = 1, title = cell_type)
  }
  pdf(file.path(resultsdir,"normalized_meta_gene_scores.pdf"))
  invisible(lapply(plots, print))
  dev.off()
}

