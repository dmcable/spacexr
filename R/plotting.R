#' Plots all doublets in space
#'
#' saves as 'all_doublets.pdf'
#'
#' @param puck an object of type \linkS4class{SpatialRNA}
#' @param doublets a dataframe of RCTD results restricted to doublets
#' @param resultsdir output directory
#' @param cell_type_names list of cell type names
#' @return returns \code{\link{ggplot2}} object
#' @export
plot_doublets <- function(puck, doublets, resultsdir, cell_type_names) {
  barcodes = rownames(doublets)
  my_table = puck@coords[barcodes,]
  my_table$class = doublets$first_type
  my_table2 = puck@coords[barcodes,]
  jitter = 15
  my_table2$x = my_table2$x + jitter
  my_table2$class = doublets$second_type
  my_table = rbind(my_table, my_table2)
  n_levels = length(cell_type_names)
  my_pal = pals::kelly(n_levels+1)[2:(n_levels+1)]
  pres = unique(as.integer(my_table$class))
  pres = pres[order(pres)]
  if(n_levels > 21)
    my_pal = pals::polychrome(n_levels)
  if(n_levels > 36)
    stop("Plotting currently supports at most 36 cell types as colors")
  plot <- ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) + ggplot2::geom_point(ggplot2::aes(size = .15, shape=19,color=class)) +
    ggplot2::scale_color_manual(values = my_pal[pres])+ ggplot2::scale_shape_identity() + ggplot2::theme_bw() + ggplot2::scale_size_identity()
  pdf(file.path(resultsdir,"all_doublets.pdf"))
  invisible(print(plot))
  dev.off()
  return(plot)
}


#' Plots all cell types in space
#'
#' Plots the first cell type in doublet mode. Saves as 'all_cell_types.pdf'
#'
#' @param coords a dataframe of coordinates of each pixel
#' @param results_df a dataframe of RCTD results (see \code{\link{gather_results}})
#' @param resultsdir output directory
#' @param cell_type_names list of cell type names
#' @return returns \code{\link{ggplot2}} object
#' @export
plot_all_cell_types <- function(results_df, coords, cell_type_names, resultsdir) {
  barcodes = rownames(results_df[results_df$spot_class != "reject" & results_df$first_type %in% cell_type_names,])
  my_table = coords[barcodes,]
  my_table$class = results_df[barcodes,]$first_type
  n_levels = length(levels(my_table$class))
  my_pal = pals::kelly(n_levels+1)[2:(n_levels+1)]
  pres = unique(as.integer(my_table$class))
  pres = pres[order(pres)]
  if(n_levels > 21)
    my_pal = pals::polychrome(n_levels)
  if(n_levels > 36)
    stop("Plotting currently supports at most 36 cell types as colors")
  plot <- ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) + ggplot2::geom_point(ggplot2::aes(size = .15, shape=19,color=class)) +
    ggplot2::scale_color_manual(values = my_pal[pres])+ ggplot2::scale_shape_identity() + ggplot2::theme_bw() + ggplot2::scale_size_identity()
  pdf(file.path(resultsdir,"all_cell_types.pdf"))
  invisible(print(plot))
  dev.off()
  return(plot)
}


#' Plots doublets of each cell type individually
#'
#' Plots the first cell type in doublet mode. Saves as 'all_doublets_type.pdf'
#'
#' @param puck an object of type \linkS4class{SpatialRNA}
#' @param doublets_base a dataframe of RCTD results restricted to doublets
#' @param resultsdir output directory
#' @param cell_type_names list of cell type names
#' @export
plot_doublets_type <- function(puck, doublets_base, resultsdir, cell_type_names) {
  plots <- vector(mode = "list", length = length(cell_type_names))
  i = 1
  for (cell_type in cell_type_names) {
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
      n_levels = length(cell_type_names)
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
  pdf(file.path(resultsdir,"all_doublets_type.pdf"))
  invisible(lapply(plots, print))
  dev.off()
}


#' Plots a continuous value over locations on the puck
#'
#' Colors points based on value of the function
#'
#' @param puck an object of type \linkS4class{SpatialRNA}
#' @param barcodes a list of barcodes to include in the plot
#' @param plot_val a named (by barcode) list of values to plot
#' @param ylimit minimum and maximum values for the range of plot as a numeric list
#' @param ylim (optional) minimum and maximum value for y coordinate as a numeric list
#' @param xlim (optional) minimum and maximum value for x coordinate as a numeric list
#' @param size numeric size of points
#' @param cell_type_info cell type information and profiles (see \code{\link{get_cell_type_info}})
#'
#' @return Returns a \link{ggplot} object
#' @export
plot_puck_continuous <- function(puck, barcodes, plot_val, ylimit = c(0,1),title = NULL, counter_barcodes = NULL, label = F, my_pal = NULL, xlim = NULL, ylim = NULL, size=0.15, alpha = 1, small_point = F) {
  if(is.null(my_pal))
    my_pal = pals::kovesi.rainbow(20)
  my_table = puck@coords[barcodes,]
  plot_val <- pmax(pmin(plot_val, ylimit[2] -1e-8),ylimit[1] +1e-8)
  my_table$value = plot_val[barcodes]
  if(!is.null(ylimit))
    sc <- ggplot2::scale_colour_gradientn(colors = my_pal, limits = ylimit)
  else
    sc <- ggplot2::scale_colour_gradientn(colors = my_pal)
  plot <- ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y))
  if(small_point)
    plot <- plot + ggplot2::geom_point(ggplot2::aes(size = size, shape=16,color=value, stroke = 0),alpha = alpha)
  else
    plot <- plot + ggplot2::geom_point(ggplot2::aes(size = size, shape=19,color=value),alpha = alpha)
  plot <- plot + sc + ggplot2::scale_shape_identity() + ggplot2::theme_classic() + ggplot2::scale_size_identity()
  if(label)
    plot <- plot + ggplot2::aes(label = which(rownames(puck@coords) %in% inter_barcodes[my_class==i])) + ggplot2::geom_text()
  if(!is.null(counter_barcodes)) {
    my_table = puck@coords[counter_barcodes,]
    plot <- plot + ggplot2::geom_point(data=my_table,ggplot2::aes(x=x, y=y,size=0.15),alpha = 0.1)
  }
  if(is.null(xlim))
    xlim <- c(min(puck@coords$x) - 1,max(puck@coords$x) + 1)
  if(is.null(ylim))
    ylim <- c(min(puck@coords$y) - 1,max(puck@coords$y) + 1)
  plot <- plot + ggplot2::coord_fixed() + ggplot2::xlim(xlim) + ggplot2::ylim(ylim)
  if(!is.null(title))
    plot <- plot + ggplot2::ggtitle(title)
  plot
}

#' Plots a continuous value over filtered locations on the puck
#'
#' Colors points based on value of the function, filtered for e.g. UMI and cell type
#'
#' @param puck an object of type \linkS4class{SpatialRNA}
#' @param cell_type string specifying cell type to plot. if NULL, then all cell types are plotted
#' @param plot_val a named (by barcode) list of values to plot
#' @param min_val numeric, minimum value for the range of plot as a numeric list
#' @param max_val numeric, maximum value for the range of plot as a numeric list
#' @param minUMI numeric, minimum value for total UMIs to filter pixels
#' @param maxUMI numeric, maximum value for total UMIs to filter pixels
#' @param cell_type_info cell type information and profiles (see \code{\link{get_cell_type_info}})
#'
#' @return Returns a \link{ggplot} object
#' @export
plot_puck_wrapper <- function(puck, plot_val, cell_type = NULL, minUMI = 0, maxUMI = 200000, min_val = NULL, max_val = NULL, title = NULL, my_cond = NULL) {
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
    epsilon = 0.00001
    plot_val[plot_val >= max_val - epsilon] = max_val - epsilon
    if(!is.null(min_val))
      ylimit = c(min_val, max_val)
  }
  plot_puck_continuous(puck, names(which(my_cond)), plot_val, title = title, ylimit = ylimit)
}


#' Spatially plot the confident weights for each cell type
#'
#' Plots the confident weights for each cell type as in full_mode. Saves as 'cell_type_weights.pdf'
#'
#' @param puck an object of type \linkS4class{SpatialRNA}
#' @param weights a dataframe of RCTD output weights (see \code{\link{gather_results}})
#' @param resultsdir output directory
#' @param cell_type_names list of cell type names
#' @export
plot_weights <- function(cell_type_names, puck, resultsdir, weights) {
  plots <- vector(mode = "list", length = length(cell_type_names))
  for (i in 1:length(cell_type_names)) {
    cell_type = cell_type_names[i]
    my_cond = weights[,cell_type] > UMI_cutoff(puck@nUMI) # pmax(0.25, 0.7 - puck@nUMI / 1000)
    plot_var <- weights[,cell_type]; names(plot_var) = rownames(weights)
    if(sum(my_cond) > 0)
      plots[[i]] <- plot_puck_wrapper(puck, plot_var, NULL, minUMI = 100,maxUMI = 200000,min_val = 0, max_val = 1, title = cell_type, my_cond = my_cond)
  }
  pdf(file.path(resultsdir,"cell_type_weights.pdf"))
  invisible(lapply(plots, print))
  dev.off()
}

#' Spatially plot all weights for each cell type
#'
#' Plots all weights for each cell type as in full_mode. Saves as 'cell_type_weights_unthreshold.pdf'
#'
#' @param puck an object of type \linkS4class{SpatialRNA}
#' @param weights a dataframe of RCTD output weights (see \code{\link{gather_results}})
#' @param resultsdir output directory
#' @param cell_type_names list of cell type names
#' @export
plot_weights_unthreshold <- function(cell_type_names, puck, resultsdir, weights) {
  plots <- vector(mode = "list", length = length(cell_type_names))
  for (i in 1:length(cell_type_names)) {
    cell_type = cell_type_names[i]
    plot_var <- weights[,cell_type]; names(plot_var) = rownames(weights)
    if(sum(weights[,cell_type]) > 0)
      plots[[i]] <- plot_puck_wrapper(puck, plot_var, NULL, minUMI = 100,maxUMI = 200000,min_val = 0, max_val = 1, title = cell_type)
  }
  pdf(file.path(resultsdir,"cell_type_weights_unthreshold.pdf"))
  invisible(lapply(plots, print))
  dev.off()
}


#' Barplot of the confident counts for each cell type
#'
#' Plots the number of confident labels in 'full_mode'. Saves as 'cell_type_occur.pdf'
#'
#' @param weights a dataframe of RCTD output weights (see \code{\link{gather_results}})
#' @param resultsdir output directory
#' @param puck an object of type \linkS4class{SpatialRNA}
#' @param cell_type_names list of cell type names
#' @return returns \code{\link{ggplot2}} object
#' @export
plot_cond_occur <- function(cell_type_names, resultsdir, weights, puck) {
  occur <- numeric(length(cell_type_names))
  names(occur) = cell_type_names
  for (i in 1:length(cell_type_names)) {
    cell_type = cell_type_names[i]
    my_cond = weights[,cell_type] > UMI_cutoff(puck@nUMI)
    occur[cell_type] = sum(my_cond)
  }
  df<- reshape2::melt(as.list(occur)); colnames(df) = c('Count','Cell_Type')
  pdf(file.path(resultsdir,"cell_type_occur.pdf"))
  plot<-ggplot2::ggplot(df, ggplot2::aes(x=Cell_Type, y=Count, fill=Cell_Type)) +
    ggplot2::geom_bar(stat="identity")+ggplot2::theme_minimal()
  invisible(print(plot))
  dev.off()
  return(plot)
}

#' Barplot of the counts for each cell type
#'
#' Plots the number of (including unconfident) labels in 'full_mode'. Saves as 'cell_type_occur_unthreshold.pdf'
#'
#' @param weights a dataframe of RCTD output weights (see \code{\link{gather_results}})
#' @param resultsdir output directory
#' @param cell_type_info cell type information and profiles (see \code{\link{get_cell_type_info}})
#' @export
plot_occur_unthreshold <- function(cell_type_info, resultsdir, weights) {
  occur <- table(apply(weights, 1, which.max))[as.character(1:cell_type_info[[3]])]
  names(occur) = cell_type_info[[2]]
  df<- reshape2::melt(as.list(occur)); colnames(df) = c('Count','Cell_Type')
  pdf(file.path(resultsdir,"cell_type_occur_unthreshold.pdf"))
  plot<-ggplot2::ggplot(df, ggplot2::aes(x=Cell_Type, y=Count, fill=Cell_Type)) +
    ggplot2::geom_bar(stat="identity")+ggplot2::theme_minimal()
  invisible(print(plot))
  dev.off()
}

#' Spatially plot the weights for each cell type in doublet_mode
#'
#' Plots the weights for each cell type as in doublet_mode. Saves as 'cell_type_weights_doublet.pdf'
#'
#' @param puck an object of type \linkS4class{SpatialRNA}
#' @param weights_doublet a dataframe of RCTD output weights for doublets (see \code{\link{gather_results}})
#' @param resultsdir output directory
#' @param cell_type_names list of cell type names
#' @param results_df dataframe of RCTD results (see \code{\link{gather_results}})
#' @export
plot_weights_doublet <- function(cell_type_names, puck, resultsdir, weights_doublet, results_df) {
  plots <- vector(mode = "list", length = length(cell_type_names))
  for (i in 1:length(cell_type_names)) {
    cell_type = cell_type_names[i]
    all_weights <- weights_doublet[results_df$spot_class == "doublet_certain" & results_df$second_type == cell_type,2, drop=FALSE]
    all_weights <- rbind(all_weights, weights_doublet[!(results_df$spot_class == "reject") & results_df$first_type == cell_type,1,drop=FALSE])
    all_weights_vec <- as.vector(all_weights); names(all_weights_vec) <- rownames(all_weights)
    if(length(all_weights) > 0)
      plots[[i]] <- plot_puck_continuous(puck, rownames(all_weights), all_weights_vec, title = cell_type, ylimit = c(0,1))
    else
      plots[[i]] <- NULL
  }
  pdf(file.path(resultsdir,"cell_type_weights_doublets.pdf"))

  invisible(lapply(plots, print))
  dev.off()
}

#' Plots doublet co-occurances
#'
#' Plots the doublet co-occurances. Saves as 'doublet_stacked_bar.pdf'
#'
#' @param doub_occur a table of occurances of doublets
#' @param resultsdir output directory
#' @param cell_type_names list of cell type names
#' @return returns \code{\link{ggplot2}} object
#' @export
plot_doub_occur_stack <- function(doub_occur, resultsdir, cell_type_names) {
  cell_type_names <- unlist(lapply(cell_type_names, function(x) paste0('celltype_',x)))
  rownames(doub_occur) <- cell_type_names; colnames(doub_occur) <- cell_type_names;
  data <- reshape2::melt(doub_occur)
  colnames(data) = c('second_type','first_type','count')
  n_levels = length(cell_type_names)
  my_pal = pals::kelly(n_levels+1)[2:(n_levels+1)]
  names(my_pal) = cell_type_names
  pres = cell_type_names
  pres = pres[order(pres)]
  pdf(file.path(resultsdir,'doublet_stacked_bar.pdf'))
  plot <- ggplot2::ggplot(data, ggplot2::aes(fill=second_type, y=count, x=first_type)) +ggplot2::geom_bar(position="stack", stat="identity") + ggplot2::scale_fill_manual(values = my_pal[pres]) + ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1, angle = 45))
  invisible(print(plot))
  dev.off()
  return(plot)
}

#' Plots a factor variable in space on the puck
#'
#' Colors points based on class
#'
#' @param puck an object of type \linkS4class{SpatialRNA}
#' @param barcodes_cur a list of barcodes to include in the plot
#' @param my_class a named (by barcode) factor list for the coloring
#'
#' @return Returns a \link{ggplot} object
#' @export
plot_class <- function(puck, barcodes_cur, my_class, counter_barcodes = NULL, title = NULL) {
  my_table = puck@coords[barcodes_cur,]
  my_table$class = my_class[barcodes_cur]
  n_levels <- length(unique(my_class))
  if(n_levels > 36)
    cols = rainbow(n_levels)[sample(1:n_levels, n_levels, replace = F)]
  else {
    if(n_levels > 21)
      cols = unname(pals::polychrome(n_levels+2))[3:(n_levels+2)]
    else
      cols = pals::kelly(n_levels+2)[3:(n_levels+2)]
  }
  plot <- ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) + ggplot2::geom_point(ggplot2::aes(size = 0.5, shape=19,color=class)) +
    ggplot2::scale_shape_identity() + ggplot2::theme_classic() + ggplot2::scale_size_identity() + ggplot2::scale_color_manual(values = cols)
  if(!is.null(counter_barcodes)) {
    my_table = puck@coords[counter_barcodes,]
    plot <- plot + ggplot2::geom_point(data=my_table,ggplot2::aes(x=x, y=y,size=0.15),alpha = 0.1)
  }
  if(!is.null(title))
    plot <- plot + ggplot2::ggtitle(title)
  plot
}

#' Create all plots for an RCTD object after cell types have been assigned
#'
#' @param myRCTD a \code{\linkS4class{RCTD}} object with annotated cell types e.g. from the \code{\link{run.RCTD}} function.
#' @param datadir directory where plots should be saved
#' @export
create_RCTD_plots <- function(myRCTD, datadir) {
  results <- myRCTD@results
  norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/')
  cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
  spatialRNA <- myRCTD@spatialRNA
  resultsdir <- file.path(datadir,'RCTD_Plots') ## you may change this to a more accessible directory on your computer.
  if(!dir.exists(resultsdir))
    dir.create(resultsdir)
  plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
  # Plots all weights for each cell type as in full_mode. (saved as
  # 'results/cell_type_weights.pdf')
  plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights)
  # Plots the weights for each cell type as in doublet_mode. (saved as
  # 'results/cell_type_weights_doublets.pdf')
  plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet,
                       results$results_df)
  # Plots the number of confident pixels of each cell type in 'full_mode'. (saved as
  # 'results/cell_type_occur.pdf')
  plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
  # makes a map of all cell types, (saved as
  # 'results/all_cell_types.pdf')
  plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir)
}
