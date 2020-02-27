library(gplots)
library("RColorBrewer")

plot_heat_map <- function(data, file_loc = NULL, save.file = F, normalize = T) {
  if(normalize) {
    norm_conf = sweep(data, 2, colSums(data), '/')
    norm_conf[norm_conf < 0] = 0
  } else
    norm_conf = data
  if(save.file) {
    png(file = file_loc)
    heatmap.2(norm_conf,col=rev(heat.colors(100)), breaks=seq(0,1,0.01),scale = "none",trace="none", Rowv=FALSE, Colv=FALSE,dendrogram='none')
    dev.off()
  } else {
    heatmap.2(norm_conf,col=rev(heat.colors(100)), breaks=seq(0,1,0.01), scale = "none",trace="none", Rowv=FALSE, Colv=FALSE,dendrogram='none')
  }
}

get_id_plot <- function(plot_df) {
  ggplot2::ggplot() +
    ggplot2::geom_line(data = plot_df, aes(x = nUMI, y = first_type), color = "blue") +
    ggplot2::geom_line(data = plot_df, aes(x = nUMI, y = second_type), color = "red") +
    ggplot2::xlab('nUMI') + ggplot2::ylab('percent present') + ggplot2::ylim(0,1)
}

plot_doublet_identification <- function(meta_data, common_cell_types, resultsdir, meta_df, results_df,
                                        class_df, UMI_list, use_class = T, toss_uncertain = T) {
  n_levels = choose(meta_data$n_cell_types,2)
  plots <- vector(mode = "list", length = n_levels)
  index = 1
  square_results = Matrix(0, nrow = meta_data$n_cell_types, ncol = meta_data$n_cell_types)
  rownames(square_results) = common_cell_types[1:meta_data$n_cell_types]
  colnames(square_results) = common_cell_types[1:meta_data$n_cell_types]
  plot_df_avg = NULL
  for (i in (1:(meta_data$n_cell_types-1))) {
    for (j in ((i + 1):meta_data$n_cell_types)) {
      type1 = common_cell_types[i]
      type2 = common_cell_types[j]
      curr_barcodes = meta_df$first_type == type1 & meta_df$second_type == type2
      if(toss_uncertain) {
        curr_barcodes = curr_barcodes & (results_df$spot_class != "reject")
      }
      if(use_class) {
        first_type_class = class_df[results_df[curr_barcodes,]$first_type,"class"]
        second_type_class = class_df[results_df[curr_barcodes,]$second_type,"class"]
        type1_class = class_df[type1,"class"]; type2_class = class_df[type2,"class"]
        first_type_pres <- first_type_class == type1_class | second_type_class == type1_class
        second_type_pres <- first_type_class == type2_class | second_type_class == type2_class
      } else {
        first_type_pres <- results_df[curr_barcodes,]$first_type == type1 | results_df[curr_barcodes,]$second_type == type1
        second_type_pres <- results_df[curr_barcodes,]$first_type == type2 | results_df[curr_barcodes,]$second_type == type2
      }
      plot_df <- aggregate(first_type_pres, list(meta_df[curr_barcodes, "first_UMI"]), mean)
      colnames(plot_df) = c("nUMI", "first_type")
      rownames(plot_df) = plot_df$nUMI
      plot_df$second_type = aggregate(second_type_pres, list(meta_df[curr_barcodes, "first_UMI"]), mean)$x
      if(is.null(plot_df_avg))
        plot_df_avg = plot_df
      else
        plot_df_avg = plot_df + plot_df_avg
      plots[[index]] = get_id_plot(plot_df) + labs(title = paste(type1,type2))
      doub_UMI_list <- as.character(UMI_list[2:(length(UMI_list)-1)])
      square_results[type1,type2] = colMeans(plot_df[doub_UMI_list,])["first_type"]
      square_results[type2, type1] = colMeans(plot_df[doub_UMI_list,])["second_type"]
      index = index + 1
    }
  }
  if(use_class) {
    pdf(file.path(resultsdir,"doublet_classification_class.pdf"))
  } else {
    pdf(file.path(resultsdir,"doublet_classification.pdf"))
  }
  invisible(lapply(plots, print))
  dev.off()
  if(use_class) {
    pdf(file.path(resultsdir,"avg_doublet_classification_class.pdf"))
  } else {
    pdf(file.path(resultsdir,"avg_doublet_classification.pdf"))
  }
  plot_df_avg = plot_df_avg / (index - 1)
  print(get_id_plot(plot_df_avg) + labs(title = paste("First Type", "Second Type")))
  dev.off()
  return(square_results)
}

#eg: equals_class(results_df, curr_barcodes, use_class, type1, "first_type")
equals_class <- function(results_df, curr_barcodes, use_class, my_type, feature) {
  if(use_class) {
    type_class = class_df[results_df[curr_barcodes,feature],"class"]
    my_type_class = class_df[my_type,"class"];
    return(type_class == my_type_class)
  } else {
    return(results_df[curr_barcodes,feature] == my_type)
  }
}

#only certain ones
plot_doublet_identification_certain <- function(meta_data, common_cell_types, resultsdir, meta_df, results_df,
                                        class_df, UMI_list, use_class = T) {
  n_levels = choose(meta_data$n_cell_types,2)
  plots <- vector(mode = "list", length = n_levels)
  index = 1
  square_results = Matrix(0, nrow = meta_data$n_cell_types, ncol = meta_data$n_cell_types)
  rownames(square_results) = common_cell_types[1:meta_data$n_cell_types]
  colnames(square_results) = common_cell_types[1:meta_data$n_cell_types]
  plot_df_avg = NULL
  for (i in (1:(meta_data$n_cell_types-1))) {
    for (j in ((i + 1):meta_data$n_cell_types)) {
      plot_df = Matrix(1, nrow = length(UMI_list), ncol = 3)
      rownames(plot_df) = UMI_list; colnames(plot_df) = c("nUMI", "first_type", "second_type")
      plot_df[,"nUMI"] = UMI_list
      type1 = common_cell_types[i]
      type2 = common_cell_types[j]
      curr_barcodes_base = meta_df$first_type == type1 & meta_df$second_type == type2
      curr_barcodes_base = curr_barcodes_base & (results_df$spot_class != "reject")
      for(first_UMI in UMI_list) {
        first_type_found = 0; second_type_found = 0;
        curr_barcodes = curr_barcodes_base & (meta_df$first_UMI == first_UMI)
        both_barcodes = curr_barcodes & (results_df$spot_class == "doublet_certain")
        singlet_barcodes = curr_barcodes & (results_df$spot_class != "doublet_certain")
        equals_class(results_df, curr_barcodes, use_class, type1, "first_type")
        first_type_total = sum(both_barcodes); second_type_total = sum(both_barcodes);
        first_type_found = first_type_found + sum(equals_class(results_df, both_barcodes, use_class, type1, "first_type") | equals_class(results_df, both_barcodes, use_class, type1, "second_type"))
        second_type_found = second_type_found + sum(equals_class(results_df, both_barcodes, use_class, type2, "first_type") | equals_class(results_df, both_barcodes, use_class, type2, "second_type"))
        sing_first = equals_class(results_df, singlet_barcodes, use_class, type1, "first_type")
        sing_second = equals_class(results_df, singlet_barcodes, use_class, type2, "first_type")
        first_type_found = first_type_found + sum(sing_first)
        second_type_found = second_type_found + sum(sing_second)
        failures = sum((!sing_first) & (!sing_second))
        first_type_total = first_type_total + sum(sing_first) + failures/2
        second_type_total = second_type_total + sum(sing_second) + failures/2
        if(first_type_total != 0)
          plot_df[as.character(first_UMI), "first_type"] = first_type_found / first_type_total
        else
          plot_df[as.character(first_UMI), "first_type"] = 1
        if(second_type_total != 0)
          plot_df[as.character(first_UMI), "second_type"] = second_type_found / second_type_total
        else
          plot_df[as.character(first_UMI), "second_type"] = 1
      }
      if(type1 == "Bergmann" && type2 == "Granule") {
        ca = 1
      }
      if(is.null(plot_df_avg))
        plot_df_avg = plot_df
      else
        plot_df_avg = plot_df + plot_df_avg
      plots[[index]] = get_id_plot(as.data.frame(plot_df)) + labs(title = paste(type1,type2))
      doub_UMI_list <- as.character(UMI_list[2:(length(UMI_list)-1)])
      square_results[type1,type2] = colMeans(plot_df[doub_UMI_list,])["first_type"]
      square_results[type2, type1] = colMeans(plot_df[doub_UMI_list,])["second_type"]
      index = index + 1
    }
  }
  if(use_class) {
    pdf(file.path(resultsdir,"doublet_classification_class.pdf"))
  } else {
    pdf(file.path(resultsdir,"doublet_classification.pdf"))
  }
  invisible(lapply(plots, print))
  dev.off()
  if(use_class) {
    pdf(file.path(resultsdir,"avg_doublet_classification_class.pdf"))
  } else {
    pdf(file.path(resultsdir,"avg_doublet_classification.pdf"))
  }
  plot_df_avg = plot_df_avg / (index - 1)
  print(get_id_plot(as.data.frame(plot_df_avg)) + labs(title = paste("First Type", "Second Type")))
  dev.off()
  return(square_results)
}
