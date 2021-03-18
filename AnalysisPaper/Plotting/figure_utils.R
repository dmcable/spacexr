library(gplots)
library("RColorBrewer")

plot_heat_map <- function(data, file_loc = NULL, save.file = F, normalize = T, dendrogram = 'none') {
  if(normalize) {
    norm_conf = sweep(data, 2, colSums(data), '/')
    norm_conf[norm_conf < 0] = 0
  } else
    norm_conf = data
  if(save.file) {
    png(file = file_loc)
    heatmap.2(norm_conf,col=rev(heat.colors(100)), breaks=seq(0,1,0.01),scale = "none",trace="none", Rowv=FALSE, Colv=FALSE,dendrogram=dendrogram)
    dev.off()
  } else {
    heatmap.2(norm_conf,col=rev(heat.colors(100)), breaks=seq(0,1,0.01), scale = "none",trace="none", Rowv=FALSE, Colv=FALSE,dendrogram=dendrogram)
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
      if(type1 == "Granule" && type2 == "UBCs") {
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

get_decompose_plots <- function(meta_df, type1, type2, weights_doublet, meta_data, iv, expect1, expect2, var, first_beads, second_beads, beads) {
  cur_df <- meta_df[meta_df$first_type == type1 & meta_df$second_type == type2,]
  plot_df <- aggregate(weights_doublet[rownames(cur_df),1], list(cur_df$first_UMI/meta_data$UMI_tot), mean)
  colnames(plot_df) = c('proportion','type1_avg')
  plot_df[,"st_dev"] = aggregate(weights_doublet[rownames(cur_df),1], list(cur_df$first_UMI/meta_data$UMI_tot), sd)$x*1.96
  weight_plot<- ggplot2::ggplot(plot_df, ggplot2::aes(x=proportion, y=type1_avg, colour = "type1_avg")) +
    ggplot2::geom_line() +
    ggplot2::geom_point()+
    ggplot2::geom_line(ggplot2::aes(y=proportion,colour = "proportion")) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=type1_avg-st_dev, ymax=type1_avg+st_dev), width=.05,
                           position=ggplot2::position_dodge(0.05)) + ggplot2::labs(title=paste(type1,type2)) + theme_classic()
  plot_df_weight <- plot_df
  #RMSE
  RMSE <- sqrt(mean((weights_doublet[rownames(cur_df),1] - cur_df$first_UMI/meta_data$UMI_tot)^2))

  #
  v = as.vector(var[rownames(cur_df),])
  x = as.vector(expect1[rownames(cur_df),])
  y = as.vector(t(first_beads[iv$gene_list,rownames(cur_df)]))
  z = x + as.vector(expect2[rownames(cur_df),]) #total
  #cor(x[z>0]/z[z>0],y[z>0]/z[z>0])

  pred_por = x[z>0]/z[z>0]
  true_por = y[z>0]/z[z>0]
  my_var = v[z > 0]/z[z > 0]^2
  #true vs expected mean squared error.
  num_bins = 10
  mean_res = numeric(num_bins)
  mean_pred = numeric(num_bins)
  exp_mse = numeric(num_bins)
  true_mse = numeric(num_bins)
  null_mse = numeric(num_bins)
  for(ind in 1:num_bins) {
    cur_ind = pred_por > (ind-1)/(num_bins) & pred_por < (ind)/(num_bins)
    mean_res[ind] = mean(true_por[cur_ind]); mean_pred[ind] = mean(pred_por[cur_ind])
    true_mse[ind] = mean((true_por[cur_ind] - pred_por[cur_ind])^2)
    exp_mse[ind] = mean(my_var[cur_ind])
    null_mse[ind] = mean((true_por[cur_ind] - 0.5)^2)
  }
  err_df <- data.frame(mean_pred, exp_mse, true_mse)
  err_df <- data.frame(mean_pred, sqrt(exp_mse), sqrt(true_mse), sqrt(null_mse))
  colnames(err_df) = c('mean_pred', 'exp_mse', 'true_mse', 'null_mse')
  err_plot<- ggplot2::ggplot(err_df, ggplot2::aes(x=mean_pred, y=exp_mse, colour = "exp_mse")) +
    ggplot2::geom_line() +
    ggplot2::geom_point()+
    ggplot2::geom_line(ggplot2::aes(y=true_mse,colour = "true_mse")) +
    ggplot2::geom_line(ggplot2::aes(y=null_mse,colour = "null_mse")) + ggplot2::ylim(c(0,1)) + ggplot2::labs(title=paste(type1,type2))+ theme_classic()
  plot_df = data.frame(mean_pred, mean_res)
  bias_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x=mean_pred,y=mean_res)) + geom_line() + ggplot2::labs(title=paste(type1,type2))+ theme_classic()
  #overall predictions distribution
  plot_df <- as.data.frame(pred_por)
  hist_plot <- ggplot(plot_df, aes(x=pred_por)) + geom_histogram(bins=30) + ggplot2::labs(title=paste(type1,type2))+ theme_classic()
  overall_mse <- mean((true_por - pred_por)^2)
  null_mse <- mean((true_por - 0.5)^2)
  samp_mse <- mean(my_var)
  R2 <- (null_mse - overall_mse)/(null_mse - samp_mse)
  #DE genes
  epsilon = 1e-9
  de_score <- log(iv$cell_type_info[[1]][,type1] + epsilon) - log(iv$cell_type_info[[1]][,type2] + epsilon)
  names(de_score) <- iv$gene_list
  cur_df_mid <- meta_df[meta_df$first_type == type1 & meta_df$second_type == type2 & meta_df$first_UMI >= meta_data$UMI_tot*0.4 & meta_df$first_UMI <= meta_data$UMI_tot*0.6,]
  de_score <- de_score[colSums(beads[rownames(cur_df_mid), ] > 0) >= 20] # highly expressed
  N_genes <- 10
  de_thresh = 3
  de_one <- de_score[de_score > de_thresh]
  if(length(de_one) > 0)
    de_one <- tail(de_one[order(de_one)],N_genes)
  else
    de_one <- tail(de_score[order(de_score)],1)
  de_two <- de_score[de_score < -de_thresh]
  if(length(de_two) > 0)
    de_two <- head(de_two[order(de_two)],N_genes)
  else
    de_two <- head(de_score[order(de_score)],1)
  if(length(de_one))
  de_gene_names <- c(names(de_one),names(de_two))
  e_por <- colMeans(expect1[rownames(cur_df_mid), de_gene_names] / (beads[rownames(cur_df_mid), de_gene_names]), na.rm=T)
  t_por <- colMeans((t(first_beads[de_gene_names, rownames(cur_df_mid)]) / (beads[rownames(cur_df_mid), de_gene_names])), na.rm=T)
  plot_df <- data.frame(e_por, t_por)
  score_one <- mean(abs(e_por[names(de_one)] - t_por[names(de_one)])/ t_por[names(de_one)])
  score_two <- mean(abs(e_por[names(de_two)] - t_por[names(de_two)])/ (1-t_por[names(de_two)]))
  de_gene_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x=e_por,y=t_por)) + geom_point() + ggplot2::labs(title=paste(type1,type2))
  plot_df <- data.frame(e_por, t_por, factor(de_gene_names, levels = de_gene_names))
  colnames(plot_df) = c('predicted', 'true', 'gene')
  de_plot_df <- plot_df
  de_ind_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x=gene,y=predicted)) + geom_point(color='darkblue') + geom_point(aes(x = gene, y = true), color='darkred') +
    ggplot2::labs(title=paste(type1,type2))+ theme(axis.text.x = element_text(hjust = 1, angle = 45))+ theme_classic()
  return(list(weight_plot = weight_plot, RMSE = RMSE, bias_plot = bias_plot, err_plot = err_plot, hist_plot = hist_plot, R2 = R2, de_gene_plot = de_gene_plot, de_ind_plot = de_ind_plot,
              plot_df_weight = plot_df_weight, de_scores = c(score_one, score_two), de_plot_df = de_plot_df))
}
