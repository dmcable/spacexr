plot_to_results <- function(n_bins, plot_df, region_range = 0:1) {
  cutoffs <- (0:n_bins) / n_bins
  results_df <- matrix(0, nrow=n_bins*length(region_range),6)
  for(region in region_range) {
    for(bin in 1:n_bins) {
      in_bin <- which(plot_df$berg_weight <= cutoffs[bin+1] & plot_df$berg_weight >= cutoffs[bin] & plot_df$region == region)
      mean_Y <- (mean(plot_df$Y[in_bin]))
      mean_pred <- mean(plot_df$pred[in_bin])
      mean_weight <- mean(plot_df$berg_weight[in_bin])
      sd_Y <- sd(plot_df$Y[in_bin]) / sqrt(length(in_bin))
      results_df[(region - min(region_range))*n_bins + bin,] <- c(region, bin, mean_Y,mean_pred,mean_weight,sd_Y)
    }
  }
  results_df <- data.frame(results_df)
  colnames(results_df) <- c('region', 'bin', 'Y', 'pred', 'berg_weight','se')
  results_df$region <- factor(results_df$region)
  return(results_df)
}

generate_results_df <- function(gene, cur_cell_types, n_bins, RCTDde_list, de_results_list, cell_types, merge_reps = T) {
  plot_df_list <- list()
  for(rep_ind in 1:length(RCTDde_list)) {
    RCTDde <- RCTDde_list[[rep_ind]]
    gene_fits <- de_results_list[[rep_ind]]$gene_fits
    plot_df <- get_quant_df(RCTDde, gene_fits, cell_types, cur_cell_types, gene)
    plot_df_list[[rep_ind]] <- plot_df
  }
  if(merge_reps) {
    plot_df <- do.call("rbind", plot_df_list)
    results_df <- plot_to_results(n_bins, plot_df)
  } else {
    results_df <- lapply(1:length(RCTDde_list), function(rep_ind) cbind(rep_ind, plot_to_results(n_bins,plot_df_list[[rep_ind]])))
    results_df <- data.frame(do.call("rbind",results_df))
    colnames(results_df)[1] <- 'rep'
    results_df$rep <- factor(results_df$rep)
  }
  return(results_df)
}

plot_single_gene_reps <- function(gene, cur_cell_types, RCTDde_list, de_results_list, cell_types) {
  n_bins <- 1
  results_df <- generate_results_df(gene, cur_cell_types, n_bins, RCTDde_list, de_results_list, cell_types, merge_reps = F)
  jitter_obj <- position_jitter(width = 0.05, height = 0, seed = 123)
  ggplot(results_df) + geom_point(aes(x = region, y = Y, color = rep), position = jitter_obj) + 
    geom_line(aes(x = region, y = pred, color = rep, group = rep), position = jitter_obj) + 
    geom_errorbar(aes(x = region, y = Y, ymin = Y - 1.96*se, ymax = Y + 1.96*se, width = .05, color = rep), position = jitter_obj) + 
    theme_classic()
}