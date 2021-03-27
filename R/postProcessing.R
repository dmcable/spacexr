#functions for processing data after RCTD is fit to the data

# Collects RCTD results
gather_results <- function(RCTD, results) {
  cell_type_names = RCTD@cell_type_info$renorm[[2]]
  barcodes <- colnames(RCTD@spatialRNA@counts)
  N <- length(results)
  weights = Matrix(0, nrow = N, ncol = length(cell_type_names))
  weights_doublet = Matrix(0, nrow = N, ncol = 2)
  rownames(weights) = barcodes; rownames(weights_doublet) = barcodes
  colnames(weights) = cell_type_names; colnames(weights_doublet) = c('first_type', 'second_type')
  empty_cell_types = factor(character(N),levels = cell_type_names)
  spot_levels <- c("reject", "singlet", "doublet_certain", "doublet_uncertain")
  results_df <- data.frame(spot_class = factor(character(N),levels=spot_levels),
                           first_type = empty_cell_types, second_type = empty_cell_types,
                           first_class = logical(N), second_class = logical(N),
                           min_score = numeric(N), singlet_score = numeric(N),
                           conv_all = logical(N), conv_doublet = logical(N))
  score_mat <- list()
  for(i in 1:N) {
    if(i %% 1000 == 0)
      print(paste("gather_results: finished",i))
    weights_doublet[i,] = results[[i]]$doublet_weights
    weights[i,] = results[[i]]$all_weights
    results_df[i, "spot_class"] = results[[i]]$spot_class
    results_df[i, "first_type"] = results[[i]]$first_type
    results_df[i, "second_type"] = results[[i]]$second_type
    results_df[i, "first_class"] = results[[i]]$first_class
    results_df[i, "second_class"] = results[[i]]$second_class
    results_df[i, "min_score"] = results[[i]]$min_score
    results_df[i, "singlet_score"] = results[[i]]$singlet_score
    results_df[i, "conv_all"] = results[[i]]$conv_all
    results_df[i, "conv_doublet"] = results[[i]]$conv_doublet
    score_mat[[i]] <- results[[i]]$score_mat
  }
  rownames(results_df) = barcodes
  RCTD@results <- list(results_df = results_df, weights = weights, weights_doublet = weights_doublet, score_mat = score_mat)
  return(RCTD)
}

