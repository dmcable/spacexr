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


#' @export
get_decomposed_data_full_doublet <- function(gene_list, puck, weights, ct_info) {
  first_DGE <- Matrix(0, nrow = dim(weights)[1], ncol = length(gene_list))
  second_DGE <- Matrix(0, nrow = dim(weights)[1], ncol = length(gene_list))
  rownames(first_DGE) = rownames(weights); rownames(second_DGE) = rownames(weights)
  colnames(first_DGE) = gene_list; colnames(second_DGE) = gene_list
  for(ind in 1:dim(weights)[1]) {
    barcode = rownames(weights)[ind]
    doub_res <- decompose_doublet_fast(puck@counts[gene_list,barcode], weights[barcode,], gene_list, ct_info, colnames(weights)[1],colnames(weights)[2])
    first_DGE[barcode,] <- doub_res$expect_1; second_DGE[barcode,] <- doub_res$expect_2
  }
  norm1 <- sweep(first_DGE, 1, weights[rownames(weights),1] * puck@nUMI[rownames(first_DGE)], '/')
  norm2 <- sweep(second_DGE, 1, weights[rownames(weights),2] * puck@nUMI[rownames(second_DGE)], '/')
  all_DGE <- rbind(norm1, norm2)
  cell_type_labels <- unlist(list(rep(colnames(weights)[1], dim(weights)[1]), rep(colnames(weights)[2], dim(weights)[1])))
  nUMI <- c(weights[rownames(weights),1] *puck@nUMI[rownames(first_DGE)], weights[rownames(weights),2]*puck@nUMI[rownames(second_DGE)])
  rownames(all_DGE) = 1:dim(all_DGE)[1]; names(cell_type_labels) <- 1:dim(all_DGE)[1]; names(nUMI) <- 1:dim(all_DGE)[1]
  ref_d <- Reference(t(all_DGE), factor(cell_type_labels), nUMI, require_int = F)
  return(ref_d)
}

#decompose a doublet into two cells
decompose_doublet_fast <- function(bead, weights, gene_list, cell_type_info, type1, type2) {
  N_genes = length(gene_list)
  expect_1 = vector(mode="numeric",length = N_genes)
  expect_2 = vector(mode="numeric",length = N_genes)
  variance = vector(mode="numeric",length = N_genes)
  names(expect_1) = gene_list; names(expect_2) = gene_list; names(variance) = gene_list
  epsilon = 1e-10
  for(ind in which(bead > 0)) {
    gene = gene_list[ind]
    denom = weights[1] * cell_type_info[[1]][gene,type1] + weights[2] * cell_type_info[[1]][gene,type2] + epsilon
    posterior_1 = (weights[1] * cell_type_info[[1]][gene,type1] + epsilon / 2) / denom
    expect_1[[ind]] = posterior_1 * bead[gene]
    expect_2[[ind]] = bead[gene] - posterior_1 * bead[gene]
    variance[[ind]] = posterior_1 * bead[gene] * (1 - posterior_1)
  }
  return(list(expect_1 = expect_1, expect_2 = expect_2, variance = variance))
}

