library(Matrix)
library(RCTD)
index = 1
for(barcode in doublet_barcodes) {
  print(index)
  index = index + 1
  score_best <- 100000
  score_second <- 100000
  best_type <- NULL
  for (type in inter_names) {
    score <- decompose_sparse(cell_type_info[[1]], gene_list, puck@nUMI[barcode], puck@counts[gene_list,barcode], type1=type, type2=as.character(second_type_list[barcode]), score_mode = T, constrain = F)
    if(score < score_best) {
      score_second <- score_best
      score_best <- score
      best_type <- type
    } else if(score < score_second) {
      score_second <- score
    }
    inter_df[barcode,type] <- score
  }
  inter_df[barcode,"confident"] <- (score_second - score_best) > log_l_thresh
  inter_df[barcode,"score_diff"] <- (score_second - score_best)
  inter_df[barcode,"best_type"] <- best_type
}
saveRDS(inter_df, file.path(resultsdir,'inter_df_process_doublet.RDS'))
