normalize_de_estimates <- function(myRCTD, normalize_expr, remove_junk = T, param_position = 2) {
  if(remove_junk) {
    all_genes <- rownames(myRCTD@de_results$gene_fits$con_mat)
    junk_genes <- all_genes[c(grep("MT-",all_genes),grep("RPS",all_genes),
                              grep("RPL",all_genes))]
    junk_genes <- c(junk_genes, 'MALAT1')
  } else {
    junk_genes <- c()
  }
  for(cell_type in myRCTD@internal_vars_de$cell_types) {
    con_genes_all <- names(which(myRCTD@de_results$gene_fits$con_mat[,cell_type]))
    con_genes <- con_genes_all[!(tolower(con_genes_all) %in% tolower(junk_genes))]
    if(!('mean_val_cor' %in% names(myRCTD@de_results$gene_fits)))
      myRCTD@de_results$gene_fits$mean_val_cor <- list()
    if(normalize_expr) {
      reg_1 <- myRCTD@de_results$gene_fits$all_vals[con_genes_all,1,cell_type]
      reg_2 <- reg_1 + myRCTD@de_results$gene_fits$all_vals[con_genes_all,param_position,cell_type]
      reg_1_cor <- reg_1 - log(sum(exp(reg_1)[con_genes]))
      reg_2_cor <- reg_2 - log(sum(exp(reg_2)[con_genes]))
      mean_val_cor <- reg_2_cor - reg_1_cor
      myRCTD@de_results$gene_fits$mean_val_cor[[cell_type]] <- mean_val_cor
    } else {
      myRCTD@de_results$gene_fits$mean_val_cor[[cell_type]] <- myRCTD@de_results$gene_fits$mean_val[,cell_type]
    }
  }
  return(myRCTD)
}
