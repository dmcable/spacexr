get_means_sds <- function(cell_type, gene, de_results_list) {
  ct_ind <- which(colnames(de_results_list[[1]]$gene_fits$mean_val) == cell_type)*2
  means <- unlist(lapply(de_results_list, function(x) x$gene_fits$mean_val[gene,cell_type]))
  sds <- unlist(lapply(de_results_list, function(x) x$gene_fits$s_mat[gene,ct_ind]))
  return(list(means = means, sds = sds))
}

get_de_pop <- function(cell_type, de_results_list) {
  ct_ind <- which(colnames(de_results_list[[1]]$gene_fits$mean_val) == cell_type)*2
  #gene_list <- intersect(intersect(rownames(de_gene_results_08$gene_fits$mean_val), rownames(de_gene_results_09$gene_fits$mean_val)),rownames(der_11$gene_fits$mean_val))
  gene_list <- Reduce(intersect, lapply(de_results_list, function(x) names(which(x$gene_fits$con_mat[,cell_type]))))
  de_pop <- matrix(0, nrow = length(gene_list), ncol = 4)
  colnames(de_pop) <- c('tau', 'log_fc_est', 'sd_est', 'Z_est')
  rownames(de_pop) <- gene_list
  for(gene in gene_list) {
    #means <- c(de_gene_results_08$gene_fits$mean_val[gene,cell_type], de_gene_results_09$gene_fits$mean_val[gene,cell_type], der_11$gene_fits$mean_val[gene,cell_type])
    means <- unlist(lapply(de_results_list, function(x) x$gene_fits$mean_val[gene,cell_type]))
    #sds <- c(de_gene_results_08$gene_fits$s_mat[gene,ct_ind], de_gene_results_09$gene_fits$s_mat[gene,ct_ind], der_11$gene_fits$s_mat[gene,ct_ind])
    sds <- unlist(lapply(de_results_list, function(x) x$gene_fits$s_mat[gene,ct_ind]))
    sds[is.na(sds)] <- 1000
    sig_p <- sqrt(estimate_var(means, sds))
    var_t <- sds^2 + sig_p^2
    var_est <- 1/sum(1 / var_t)
    mean_est <- sum(means / var_t)*var_est

    sd_est <- sqrt(var_est)
    Z_est <- mean_est / sd_est
    de_pop[gene, ] <- c(sig_p, mean_est, sd_est, Z_est)
  }
  de_pop <- as.data.frame(de_pop)
  return(de_pop)
}

one_ct_genes <- function(cell_type, myRCTD_list, de_results_list, resultsdir, cell_types_present, q_thresh = .01, p_thresh = 1, filter = T, order_gene = F, plot_results = T) {
  de_pop <- get_de_pop(cell_type, de_results_list)
  myRCTD <- myRCTD_list[[1]]
  gene_big <- Reduce(intersect, lapply(myRCTD_list,
                              function(myRCTD) get_gene_list_type_wrapper(myRCTD, cell_type, cell_types_present)))
  cell_type_means <- myRCTD@cell_type_info$info[[1]][gene_big,cell_types_present]
  cell_prop <- sweep(cell_type_means,1,apply(cell_type_means,1,max),'/')
  p_vals <- 2*(1-pnorm(abs(de_pop[gene_big,'Z_est'])))
  names(p_vals) <- gene_big
  q_vals<- p.adjust(p_vals,'BH')
  if(filter)
    gene_final <- intersect(gene_big[which(q_vals < q_thresh & p_vals < p_thresh)],
                       gene_big[which(abs(de_pop[gene_big,'log_fc_est']) > 0.4)])
  else
    gene_final <- gene_big
  final_df <- cbind(de_pop[gene_final,],cell_prop[gene_final,c(cell_type)],
                    cell_type_means[gene_final,cell_type], q_vals[gene_final])
  colnames(final_df) <- c( 'tau',    'log_fc_est'  ,   'sd_est'    ,  'Z_est'    , 'ct_prop' ,'expr' ,'q_val')
  final_df$p <- 2*(1 - pnorm(abs(final_df$Z_est)))
  L <- length(myRCTD_list)
  mean_sd_df <- matrix(0, nrow = length(gene_final), ncol = L*2)
  rownames(mean_sd_df) <- gene_final
  colnames(mean_sd_df) <- c(unlist(lapply(1:L, function(x) paste('mean', x))), unlist(lapply(1:L, function(x) paste('sd', x))))
  for(gene in gene_final) {
    m_sd <- get_means_sds(cell_type, gene, de_results_list)
    mean_sd_df[gene,] <- c(m_sd$means, m_sd$sds)
  }
  final_df <- cbind(final_df, mean_sd_df)
  if(order_gene)
    final_df <- final_df[order(gene_final), ]
  else
    final_df <- final_df[order(-abs(final_df$mean_est)),]
  #plot(log(final_df$expr,10), log(final_df$p,10))
  if(plot_results) {
    print('writing')
    write.csv(final_df,file.path(resultsdir,paste0(cell_type,'_cell_type_genes.csv')))
  }
  print('done')
  return(list(de_pop = de_pop, gene_final = gene_final, final_df = final_df))
}

estimate_var <- function(x, s) {
  return(max(var(x) - mean(s^2),0))
}
