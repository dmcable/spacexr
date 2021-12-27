get_n_sig_genes_Z_test <- function(cell_type, cur_cell_types, myRCTDde, X2,puck, fdr = .01) {
  final_df <- get_Z_test_res_ct(cell_type, cur_cell_types, myRCTDde, X2,puck, fdr = fdr)
  return(dim(final_df)[1])
}

get_Z_test_res_ct <- function(cell_type, cur_cell_types, myRCTDde, X2,puck, fdr = .01, log_fc_thresh = 0.4, p_thresh = 1) {
  cell_type_ind <- which(cur_cell_types == cell_type)
  n_cell_types <- length(cur_cell_types)
  sing_barcodes <- names(which(myRCTDde@internal_vars_de$my_beta[,cell_type_ind] > 0.9))
  if(length(sing_barcodes) < 4)
    return(0)
  region_ids <- unlist(lapply(sing_barcodes, 
                              function(barcode) which(as.logical(X2[barcode,]))))
  names(region_ids) <- sing_barcodes
  gene_list_type <- get_gene_list_type_wrapper(myRCTDde, cell_type, cell_types_present)
  p_val_list <- numeric(length(gene_list_type)); names(p_val_list) <- gene_list_type
  log_fc_list <- numeric(length(gene_list_type)); names(log_fc_list) <- gene_list_type
  log_fc_12_list <- numeric(length(gene_list_type)); names(log_fc_12_list) <- gene_list_type
  n_regions_con <- length(table(region_ids))
  if(n_regions_con < 2)
    return(0)
  region_id_names <- names(table(region_ids))
  for(gene in gene_list_type) {
    Y <- puck@counts[gene,sing_barcodes] / puck@nUMI[sing_barcodes]
    obs_df <- data.frame(unlist(lapply(region_id_names, function(region) mean(Y[sing_barcodes[region_ids == region]]))),
                         unlist(lapply(region_id_names, function(region) length(Y[sing_barcodes[region_ids == region]]))),
                         unlist(lapply(region_id_names, function(region) sd(Y[sing_barcodes[region_ids == region]]))))
    colnames(obs_df) <- c('mean','n','sd')
    obs_df$se <- obs_df$sd / sqrt(obs_df$n)
    
    ovr_best_p_val <- 1
    best_log_fc <- 0
    best_sd <- 0
    for(i1 in 1:(n_regions_con-1))
      for(i2 in (i1+1):n_regions_con) {
        if(obs_df$mean[i1] > 0 & obs_df$mean[i2] > 0) {
          log_fc_lin <- abs(obs_df$mean[i1] - obs_df$mean[i2])
          sd_cur <- sqrt(obs_df$se[i1]^2 + obs_df$se[i2]^2)
          z_score <- log_fc_lin / sd_cur
          log_fc <- abs(log(obs_df$mean[i1]+1e-20) - log(obs_df$mean[i2]+1e-20))
          if(i1 == 1 && i2 == 2)
            log_fc_12_list[gene] <- log(obs_df$mean[i1]+1e-20) - log(obs_df$mean[i2]+1e-20)
          p_val <- 2*(1-pnorm(z_score)) * choose(n_regions_con, 2)
          if(p_val < min(ovr_best_p_val, p_thresh)) {
            ovr_best_p_val <- p_val
            best_log_fc <- log_fc
            best_sd <- sd_cur
          }
        }
      }
    p_val_list[gene] <- min(1, ovr_best_p_val)
    log_fc_list[gene] <- best_log_fc
  }
  gene_list_sig <- fdr_sig_genes(gene_list_type, p_val_list, fdr)
  gene_list_sig <- gene_list_sig[log_fc_list[gene_list_sig] > log_fc_thresh]
  final_df <- cbind(p_val_list[gene_list_sig], log_fc_list[gene_list_sig], log_fc_12_list[gene_list_sig])
  rownames(final_df) <- gene_list_sig
  colnames(final_df) <- c('p', 'log_fc', 'log_fc_12')
  return(final_df)
}
cor_ct_patterns <- function(cell_type_1, cell_type_2, myRCTDde, cell_types_present, X2, gene_fits, cur_cell_types) {
  n_cell_types <- length(cur_cell_types)
  n_regions <- dim(X2)[2]
  same_genes <- intersect(get_gene_list_type_wrapper(myRCTDde, cell_type_1, cell_types_present),
                          get_gene_list_type_wrapper(myRCTDde, cell_type_2, cell_types_present))
  cell_ind_1 <- which(cur_cell_types == cell_type_1); cell_ind_2 <- which(cur_cell_types == cell_type_2)
  I_mat_ind_1 <- (1:n_regions) + (n_regions*(cell_ind_1 - 1))
  I_mat_ind_2 <- (1:n_regions) + (n_regions*(cell_ind_2 - 1))
  res_corr <- numeric(length(same_genes))
  names(res_corr) <- same_genes
  res_n <- numeric(length(same_genes))
  names(res_n) <- same_genes
  p_val_thresh <- 0.001
  for(gene in same_genes) {
    con_regions <- get_con_regions(gene_fits, gene, n_regions, cell_ind_1, n_cell_types) & 
      get_con_regions(gene_fits, gene, n_regions, cell_ind_2, n_cell_types)
    n_regions_con <- sum(con_regions)
    x1 <- gene_fits$all_vals[gene, con_regions,cell_ind_1]
    x2 <- gene_fits$all_vals[gene, con_regions,cell_ind_2]
    I_mat_ind_cur_1 <- I_mat_ind_1[con_regions]; I_mat_ind_cur_2 <- I_mat_ind_2[con_regions]
    var_vals_1 <- (gene_fits$I_mat[gene, I_mat_ind_cur_1])^2
    var_vals_2 <- (gene_fits$I_mat[gene, I_mat_ind_cur_2])^2
    ovr_best_p_val <- 1
    best_log_fc <- 0
    best_sd <- 0
    trials <- 0;
    success <- 0;
    for(i1 in 1:(n_regions_con-1))
      for(i2 in (i1+1):n_regions_con) {
        log_fc_1 <- x1[i1] - x1[i2]
        sd_cur_1 <- sqrt(var_vals_1[i1] + var_vals_1[i2])
        z_score_1 <- abs(log_fc_1) / sd_cur_1
        p_val_1 <- 2*(1-pnorm(z_score_1))
        log_fc_2 <- x2[i1] - x2[i2]
        sd_cur_2 <- sqrt(var_vals_2[i1] + var_vals_2[i2])
        z_score_2 <- abs(log_fc_2) / sd_cur_2
        p_val_2 <- 2*(1-pnorm(z_score_2))
        if(max(p_val_1, p_val_2) < p_val_thresh) {
          trials <- trials + 1
          if(sign(log_fc_1) == sign(log_fc_2))
            success <- success + 1
        }
      }
    if(trials > 0) {
      res_corr[gene] <- success / trials
    } else {
      res_corr[gene] <- -1
    }
    res_n[gene] <- trials
  }
  return(list(res_n = res_n, res_corr = res_corr))
}

cor_ct_patterns_quant <- function(cell_type_1, cell_type_2, myRCTDde, cell_types_present, X2, gene_fits, cur_cell_types) {
  n_regions <- dim(X2)[2]
  same_genes <- intersect(get_gene_list_type_wrapper(myRCTDde, cell_type_1, cell_types_present),
                          get_gene_list_type_wrapper(myRCTDde, cell_type_2, cell_types_present))
  cell_ind_1 <- which(cur_cell_types == cell_type_1); cell_ind_2 <- which(cur_cell_types == cell_type_2)
  I_mat_ind_1 <- (1:n_regions) + (n_regions*(cell_ind_1 - 1))
  I_mat_ind_2 <- (1:n_regions) + (n_regions*(cell_ind_2 - 1))
  cor_results <- matrix(0,nrow = n_regions, ncol = n_regions)
  for(region_1 in 1:(n_regions-1))
    for(region_2 in (region_1+1):n_regions) {
      con_genes <- names(which(
        gene_fits$con_all[same_genes,I_mat_ind_1[region_1]] & gene_fits$con_all[same_genes,I_mat_ind_2[region_1]]
        & gene_fits$con_all[same_genes,I_mat_ind_1[region_2]] & gene_fits$con_all[same_genes,I_mat_ind_2[region_2]]))
      diff_1 <- gene_fits$all_vals[same_genes,region_1,cell_ind_1] - gene_fits$all_vals[same_genes,region_2,cell_ind_1]
      diff_2 <- gene_fits$all_vals[same_genes,region_1,cell_ind_2] - gene_fits$all_vals[same_genes,region_2,cell_ind_2]
      cor_results[region_1,region_2] <- cor(diff_1,diff_2)
    }
  return(cor_results)
}

get_marker_by_region <- function(all_ct, cell_type, n_regions, myRCTDde, cell_types_present, gene_fits, p_val_thresh = .001, trials = 5000, log_fc_thresh = 0.4) {
  other_cell_types <- all_ct[all_ct != cell_type]
  cell_ind <- which(cell_type == cur_cell_types)
  I_mat_ind <- (1:n_regions) + (n_regions*(cell_ind - 1))
  # for cell type 1, log_fc_thresh = 0.4 -> marker genes for all regions. For ct 2, marker genes only for region 4
  # for cell type 2, log_fc_thresh = 0.4 -> marker genes for all regions
  # for cell type 4, log_fc_thresh = 0.4 -> marker genes for regions 1, 3, and 4. If lower p val thresh to 0.01, then we we get all regions
  # if we use our alternative monte carlo p value test, we get significant for all cases
  # for cell type 4, best genes are (Aurka, Syt11, Piwil1, Smg9)
  ct_genes <- get_gene_list_type_wrapper(myRCTDde, cell_type, cell_types_present)
  marker_ratio <- log(myRCTDde@cell_type_info$info[[1]][ct_genes,cell_type] / 
                        apply(myRCTDde@cell_type_info$info[[1]][ct_genes,other_cell_types] ,1,max)) 
  marker_genes <- names(which(marker_ratio > log_fc_thresh))
  Z_thresh_vals <- unlist(lapply(marker_genes, function(gene) {
    print(gene)
    z <- matrix(0, nrow = trials, ncol = n_regions)
    z[] <- rnorm(n_regions*trials)
    sd_vals <- diag(gene_fits$I_val[[which(rownames(gene_fits$I_mat) == gene)]][I_mat_ind,I_mat_ind])
    z <- sweep(z, 2, sd_vals, '*')
    diff_sim <- apply(z,1,max) - apply(z,1,function(x) x[order(-x)][2])
    first_region_sim <- apply(z,1,which.max)
    second_region_sim <- apply(z,1,function(x) order(-x)[2])
    Z_sim <- diff_sim / sqrt(sd_vals[first_region_sim]^2 + sd_vals[second_region_sim]^2)
    Z_thresh <- Z_sim[order(-Z_sim)[round(trials * p_val_thresh)]]
    return(Z_thresh)
  }))
  names(Z_thresh_vals) <- marker_genes
  
  #table(gene_fits$con_all[marker_genes,I_mat_ind]) #all converged for all cell types wooo!!!
  diff <- apply(gene_fits$all_vals[marker_genes, , cell_ind],1,max) - 
    apply(gene_fits$all_vals[marker_genes, , cell_ind],1,function(x) x[order(-x)][2])
  first_region <- apply(gene_fits$all_vals[marker_genes, , cell_ind],1,which.max)
  second_region <- apply(gene_fits$all_vals[marker_genes, , cell_ind],1,function(x) order(-x)[2])
  Z_scores <- unlist(lapply(marker_genes, function(gene) diff[gene] / sqrt(sum(gene_fits$I_mat[gene,I_mat_ind[c(first_region[gene],second_region[gene])]]^2))))
  p_vals <- 2*(1-pnorm(Z_scores))
  #sig_genes <- names(which(p_vals < p_thresh & diff > log_fc_thresh))
  sig_genes <- names(which(Z_scores > Z_thresh_vals & diff > log_fc_thresh))
  table(first_region[sig_genes])
  hist(p_vals)
  info_df <- data.frame(sig_genes, first_region[sig_genes], Z_scores[sig_genes], Z_thresh_vals[sig_genes], 
                        p_vals[sig_genes], diff[sig_genes], second_region[sig_genes], marker_ratio[sig_genes],
                        log(myRCTDde@cell_type_info$info[[1]][sig_genes,cell_type]))
  colnames(info_df) <- c('gene','first_region','Z', 'Z_thresh','p','diff','second_region','marker_ratio', 'expr')
  return(info_df)
}

make_ct_region_plot <- function(cur_cell_types, cell_type, region, my_beta, X2, gene_thresh, puck, gene, sing_thresh = 0.8) {
  ct_ind <- which(cur_cell_types == cell_type)
  prop_list <- my_beta[intersect(rownames(my_beta),names(which(X2[,region] == 1))),ct_ind]
  ct_reg_barc <- names(which(prop_list > sing_thresh))
  MULT <- 500
  if(length(gene) > 1)
    Y_plot <- MULT*colSums(puck@counts[gene,])/puck@nUMI
  else
    Y_plot <- MULT*puck@counts[gene,]/puck@nUMI
  my_title = paste(gene, gene_thresh)
  barc_plot <- c(names(which(prop_list > sing_thresh | prop_list < 0.2)),
                 rownames(my_beta[intersect(rownames(my_beta),names(which(X2[,region] == 0))),]))# exclude the middle ones
  my_class <- rep(0,length(barc_plot)); names(my_class) <- barc_plot
  my_class[!(barc_plot %in% ct_reg_barc) & (Y_plot[barc_plot] <= gene_thresh)] <- 1
  my_class[!(barc_plot %in% ct_reg_barc) & (Y_plot[barc_plot] > gene_thresh)] <- 3
  my_class[(barc_plot %in% ct_reg_barc) & (Y_plot[barc_plot] <= gene_thresh)] <- 2
  my_class[(barc_plot %in% ct_reg_barc) & (Y_plot[barc_plot] > gene_thresh)] <- 4
  p3 <- plot_class(puck, barc_plot[order(my_class[barc_plot])], factor(my_class), title = my_title)
  suppressMessages(p3 <- p3 + scale_color_manual(values=c("#CCE2EF","#F6DECC","#0072B2","#D55E00")))
  return(p3)
}

get_Z_test_res <- function(cur_cell_types, myRCTDde, X2, puck) {
  Z_test_results <- unlist(lapply(cur_cell_types,
                                  function(cell_type) get_n_sig_genes_Z_test(cell_type, cur_cell_types, myRCTDde, X2,puck, fdr = .01)))
  names(Z_test_results) <- cur_cell_types
  #saveRDS(Z_test_results, file = file.path(datadir,'z_test.rds'))
  RCTDde_test <- unlist(lapply(cur_cell_types, 
                               function(cell_type) length(rownames(myRCTDde@de_results$res_gene_list[[cell_type]]))))
  names(RCTDde_test) <- cur_cell_types
  n_singlets <- unlist(lapply(1:length(cur_cell_types),
                              function(cell_type_ind) sum(myRCTDde@internal_vars_de$my_beta[,cell_type_ind] > 0.9)))
  names(n_singlets) <- cur_cell_types
  n_total <- colSums(myRCTDde@internal_vars_de$my_beta[,cur_cell_types])
  df1 <- data.frame(cur_cell_types, n_singlets, n_total, Z_test_results,'Z')
  df2 <- data.frame(cur_cell_types, n_singlets, n_total, RCTDde_test, 'RCTD')
  colnames(df1)[4:5] <- c('count','method'); colnames(df2)[4:5] <- c('count','method')
  plot_df <- rbind(df1,df2)
  plot_df[plot_df$method == 'Z','n_total'] <- plot_df[plot_df$method == 'Z','n_singlets']
  return(plot_df)
}

is.maxcyclic.thresh.large <- function(x, D_thresh = -0.25) {
  my_ind <- (which.max(x) + 1) %% 4 + 1
  return(x[my_ind] - x[order(x)][2] <= D_thresh)
}

breathing.room <- function(x) {
  x <- x[order(x)]
  return(min(x[2] - x[1], x[4] - x[3]))
}

shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}
norm_vec <- function(x) {
  x <- x - min(x)
  x <- x / max(x)
  max_ind <- which.max(x)
  return(shifter(x, max_ind - 2))
}

is.cyclic <- function(x) {
  if(x[2] < x[1]) {
    return(!(x[3] > x[2] & x[4] < x[3]))
  } else {
    return(!(x[3] < x[2] & x[4] > x[3]))
  }
}

is.maxcyclic <- function(x) {
  return(abs(which.min(x) - which.max(x)) == 2)
}

D_thresh <- 0.25

is.maxcyclic.thresh <- function(x) {
  my_ind <- (which.max(x) + 1) %% 4 + 1
  return(x[my_ind] - min(x) <= D_thresh)
}

maxcyclic.thresh.prob <- function(x) {
  return(min(sum(x - min(x) < D_thresh),3)/3)
}
