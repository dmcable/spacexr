get_means_sds <- function(cell_type, gene, de_results_list, params_to_test) {
  de_results <- de_results_list[[1]]
  ct_ind <- which(colnames(de_results$gene_fits$mean_val) == cell_type)
  L <- dim(de_results$gene_fits$s_mat)[2] / dim(de_results$gene_fits$mean_val)[2]
  ct_ind <- L*(ct_ind - 1) + params_to_test
  means <- rep(0, length(de_results_list))
  sds <- rep(-1, length(de_results_list))
  con <- unlist(lapply(de_results_list, function(x)
    ifelse(gene %in% rownames(x$gene_fits$con_mat),
           x$gene_fits$con_mat[gene,cell_type], FALSE)))
  means[con] <- unlist(lapply(de_results_list[con], function(x) x$gene_fits$mean_val_cor[[cell_type]][gene]))
  sds[con] <- unlist(lapply(de_results_list[con], function(x) x$gene_fits$s_mat[gene,ct_ind]))
  return(list(means = means, sds = sds))
}

get_de_pop <- function(cell_type, de_results_list, cell_prop, params_to_test, use.groups = F, group_ids = NULL,
                       MIN.CONV.REPLICATES = 2, MIN.CONV.GROUPS = 2, CT.PROP = 0.5, S.MAX = 4) {
  if(!use.groups)
    group_ids <- NULL
  de_results <- de_results_list[[1]]
  ct_ind <- which(colnames(de_results$gene_fits$mean_val) == cell_type)
  L <- dim(de_results$gene_fits$s_mat)[2] / dim(de_results$gene_fits$mean_val)[2]
  ct_ind <- L*(ct_ind - 1) + params_to_test
  gene_list <- Reduce(union, lapply(de_results_list, function(x) names(which(x$gene_fits$con_mat[,cell_type]))))
  gene_list <- intersect(gene_list, rownames(cell_prop)[(which(cell_prop[,cell_type] >= CT.PROP))])
  if(!use.groups) {
    de_pop <- matrix(0, nrow = length(gene_list), ncol = 5)
    colnames(de_pop) <- c('tau', 'log_fc_est', 'sd_est', 'Z_est', 'p_cross')
  } else {
    group_names <- unique(group_ids)
    n_groups <- length(group_names)
    de_pop <- matrix(0, nrow = length(gene_list), ncol = 6 + 2*n_groups)
    colnames(de_pop) <- c('tau', 'log_fc_est', 'sd_est', 'Z_est', 'p_cross','delta',
                          unlist(lapply(group_names, function(x) paste0(x,'_group_mean'))),
                          unlist(lapply(group_names, function(x) paste0(x,'_group_sd'))))
  }
  rownames(de_pop) <- gene_list
  ii <- 1
  for(gene in gene_list) {
    ii <- ii + 1
    if(ii %% 1000 == 0)
      message(paste('get_de_pop: testing gene,', gene,', of index:', ii))
    #con <- unlist(lapply(de_results_list, function(x) gene %in%
    #         names(which(x$gene_fits$con_mat[,cell_type]))))
    check_con <- function(x) {
      ifelse(gene %in% rownames(x$gene_fits$con_mat),
             x$gene_fits$con_mat[gene,cell_type] && !is.na(x$gene_fits$s_mat[gene, ct_ind]) &&
               (x$gene_fits$s_mat[gene, ct_ind] < S.MAX), FALSE)
    }
    con <- unlist(lapply(de_results_list, check_con))
    if(use.groups)
      con <- unname(con & table(group_ids[con])[as.character(group_ids)] >= 2)
    used_groups <- names(table(group_ids[con]))
    if(sum(con) < MIN.CONV.REPLICATES || (use.groups && length(used_groups) < MIN.CONV.GROUPS)) {
      if(use.groups)
        de_pop[gene, ] <- c(-1, 0, 0, 0, 0, 0, rep(0, n_groups), rep(-1, n_groups))
      else
        de_pop[gene, ] <- c(-1, 0, 0, 0, 0)
    } else {
      means <- unlist(lapply(de_results_list[con], function(x) x$gene_fits$mean_val_cor[[cell_type]][gene]))
      sds <- unlist(lapply(de_results_list[con], function(x) x$gene_fits$s_mat[gene,ct_ind]))
      sds[is.na(sds)] <- 1000
      if(is.null(group_ids))
        gid <- NULL
      else
        gid <- group_ids[con]
      sig_p <- estimate_tau_group(means, sds, group_ids = gid)
      var_t <- sds^2 + sig_p^2
      if(!use.groups) {
        var_est <- 1/sum(1 / var_t)
        mean_est <- sum(means / var_t)*var_est
        p_cross <- get_p_qf(means, sds)
      } else {
        S2 <- 1/(aggregate(1/var_t,list(group_ids[con]),sum)$x)
        E <- (aggregate(means/var_t,list(group_ids[con]),sum)$x)*S2
        Delta <- estimate_tau_group(E, sqrt(S2))
        var_T <- (Delta^2 + S2)
        var_est <- 1/sum(1/var_T) # A_var
        mean_est <- sum(E / var_T) * var_est # A_est
        p_cross <- get_p_qf(E, sqrt(S2))
        E_all <- rep(0, n_groups); s_all <- rep(-1, n_groups)
        names(E_all) <- group_names; names(s_all) <- group_names
        E_all[used_groups] <- E; s_all[used_groups] <- sqrt(S2)
      }
      sd_est <- sqrt(var_est)
      Z_est <- mean_est / sd_est
      if(use.groups)
        de_pop[gene, ] <- c(sig_p, mean_est, sd_est, Z_est, p_cross, Delta, E_all, s_all)
      else
        de_pop[gene, ] <- c(sig_p, mean_est, sd_est, Z_est, p_cross)
    }
  }
  de_pop <- as.data.frame(de_pop)
  return(de_pop)
}

one_ct_genes <- function(cell_type, myRCTD_list, de_results_list, resultsdir, cell_types_present, params_to_test,
                         q_thresh = .01, p_thresh = 1, filter = T, order_gene = F, plot_results = T,
                         use.groups = F, group_ids = NULL, MIN.CONV.REPLICATES = 2,
                         MIN.CONV.GROUPS = 2, CT.PROP = 0.5, log_fc_thresh = 0.4, normalize_expr = F) {
  print(paste0('one_ct_genes: population inference on cell type, ', cell_type))
  myRCTD <- myRCTD_list[[1]]
  cell_type_means <- myRCTD@cell_type_info$info[[1]][,cell_types_present]
  cell_prop <- sweep(cell_type_means,1,apply(cell_type_means,1,max),'/')
  de_pop <- get_de_pop(cell_type, de_results_list, cell_prop, params_to_test, use.groups = use.groups, group_ids = group_ids,
                       MIN.CONV.REPLICATES = MIN.CONV.REPLICATES, MIN.CONV.GROUPS = MIN.CONV.GROUPS, CT.PROP = CT.PROP)
  gene_big <- rownames(de_pop)[which(de_pop$tau >= 0)]
  p_vals <- 2*(pnorm(-abs(de_pop[gene_big,'Z_est'])))
  names(p_vals) <- gene_big
  q_vals<- p.adjust(p_vals,'BH')
  if(filter)
    gene_final <- intersect(gene_big[which(q_vals < q_thresh & p_vals < p_thresh)],
                       gene_big[which(abs(de_pop[gene_big,'log_fc_est']) > log_fc_thresh)])
  else
    gene_final <- gene_big
  gene_df <- cbind(de_pop[gene_big,],cell_prop[gene_big,c(cell_type)],
                    cell_type_means[gene_big,cell_type], q_vals[gene_big])
  colnames(gene_df) <- c(colnames(de_pop), 'ct_prop' ,'expr' ,'q_val')
  gene_df$p <- 2*(pnorm(-abs(gene_df$Z_est)))
  final_df <- gene_df[gene_final, ]
  L <- length(myRCTD_list)
  mean_sd_df <- matrix(0, nrow = length(gene_final), ncol = L*2)
  rownames(mean_sd_df) <- gene_final
  colnames(mean_sd_df) <- c(unlist(lapply(1:L, function(x) paste('mean', x))), unlist(lapply(1:L, function(x) paste('sd', x))))
  for(gene in gene_final) {
    m_sd <- get_means_sds(cell_type, gene, de_results_list, params_to_test)
    mean_sd_df[gene,] <- c(m_sd$means, m_sd$sds)
  }
  final_df <- cbind(final_df, mean_sd_df)
  if(length(gene_final) > 1)
    if(order_gene)
      final_df <- final_df[order(gene_final), ]
    else
      final_df <- final_df[order(-abs(final_df$log_fc_est)),]
  #plot(log(final_df$expr,10), log(final_df$p,10))
  if(plot_results) {
    print('writing')
    write.csv(final_df,file.path(resultsdir,paste0(cell_type,'_cell_type_genes.csv')))
  }
  print('done')
  return(list(de_pop = gene_df, gene_final = gene_final, final_df = final_df))
}

get_p_qf <- function(x, se, delta = 0) {
  S <- diag(se^2+delta^2)
  n <- length(x)
  A <- matrix(-1/(n*(n-1)), nrow = n, ncol =n)
  diag(A) <- 1/n
  A
  AS <- A %*% S
  lambda <- eigen(AS)$values
  max(CompQuadForm::imhof(var(x), pmax(lambda, 10^(-8)), epsabs = 10^(-8), epsrel = 10^(-8))$Qq, 0)
}

estimate_tau_group <- function(x, s, n.iter = 20, epsilon = .001, group_ids = NULL) {
  if(is.null(group_ids))
    return(estimate_tau(x,s))
  else {
    return(mean(unlist(lapply(unique(group_ids),function(val)
      estimate_tau(x[group_ids == val], s[group_ids == val])))))
  }
}

estimate_tau <- function(x, s, n.iter = 100, epsilon = .001) {
  k <- length(x)
  tau <- 0
  for(i in 1:n.iter) {
    w <- 1/(s^2 + tau^2)
    u <- sum(x*w)/sum(w)
    Q <- sum((x - u)^2 * w)
    tau_new <- sqrt(max((Q - (k - 1)) / (sum(w) - sum(w^2)/sum(w)), 0))
    tau_last <- tau
    tau <- (tau_new + tau_last) / 2
    if(abs(tau - tau_last) < epsilon) {
      return(tau)
    }
  }
  warning('estimate_tau: not converged')
  return(tau)
}
