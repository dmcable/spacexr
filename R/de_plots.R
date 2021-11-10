#'
#' @export
make_all_de_plots <- function(myRCTD, datadir, cell_types_present = NULL) {
  if(!dir.exists(datadir))
    dir.create(datadir)
  make_de_plots_quant(myRCTD, datadir)
  write_de_summary(myRCTD, datadir)
  if(myRCTD@internal_vars_de$test_mode == 'direct') {
    post_filter_genes(myRCTD, datadir, cell_types_present = cell_types_present)
    make_de_plots_spatial(myRCTD, datadir)
  } else {
    make_de_plots_regions(myRCTD, datadir) #multi
  }
}

make_de_plots_quant <- function(myRCTD, datadir) {
  if(!dir.exists(file.path(datadir,'de_plots_quant')))
    dir.create(file.path(datadir,'de_plots_quant'))
  for(cell_type in myRCTD@internal_vars_de$cell_types) {
    res_genes <- myRCTD@de_results$res_gene_list[[cell_type]]
    if(myRCTD@internal_vars_de$test_mode == 'direct') {
      plot_sig_genes_quant(cell_type,myRCTD@internal_vars_de$all_barc,myRCTD@internal_vars_de$my_beta,myRCTD@originalSpatialRNA,
                             res_genes,datadir,myRCTD@internal_vars_de$X2, myRCTD@de_results$gene_fits, 
                             which(myRCTD@internal_vars_de$cell_types == cell_type), myRCTD)
    } else { #multi
      plot_sig_genes_quant_regions(cell_type,myRCTD@internal_vars_de$all_barc,myRCTD@internal_vars_de$my_beta,myRCTD@originalSpatialRNA,
                           res_genes,datadir,myRCTD@internal_vars_de$X2, myRCTD@de_results$gene_fits, which(myRCTD@internal_vars_de$cell_types == cell_type))
    }
  }
}

make_de_plots_spatial <- function(myRCTD, datadir) {
  if(!dir.exists(file.path(datadir,'de_plots')))
    dir.create(file.path(datadir,'de_plots'))
  for(cell_type in myRCTD@internal_vars_de$cell_types) {
    res_genes <- myRCTD@de_results$res_gene_list[[cell_type]]
    #res_genes <- res_genes[1:min(75,dim(res_genes)[1]),]
    plot_sig_genes_two_regions(cell_type,myRCTD@internal_vars_de$all_barc,myRCTD@internal_vars_de$my_beta,myRCTD@originalSpatialRNA,res_genes,datadir,myRCTD@internal_vars_de$X2)
  }
}

make_de_plots_regions <- function(myRCTD, datadir) {
  if(!dir.exists(file.path(datadir,'de_plots')))
    dir.create(file.path(datadir,'de_plots'))
  for(cell_type in myRCTD@internal_vars_de$cell_types) {
    res_genes <- myRCTD@de_results$res_gene_list[[cell_type]]
    #res_genes <- res_genes[1:min(75,dim(res_genes)[1]),]
    plot_sig_genes_regions(cell_type,myRCTD@internal_vars_de$all_barc,myRCTD@internal_vars_de$my_beta,myRCTD@originalSpatialRNA,res_genes,datadir,myRCTD@internal_vars_de$X2)
  }
}

make_de_plots_predictions <- function(myRCTD, datadir, test_mode = 'direct') {
  if(!dir.exists(file.path(datadir,'prediction_plots')))
    dir.create(file.path(datadir,'prediction_plots'))
  for(cell_type in myRCTD@internal_vars_de$cell_types) {
    res_genes <- myRCTD@de_results$res_gene_list[[cell_type]]
    #res_genes <- res_genes[1:min(75,dim(res_genes)[1]),]
    plot_prediction_genes(cell_type,myRCTD@internal_vars_de$all_barc,myRCTD@internal_vars_de$my_beta,myRCTD@originalSpatialRNA,res_genes,test_mode, datadir)
  }
}

write_de_summary <- function(myRCTD, datadir) {
  if(!dir.exists(file.path(datadir,'de_summary')))
    dir.create(file.path(datadir,'de_summary'))
  for(cell_type in myRCTD@internal_vars_de$cell_types) {
    res_genes <- myRCTD@de_results$res_gene_list[[cell_type]]
    write.csv(res_genes,file.path(datadir,'de_summary',paste0(stringr::str_replace(cell_type,'/','_'),'.csv')))
  }
}

plot_prediction_genes <- function(cell_type, all_barc, my_beta, puck, res_genes, test_mode, datadir) {
  if(!dir.exists(file.path(datadir,'prediction_plots')))
    dir.create(file.path(datadir,'prediction_plots'))
  MULT <- 500
  gene_list_sig <- rownames(res_genes)
  if(test_mode == 'multi')
    sing_thresh <- 0.8
  else
    sing_thresh <- 0.999
  barcodes_sing <- names(which(my_beta[all_barc,cell_type] > sing_thresh))
  plots <- list()
  if(length(gene_list_sig) > 0 & length(barcodes_sing) > 0) {
    for(i in 1:length(gene_list_sig)) {#length(gene_list_sig)
      gene = gene_list_sig[i]
      predictions <- predict_DEGLAM(2, gene_fits, gene, myRCTDde@internal_vars_de$X2[barcodes_sing,])[,1]
      plots[[i]] <- plot_puck_continuous(puck, barcodes_sing, predictions*MULT,
                                         ylimit = c(0,quantile(predictions*MULT, 0.95)), title = gene)
    }
    pdf(file.path(datadir,paste0("prediction_plots/de_genes_",stringr::str_replace(cell_type,'/','_'),".pdf")))
    invisible(lapply(plots, print))
    dev.off()
  }
}

plot_sig_genes_two_regions <- function(cell_type, all_barc, my_beta, puck, res_genes, plotdir, X2) {
  gene_list_sig <- rownames(res_genes)
  barcodes_sing <- names(which(my_beta[all_barc,cell_type] > 0.999))
  MULT = 500
  
  barc_plot <- intersect(barcodes_sing,colnames(puck@counts)[puck@nUMI >= 200])
  
  plots <- list()
  if(length(gene_list_sig) > 0) {
    for(i in 1:length(gene_list_sig)) {#length(gene_list_sig)
      gene = gene_list_sig[i]
      Y_plot <- MULT*puck@counts[gene,]/puck@nUMI
      my_title = paste(gene,": log_fc=",round(res_genes[gene,'log_fc'],2),', p=',res_genes[gene,'p_val'])
      my_class <- rep(0,length(barc_plot)); names(my_class) <- barc_plot
      my_class[(X2[barc_plot,2] <= 0.2) & (Y_plot[barc_plot] == 0)] <- 1
      my_class[(X2[barc_plot,2] <= 0.2) & (Y_plot[barc_plot] > 0)] <- 3
      my_class[(X2[barc_plot,2] > 0.2) & (Y_plot[barc_plot] == 0)] <- 2
      my_class[(X2[barc_plot,2] > 0.2) & (Y_plot[barc_plot] > 0)] <- 4
      p3 <- plot_class(puck, barc_plot[order(my_class[barc_plot])], factor(my_class), title = my_title)
      suppressMessages(p3 <- p3 + scale_color_manual(values=c("#CCE2EF","#F6DECC","#0072B2","#D55E00")))
      plots[[i]] <- p3
    }
    pdf(file.path(plotdir,paste0("de_plots/de_genes_",stringr::str_replace(cell_type,'/','_'),".pdf")))
    invisible(lapply(plots, print))
    dev.off()
  }
}

plot_sig_genes_regions <- function(cell_type, all_barc, my_beta, puck, res_genes, plotdir, X2) {
  gene_list_sig <- rownames(res_genes)
  barcodes_sing <- names(which(my_beta[all_barc,cell_type] > 0.8))
  if(length(barcodes_sing) < 10) {
    warning(paste0('plot_sig_genes_regions: Not plotting because less than 10 singlet barcodes found for cell type ',cell_type))
    return()
  }
  MULT = 500
  NR <- dim(X2)[2]
  barc_plot <- intersect(barcodes_sing,colnames(puck@counts)[puck@nUMI >= 200])
  region_ids <- unlist(lapply(barc_plot, 
                              function(barcode) which(as.logical(X2[barcode,]))))
  names(region_ids) <- barc_plot
  plots <- list()
  if(length(gene_list_sig) > 0) {
    for(i in 1:length(gene_list_sig)) {#length(gene_list_sig)
      gene = gene_list_sig[i]
      Y_plot <- MULT*puck@counts[gene,]/puck@nUMI
      region_ids_plot <- region_ids
      region_ids_plot[Y_plot[barc_plot] > 0] <- region_ids[Y_plot[barc_plot] > 0] + NR #
      my_title = paste(gene,": log_fc=",round(res_genes[gene,'log_fc'],2),', p=',res_genes[gene,'p_val'])
      my_class <- region_ids_plot
      p3 <- plot_class(puck, barc_plot[order(my_class[barc_plot])], factor(my_class), title = my_title)
      suppressMessages(p3 <- p3 + scale_color_manual(values=c("#CCE2EF","#F6DECC","#66FFD9","#FFDFFF",
                                                              "#0072B2","#D55E00","#009E73","#CC79A7")))
      plots[[i]] <- p3
    }
    pdf(file.path(plotdir,paste0("de_plots/de_genes_",stringr::str_replace(cell_type,'/','_'),".pdf")))
    invisible(lapply(plots, print))
    dev.off()
  }
}

plot_sig_genes_quant <- function(cell_type, all_barc, my_beta, puck, res_genes, plotdir, X2, gene_fits, cell_type_ind, myRCTD) {
  gene_list_sig <- rownames(res_genes)
  barcodes_sing <- names(which(my_beta[all_barc,cell_type] > 0.999))
  MULT = 500
  NR = 5;
  cur_barc <- intersect(barcodes_sing,colnames(puck@counts)[puck@nUMI >= 200])
  D <- dim(X2)[2] 
  plots <- list()
  if(length(gene_list_sig) > 0) {
    for(i in 1:length(gene_list_sig)) {#length(gene_list_sig)
      gene = gene_list_sig[i]
      my_title = paste(gene,": log_fc=",round(res_genes[gene,'log_fc'],2),', p=',res_genes[gene,'p_val'])
      p3 <- plot_quant_df(myRCTD, gene_fits, myRCTD@internal_vars_de$cell_types, cell_type, gene, NR = NR)
      plots[[i]] <- p3 + ggplot2::ggtitle(my_title) 
    }
    pdf(file.path(plotdir,paste0("de_plots_quant/de_genes_",stringr::str_replace(cell_type,'/','_'),".pdf")))
    invisible(lapply(plots, print))
    dev.off()
  }
}

predict_DEGLAM <- function(cell_type_ind, gene_fits, gene, X2_mat) {
  sigma <- as.numeric(gene_fits$sigma_g[gene])/100
  predictions <- exp(X2_mat %*% gene_fits$all_vals[gene, ,cell_type_ind])
  return(predictions * exp(sigma^2/2))
}

plot_sig_genes_quant_regions <- function(cell_type, all_barc, my_beta, puck, res_genes, plotdir, X2, gene_fits, cell_type_ind) {
  gene_list_sig <- rownames(res_genes)
  barcodes_sing <- names(which(my_beta[all_barc,cell_type] > 0.8))
  if(length(barcodes_sing) < 10) {
    warning(paste0('plot_sig_genes_quant_regions: Not plotting because less than 10 singlet barcodes found for cell type ',cell_type))
    return()
  }
  MULT = 500
  NR = dim(X2)[2];
  cur_barc <- intersect(barcodes_sing,colnames(puck@counts)[puck@nUMI >= 200])
  region_ids <- unlist(lapply(cur_barc, 
                              function(barcode) which(as.logical(X2[barcode,]))))
  names(region_ids) <- cur_barc
  
  plots <- list()
  if(length(gene_list_sig) > 0) {
    for(i in 1:length(gene_list_sig)) {#length(gene_list_sig)
      gene = gene_list_sig[i]
      Y_plot <- MULT*puck@counts[gene,]/puck@nUMI
      my_title = paste(gene,": log_fc=",round(res_genes[gene,'log_fc'],2),', p=',res_genes[gene,'p_val'])
      obs_df <- data.frame(unlist(lapply(1:NR, function(region) mean(Y_plot[cur_barc[region_ids == region]]))),
                           unlist(lapply(1:NR, function(region) length(Y_plot[cur_barc[region_ids == region]]))))
      UMI_list <- unlist(lapply(1:NR, function(region) median(puck@nUMI[cur_barc[region_ids == region]])))
      colnames(obs_df) <- c('mean', 'length');
      obs_df$pred <- exp(unlist(lapply(1:NR, function(region) gene_fits$all_vals[gene,region,cell_type_ind])))
      sigma <- as.numeric(gene_fits$sigma_g[gene])/100
      obs_df$pred <- obs_df$pred*exp(sigma^2/2)*MULT
      ci_width <- sqrt(1/UMI_list*obs_df$pred/MULT)/sqrt(obs_df$length)*MULT*1.96
      obs_df$ub <- obs_df$pred + ci_width
      obs_df$lb <- pmax(obs_df$pred - ci_width,0)
      obs_df$bin <- rownames(obs_df)
      p3 <- ggplot2::ggplot(obs_df,mapping = ggplot2::aes(bin,mean)) + ggplot2::geom_point() + ggplot2::theme_classic() +
        ggplot2::geom_line(ggplot2::aes(bin,pred)) + ggplot2::ggtitle(my_title) + ggplot2::geom_errorbar(aes(ymin = lb, ymax = ub))
      plots[[i]] <- p3
    }
    pdf(file.path(plotdir,paste0("de_plots_quant/de_genes_",stringr::str_replace(cell_type,'/','_'),".pdf")))
    invisible(lapply(plots, print))
    dev.off()
  }
}

get_quant_df <- function(RCTDde, gene_fits, cell_types, cur_cell_types, gene, multi_region = F, prop_thresh = 0.999) {
  my_beta <- RCTDde@internal_vars_de$my_beta
  pred_tot <- numeric(dim(my_beta)[1])
  X2 <- RCTDde@internal_vars_de$X2[rownames(my_beta),2]
  if(multi_region)
    X2 <- apply(myRCTDde@internal_vars_de$X2,1,function(x) which(as.logical(x)))
  nUMI <- RCTDde@spatialRNA@nUMI[rownames(my_beta)]
  Y <- RCTDde@originalSpatialRNA@counts[gene, rownames(my_beta)] / nUMI
  for(cell_type_ind in 1:length(cell_types)) {
    sigma <- as.numeric(gene_fits$sigma_g[gene])/100
    pred_ct <- predict_DEGLAM(cell_type_ind, gene_fits, gene, RCTDde@internal_vars_de$X2[rownames(my_beta),])
    pred_tot <- pred_tot + pred_ct*my_beta[,cell_type_ind]
  }
  a_pred <- pred_tot / exp(sigma^2/2)
  var_pred <- a_pred^2 * exp(sigma^2/2) * (exp(sigma^2/2) - 1)
  var_tot <- var_pred + pred_tot / nUMI
  if(length(cur_cell_types) > 1)
    my_barc <- which(rowSums(my_beta[,cur_cell_types]) > prop_thresh)
  else
    my_barc <- which(my_beta[,cur_cell_types] > prop_thresh)
  if(!multi_region) {
    plot_df <- cbind(pred_tot[my_barc], var_tot[my_barc], my_beta[my_barc,min(which(cell_types %in% cur_cell_types))],X2[my_barc], Y[my_barc])
    colnames(plot_df) <- c('pred','var', 'berg_weight', 'region', 'Y')
  } else {
    plot_df <- cbind(pred_tot[my_barc], var_tot[my_barc], my_beta[my_barc,cur_cell_types],X2[my_barc], Y[my_barc])
    colnames(plot_df) <- c('pred', 'var',cur_cell_types, 'region', 'Y')
  }
  plot_df <- data.frame(plot_df)
}

bin_quant_df <- function(quant_df, NR = 5) {
  results_df <- data.frame(gene = character(), region = integer(), model = integer(),
                           N = integer(), Y = numeric(), pred = numeric(), se = numeric())
  model <- 1
  Threshes <- (0:NR)/NR
  for(cat in 1:NR) {
    cur_barcs <- rownames(quant_df)[quant_df$region >= Threshes[cat] & quant_df$region <= Threshes[cat + 1]]
    n <- length(cur_barcs)
    Y <- mean(quant_df[cur_barcs, 'Y'])
    pred <- mean(quant_df[cur_barcs, 'pred'])
    se <- sqrt(mean(quant_df[cur_barcs,'var'])/n)
    de<-data.frame(gene, cat, model, n, Y, pred, se)
    names(de)<-c('gene','region', 'model','N','Y','pred', 'se')
    results_df <- bind_rows(results_df, de)
  }
  return(results_df)
}

plot_quant_df <- function(myRCTD, gene_fits, cell_types, cell_type, gene, NR = 5) {
  quant_df <- get_quant_df(myRCTD, gene_fits, cell_types, cell_type, gene, multi_region = F, prop_thresh = 0.999)
  final_df <- bin_quant_df(quant_df, NR = NR)
  final_df[,c('Y','se','pred')] <- final_df[,c('Y','se','pred')] * 500 #scale to counts per 500
  final_df$region <- (final_df$region - 1) / (NR - 1)
  p <- ggplot(final_df, aes(x = region, y = pmax(Y, 10^(-8)))) + geom_point() + geom_line(aes(y = pred))+
    geom_errorbar(aes(ymin = pmax(Y - 1.96*se,10^(-8)), ymax = (Y + 1.96*se)), width = 0.05) + 
    theme_classic() + ylab('Average expression') + xlab('Predictive variable') + scale_x_continuous(breaks = c(0,0.5,1), labels = c(0,0.5,1), limits = c(-0.05,1.05))
  return(p)
}
