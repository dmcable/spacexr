#' Makes all CSIDE plots on RCTD object, after running CSIDE
#'
#' @param myRCTD \code{\linkS4class{RCTD}} object containing \code{de_results}, after running CSIDE
#' @param datadir output directory
#' @export
make_all_de_plots <- function(myRCTD, datadir) {
  if(!dir.exists(datadir))
    dir.create(datadir)
  make_de_plots_genes(myRCTD, datadir)
  write_de_summary(myRCTD, datadir)
  if(myRCTD@internal_vars_de$test_mode == 'individual') {
    make_de_plots_quant(myRCTD, datadir)
    make_de_plots_spatial(myRCTD, datadir)
  } else {
    make_de_plots_regions(myRCTD, datadir) #multi
  }
}

#' Makes quantitative CSIDE plots on RCTD object, after running CSIDE
#'
#' @param myRCTD \code{\linkS4class{RCTD}} object containing \code{de_results}, after running CSIDE
#' @param datadir output directory
#' @export
make_de_plots_quant <- function(myRCTD, datadir) {
  if(!dir.exists(file.path(datadir,'de_plots_quant')))
    dir.create(file.path(datadir,'de_plots_quant'))
  for(cell_type in myRCTD@internal_vars_de$cell_types) {
    sig_genes <- myRCTD@de_results$sig_gene_list[[cell_type]]
    if(myRCTD@internal_vars_de$test_mode == 'individual') {
      plot_sig_genes_quant(cell_type,myRCTD@internal_vars_de$barcodes,myRCTD@internal_vars_de$my_beta,myRCTD@originalSpatialRNA,
                             sig_genes,datadir,myRCTD@internal_vars_de$X2, myRCTD@de_results$gene_fits,
                             which(myRCTD@internal_vars_de$cell_types == cell_type), myRCTD)
    } else { #multi
      plot_sig_genes_quant_regions(cell_type,myRCTD@internal_vars_de$barcodes,myRCTD@internal_vars_de$my_beta,myRCTD@originalSpatialRNA,
                           sig_genes,datadir,myRCTD@internal_vars_de$X2, myRCTD@de_results$gene_fits, which(myRCTD@internal_vars_de$cell_types == cell_type))
    }
  }
}

#' Makes spatial gene CSIDE plots (colored continuously) on RCTD object, after running CSIDE
#'
#' @param myRCTD \code{\linkS4class{RCTD}} object containing \code{de_results}, after running CSIDE
#' @param datadir output directory
#' @export
make_de_plots_genes <- function(myRCTD, datadir) {
  for(cell_type in myRCTD@internal_vars_de$cell_types) {
    sig_genes <- myRCTD@de_results$sig_gene_list[[cell_type]]
    plot_sig_genes(cell_type, myRCTD@internal_vars_de$barcodes, myRCTD@internal_vars_de$my_beta,
                 myRCTD@originalSpatialRNA, sig_genes, myRCTD@internal_vars_de$doublet_mode, datadir)
  }
}

#' Makes spatial gene CSIDE plots (colored by two discrete regions) on RCTD replicates object, after running CSIDE
#'
#' Runs on genes that were identified as significant at the population level
#'
#' These plots are colored by two discrete regions based on high or low explanatory variable values.
#' Bold points represent expressed, whereas unbold points represent pixels not expressing the gene.
#'
#' @param RCTD.replicates a \code{\linkS4class{RCTD.replicates}} object after performing population-level DE inference.
#' @param datadir output directory
#' @export
make_de_plots_replicates <- function(RCTD.replicates, datadir) {
  cell_types <- RCTD.replicates@RCTD.reps[[1]]@internal_vars_de$cell_types
  for(cell_type in cell_types) {
    sig_gene_df <- RCTD.replicates@population_sig_gene_df[[cell_type]]
    i <- 0
    for(myRCTD in RCTD.replicates@RCTD.reps) {
      i <- i + 1
      datadir_ct <- file.path(datadir, paste0('rep_',i))
      if(!dir.exists(datadir_ct))
        dir.create(datadir_ct)
      if(!dir.exists(file.path(datadir_ct,'de_plots_two_regions')))
        dir.create(file.path(datadir_ct,'de_plots_two_regions'))
      plot_sig_genes_two_regions(cell_type,myRCTD@internal_vars_de$barcodes,
                                 myRCTD@internal_vars_de$my_beta,myRCTD@originalSpatialRNA,
                                 sig_gene_df,datadir_ct,myRCTD@internal_vars_de$X2,
                                 param_index_best = 2, log_fc_name = 'log_fc_est',
                                 p_val_name = 'p')
    }
  }
}

#' Makes spatial gene CSIDE plots (colored by two discrete regions) on RCTD object, after running CSIDE
#'
#' These plots are colored by two discrete regions based on high or low explanatory variable values.
#' Bold points represent expressed, whereas unbold points represent pixels not expressing the gene.
#'
#' @param myRCTD \code{\linkS4class{RCTD}} object containing \code{de_results}, after running CSIDE
#' @param datadir output directory
#' @export
make_de_plots_spatial <- function(myRCTD, datadir) {
  if(!dir.exists(file.path(datadir,'de_plots_two_regions')))
    dir.create(file.path(datadir,'de_plots_two_regions'))
  for(cell_type in myRCTD@internal_vars_de$cell_types) {
    sig_genes <- myRCTD@de_results$sig_gene_list[[cell_type]]
    #sig_genes <- sig_genes[1:min(75,dim(sig_genes)[1]),]
    plot_sig_genes_two_regions(cell_type,myRCTD@internal_vars_de$barcodes,myRCTD@internal_vars_de$my_beta,myRCTD@originalSpatialRNA,sig_genes,datadir,myRCTD@internal_vars_de$X2)
  }
}

#' Makes spatial gene CSIDE plots (colored by discrete regions) on RCTD object, after running CSIDE
#'
#' @param myRCTD \code{\linkS4class{RCTD}} object containing \code{de_results}, after running CSIDE
#' @param datadir output directory
#' @export
make_de_plots_regions <- function(myRCTD, datadir) {
  if(!dir.exists(file.path(datadir,'de_plots')))
    dir.create(file.path(datadir,'de_plots'))
  for(cell_type in myRCTD@internal_vars_de$cell_types) {
    sig_genes <- myRCTD@de_results$sig_gene_list[[cell_type]]
    #sig_genes <- sig_genes[1:min(75,dim(sig_genes)[1]),]
    plot_sig_genes_regions(cell_type,myRCTD@internal_vars_de$barcodes,myRCTD@internal_vars_de$my_beta,myRCTD@originalSpatialRNA,sig_genes,datadir,myRCTD@internal_vars_de$X2)
  }
}

make_de_plots_predictions <- function(myRCTD, datadir, test_mode = 'individual') {
  if(!dir.exists(file.path(datadir,'prediction_plots')))
    dir.create(file.path(datadir,'prediction_plots'))
  for(cell_type in myRCTD@internal_vars_de$cell_types) {
    sig_genes <- myRCTD@de_results$sig_gene_list[[cell_type]]
    #sig_genes <- sig_genes[1:min(75,dim(sig_genes)[1]),]
    plot_prediction_genes(cell_type,myRCTD@internal_vars_de$barcodes,myRCTD@internal_vars_de$my_beta,
                          myRCTD@originalSpatialRNA,myRCTD@internal_vars_de$X2,sig_genes,myRCTD@internal_vars_de$doublet_mode, datadir, myRCTD@de_results$gene_fits)
  }
}

#' Saves to csv the CSIDE significant gene dataframes after running CSIDE
#'
#' @param myRCTD \code{\linkS4class{RCTD}} object containing \code{de_results}, after running CSIDE
#' @param datadir output directory
#' @export
write_de_summary <- function(myRCTD, datadir) {
  if(!dir.exists(file.path(datadir,'de_summary')))
    dir.create(file.path(datadir,'de_summary'))
  for(cell_type in myRCTD@internal_vars_de$cell_types) {
    sig_genes <- myRCTD@de_results$sig_gene_list[[cell_type]]
    write.csv(sig_genes,file.path(datadir,'de_summary',paste0("sig_",stringr::str_replace(cell_type,'/','_'),'.csv')))
    all_genes <- myRCTD@de_results$all_gene_list[[cell_type]]
    write.csv(all_genes,file.path(datadir,'de_summary',paste0("all_",stringr::str_replace(cell_type,'/','_'),'.csv')))
  }
}

plot_sig_genes <- function(cell_type, barcodes, my_beta, puck, sig_genes, doublet_mode, datadir) {
  if(!dir.exists(file.path(datadir,'de_plots')))
    dir.create(file.path(datadir,'de_plots'))
  gene_list_sig <- rownames(sig_genes)
  if(!doublet_mode)
    sing_thresh <- 0.8
  else
    sing_thresh <- 0.999
  barcodes_sing <- names(which(my_beta[barcodes,cell_type] > sing_thresh))
  plots <- list()
  if(length(gene_list_sig) > 0 & length(barcodes_sing) > 0) {
    for(i in 1:length(gene_list_sig)) {#length(gene_list_sig)
      gene = gene_list_sig[i]
      Y <- puck@counts[gene, barcodes]
      plots[[i]] <- plot_puck_continuous(puck, barcodes_sing, Y, ylimit = c(0,5), title = paste(gene,": log_fc=",round(sig_genes[gene,'log_fc'],2),', p=',sig_genes[gene,'p_val']))
    }
    pdf(file.path(datadir,paste0("de_plots/de_genes_",stringr::str_replace(cell_type,'/','_'),".pdf")))
    invisible(lapply(plots, print))
    dev.off()
  }
}

plot_prediction_genes <- function(cell_type, barcodes, my_beta, puck, X2, sig_genes, doublet_mode, datadir, gene_fits) {
  if(!dir.exists(file.path(datadir,'prediction_plots')))
    dir.create(file.path(datadir,'prediction_plots'))
  MULT <- 500
  cell_type_ind <- which(cell_type == colnames(my_beta))
  gene_list_sig <- rownames(sig_genes)
  if(!doublet_mode)
    sing_thresh <- 0.8
  else
    sing_thresh <- 0.999
  barcodes_sing <- names(which(my_beta[barcodes,cell_type] > sing_thresh))
  plots <- list()
  if(length(gene_list_sig) > 0 & length(barcodes_sing) > 0) {
    for(i in 1:length(gene_list_sig)) {#length(gene_list_sig)
      gene = gene_list_sig[i]
      predictions <- predict_CSIDE(cell_type_ind, gene_fits, gene, X2[barcodes_sing,])[,1]
      plots[[i]] <- plot_puck_continuous(puck, barcodes_sing, predictions*MULT,
                                         ylimit = c(0,quantile(predictions*MULT, 0.95)), title = gene)
    }
    pdf(file.path(datadir,paste0("prediction_plots/de_genes_",stringr::str_replace(cell_type,'/','_'),".pdf")))
    invisible(lapply(plots, print))
    dev.off()
  }
}

#' Makes a spatial plot of CSIDE fitted gene expression
#'
#' Units counts per 500
#'
#' @param myRCTD \code{\linkS4class{RCTD}} object containing \code{de_results}, after running CSIDE
#' @param gene gene to be plotted
#' @param cell_type cell_type to be plotted (only single cell type pixels)
#' @return plot of fitted gene expression
#' @export
plot_prediction_gene <- function(myRCTD, cell_type, gene) {
  barcodes <- myRCTD@internal_vars_de$barcodes
  my_beta <- myRCTD@internal_vars_de$my_beta
  puck <- myRCTD@spatialRNA
  doublet_mode <- myRCTD@internal_vars_de$doublet_mode
  cell_type_ind <- which(colnames(my_beta) == cell_type)
  MULT <- 500
  if(!doublet_mode)
    sing_thresh <- 0.8
  else
    sing_thresh <- 0.999
  barcodes_sing <- names(which(my_beta[barcodes,cell_type] > sing_thresh))
  plots <- list()
  if(length(barcodes_sing) >= 10) {
    predictions <- predict_CSIDE(cell_type_ind, myRCTD@de_results$gene_fits, gene,
                                 myRCTD@internal_vars_de$X2[barcodes_sing,])[,1]
    p <- plot_puck_continuous(puck, barcodes_sing, predictions*MULT,
                                       ylimit = c(0,quantile(predictions*MULT, 0.95)), title = gene)
    return(p)
  } else {
    print('plot_prediction_gene: not plotting because fewer than 10 cell type singlet pixels found')
  }
}

plot_sig_genes_two_regions <- function(cell_type, barcodes, my_beta, puck, sig_genes, plotdir, X2, param_index_best = NULL,
                                       log_fc_name = 'log_fc', p_val_name = 'p_val') {
  gene_list_sig <- rownames(sig_genes)
  barcodes_sing <- names(which(my_beta[barcodes,cell_type] > 0.999))
  MULT = 500
  barc_plot <- intersect(barcodes_sing,colnames(puck@counts)[puck@nUMI >= 200])
  plots <- list()
  if(length(gene_list_sig) > 0) {
    for(i in 1:length(gene_list_sig)) {#length(gene_list_sig)
      gene = gene_list_sig[i]
      if(gene %in% rownames(puck@counts)) {
        Y_plot <- MULT*puck@counts[gene,]/puck@nUMI
        my_title = paste(gene,": log_fc=",round(sig_genes[gene,log_fc_name],2),', p=',sig_genes[gene,p_val_name])
        my_class <- rep(0,length(barc_plot)); names(my_class) <- barc_plot
        if(is.null(param_index_best))
          param_index <- sig_genes[gene,'paramindex_best']
        else
          param_index <- param_index_best
        my_class[(X2[barc_plot,param_index] <= 0.5) & (Y_plot[barc_plot] == 0)] <- 1
        my_class[(X2[barc_plot,param_index] <= 0.5) & (Y_plot[barc_plot] > 0)] <- 3
        my_class[(X2[barc_plot,param_index] > 0.5) & (Y_plot[barc_plot] == 0)] <- 2
        my_class[(X2[barc_plot,param_index] > 0.5) & (Y_plot[barc_plot] > 0)] <- 4
        p3 <- plot_class(puck, barc_plot[order(my_class[barc_plot])], factor(my_class), title = my_title)
        suppressMessages(p3 <- p3 + scale_color_manual(values=c("#CCE2EF","#F6DECC","#0072B2","#D55E00")))
        plots[[i]] <- p3
      }
    }
    pdf(file.path(plotdir,paste0("de_plots_two_regions/de_genes_",stringr::str_replace(cell_type,'/','_'),".pdf")))
    invisible(lapply(plots, print))
    dev.off()
  }
}

#' Makes a spatial plot of gene expression for a particular gene
#' This plot is colored by two discrete regions based on high or low explanatory variable values.
#' Bold points represent expressed, whereas unbold points represent pixels not expressing the gene.
#'
#' @param myRCTD \code{\linkS4class{RCTD}} object containing \code{de_results}, after running CSIDE
#' @param gene gene to be plotted
#' @param cell_type cell_type to be plotted (only single cell type pixels)
#' @param min_UMI (default 200) minimum UMI for pixels that are included
#' @param expr.thresh (default 0) the minimum expression threshold to clear to be considered to be expressed
#' @param exvar_thresh threshold of the explanatory variable in order for points to be sorted into the two regions
#' @return gene expression plot
#' @export
plot_gene_two_regions <- function(myRCTD, gene, cell_type, min_UMI = 200, expr.thresh = 0, exvar_thresh = 0.5) {
  puck <- myRCTD@spatialRNA
  barcodes <- myRCTD@internal_vars_de$barcodes
  my_beta <- myRCTD@internal_vars_de$my_beta
  X2 <- myRCTD@internal_vars_de$X2
  sig_gene_df <- myRCTD@de_results$sig_gene_list[[cell_type]]
  barcodes_sing <- names(which(my_beta[barcodes,cell_type] > 0.999))
  MULT = 500
  barc_plot <- intersect(barcodes_sing,colnames(puck@counts)[puck@nUMI >= min_UMI])
  Y_plot <- MULT*puck@counts[gene,]/puck@nUMI
  my_title = paste(gene,": log_fc=",round(myRCTD@de_results$gene_fits$mean_val[gene,cell_type],2))
  my_class <- rep(0,length(barc_plot)); names(my_class) <- barc_plot
  my_class[(X2[barc_plot,2] <= exvar_thresh) & (Y_plot[barc_plot] == expr.thresh)] <- 1
  my_class[(X2[barc_plot,2] <= exvar_thresh) & (Y_plot[barc_plot] > expr.thresh)] <- 3
  my_class[(X2[barc_plot,2] > exvar_thresh) & (Y_plot[barc_plot] == expr.thresh)] <- 2
  my_class[(X2[barc_plot,2] > exvar_thresh) & (Y_plot[barc_plot] > expr.thresh)] <- 4
  p3 <- plot_class(puck, barc_plot[order(my_class[barc_plot])], factor(my_class), title = my_title)
  suppressMessages(p3 <- p3 + scale_color_manual(values=c("#CCE2EF","#F6DECC","#0072B2","#D55E00")))
  return(p3)
}

plot_sig_genes_regions <- function(cell_type, barcodes, my_beta, puck, sig_genes, plotdir, X2) {
  gene_list_sig <- rownames(sig_genes)
  barcodes_sing <- names(which(my_beta[barcodes,cell_type] > 0.8))
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
      my_title = paste(gene,": log_fc=",round(sig_genes[gene,'log_fc'],2),', p=',sig_genes[gene,'p_val'])
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

plot_sig_genes_quant <- function(cell_type, barcodes, my_beta, puck, sig_genes, plotdir, X2, gene_fits, cell_type_ind, myRCTD) {
  gene_list_sig <- rownames(sig_genes)
  barcodes_sing <- names(which(my_beta[barcodes,cell_type] > 0.999))
  MULT = 500
  NR = 5;
  cur_barc <- intersect(barcodes_sing,colnames(puck@counts)[puck@nUMI >= 200])
  D <- dim(X2)[2]
  plots <- list()
  if(length(gene_list_sig) > 0) {
    for(i in 1:length(gene_list_sig)) {#length(gene_list_sig)
      gene = gene_list_sig[i]
      my_title = paste(gene,": log_fc=",round(sig_genes[gene,'log_fc'],2),', p=',sig_genes[gene,'p_val'])
      p3 <- plot_quant_df(myRCTD, gene_fits, myRCTD@internal_vars_de$cell_types, cell_type, gene, NR = NR,
                          param_position = sig_genes[gene, 'paramindex_best'])
      plots[[i]] <- p3 + ggplot2::ggtitle(my_title)
    }
    pdf(file.path(plotdir,paste0("de_plots_quant/de_genes_",stringr::str_replace(cell_type,'/','_'),".pdf")))
    invisible(lapply(plots, print))
    dev.off()
  }
}

predict_CSIDE <- function(cell_type_ind, gene_fits, gene, X2_mat) {
  sigma <- as.numeric(gene_fits$sigma_g[gene])/100
  predictions <- exp(X2_mat %*% gene_fits$all_vals[gene, ,cell_type_ind])
  return(predictions * exp(sigma^2/2))
}

plot_sig_genes_quant_regions <- function(cell_type, barcodes, my_beta, puck, sig_genes, plotdir, X2, gene_fits, cell_type_ind) {
  gene_list_sig <- rownames(sig_genes)
  barcodes_sing <- names(which(my_beta[barcodes,cell_type] > 0.8))
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
      my_title = paste(gene,": log_fc=",round(sig_genes[gene,'log_fc'],2),', p=',sig_genes[gene,'p_val'])
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

predict_CSIDE_all <- function(RCTDde, gene) {
  cell_types <- RCTDde@internal_vars_de$cell_types
  gene_fits <- RCTDde@de_results$gene_fits
  my_beta <- RCTDde@internal_vars_de$my_beta
  pred_tot <- numeric(dim(my_beta)[1])
  sigma <- as.numeric(gene_fits$sigma_g[gene])/100
  for(cell_type_ind in 1:length(cell_types)) {
    pred_ct <- predict_CSIDE(cell_type_ind, gene_fits, gene, RCTDde@internal_vars_de$X2[rownames(my_beta),])
    pred_tot <- pred_tot + pred_ct*my_beta[,cell_type_ind]
  }
  return(pred_tot)
}

get_quant_df <- function(RCTDde, gene_fits, cell_types, cur_cell_types, gene, multi_region = F, prop_thresh = 0.999, param_position = 2) {
  my_beta <- RCTDde@internal_vars_de$my_beta
  X2 <- RCTDde@internal_vars_de$X2[rownames(my_beta),param_position]
  if(multi_region)
    X2 <- apply(RCTDde@internal_vars_de$X2,1,function(x) which(as.logical(x)))
  nUMI <- RCTDde@spatialRNA@nUMI[rownames(my_beta)]
  Y <- RCTDde@originalSpatialRNA@counts[gene, rownames(my_beta)] / nUMI
  sigma <- as.numeric(gene_fits$sigma_g[gene])/100
  pred_tot <- predict_CSIDE_all(RCTDde, gene)
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
  gene <- 'gene'
  Threshes <- (0:NR)/NR
  for(cat in 1:NR) {
    cur_barcs <- rownames(quant_df)[quant_df$region >= Threshes[cat] & quant_df$region <= Threshes[cat + 1]]
    n <- length(cur_barcs)
    Y <- mean(quant_df[cur_barcs, 'Y'])
    pred <- mean(quant_df[cur_barcs, 'pred'])
    se <- sqrt(mean(quant_df[cur_barcs,'var'])/n)
    de<-data.frame(gene, cat, model, n, Y, pred, se)
    names(de)<-c('gene','region', 'model','N','Y','pred', 'se')
    results_df <- dplyr::bind_rows(results_df, de)
  }
  return(results_df)
}

plot_quant_df <- function(myRCTD, gene_fits, cell_types, cell_type, gene, NR = 5, param_position = 2) {
  quant_df <- get_quant_df(myRCTD, gene_fits, cell_types, cell_type, gene, multi_region = F,
                           prop_thresh = 0.999, param_position = param_position)
  final_df <- bin_quant_df(quant_df, NR = NR)
  final_df[,c('Y','se','pred')] <- final_df[,c('Y','se','pred')] * 500 #scale to counts per 500
  final_df$region <- (final_df$region - 1) / (NR - 1)
  p <- ggplot(final_df, aes(x = region, y = pmax(Y, 10^(-8)))) + geom_point() + geom_line(aes(y = pred))+
    geom_errorbar(aes(ymin = pmax(Y - 1.96*se,10^(-8)), ymax = (Y + 1.96*se)), width = 0.05) +
    theme_classic() + ylab('Average expression') + xlab('Predictive variable') + scale_x_continuous(breaks = c(0,0.5,1), labels = c(0,0.5,1), limits = c(-0.05,1.05))
  return(p)
}

#' Makes a spatial plot of gene expression for a particular gene
#' This plot is colored by several discrete regions based on a categorical design matrix.
#' Bold points represent expressed, whereas unbold points represent pixels not expressing the gene.
#'
#' @param myRCTD \code{\linkS4class{RCTD}} object containing \code{de_results}, after running CSIDE
#' @param gene gene to be plotted
#' @param cell_type cell_type to be plotted (only single cell type pixels)
#' @param pixel_weight_thresh (default 0.8) minimum cell type weight for pixels that are included
#' @param expr_thresh (default 0) the minimum expression threshold to clear to be considered to be expressed
#' @return gene expression plot
#' @export
plot_gene_regions <- function(myRCTD, cell_type, gene, pixel_weight_thresh = 0.8, expr_thresh = 0) {
  barcodes <- myRCTD@internal_vars_de$barcodes
  my_beta <- myRCTD@internal_vars_de$my_beta
  puck <- myRCTD@spatialRNA
  X2 <- myRCTD@internal_vars_de$X2

  barcodes_sing <- names(which(my_beta[barcodes,cell_type] > pixel_weight_thresh))
  if(length(barcodes_sing) < 10) {
    warning(paste0('plot_gene_regions: Not plotting because less than 10 singlet barcodes found for cell type ',cell_type))
    return()
  }
  MULT = 500
  NR <- dim(X2)[2]
  barc_plot <- intersect(barcodes_sing,colnames(puck@counts)[puck@nUMI >= 200])
  region_ids <- unlist(lapply(barc_plot,
                              function(barcode) which(as.logical(X2[barcode,]))))
  names(region_ids) <- barc_plot
  Y_plot <- MULT*puck@counts[gene,]/puck@nUMI

  color_list <- c("#CCE2EF","#F6DECC","#66FFD9","#FFDFFF",
                  "#0072B2","#D55E00","#009E73","#CC79A7")
  if(NR > length(color_list)/2)
    stop(paste0('plot_gene_regions: Not plotting because not implemented for more regions than ', length(color_list)/2))
  region_ids_plot <- region_ids
  region_ids_plot[Y_plot[barc_plot] > expr_thresh] <- region_ids[Y_plot[barc_plot] > expr_thresh] + length(color_list)/2 #
  my_title = paste(gene)
  my_class <- region_ids_plot
  p3 <- plot_class(puck, barc_plot[order(my_class[barc_plot])], factor(my_class), title = my_title)
  color_list <- color_list[c(1:NR, (1:NR) + length(color_list)/2)]
  suppressMessages(p3 <- p3 + scale_color_manual(values=color_list))
  return(p3)
}

#' Makes a spatial plot of continuous gene expression for a particular gene
#'
#' Units counts per 500
#'
#' @param myRCTD \code{\linkS4class{RCTD}} object containing \code{de_results}, after running CSIDE
#' @param gene gene to be plotted
#' @param cell_type cell_type to be plotted (only single cell type pixels)
#' @param ymax (default 10) maximum expression (in counts per 500) for color scale
#' @return gene expression plot
#' @export
plot_gene_raw <- function(myRCTD, gene, cell_type, ymax = 10) {
  my_beta <- myRCTD@internal_vars_de$my_beta
  barcodes <- myRCTD@internal_vars_de$barcodes
  sing_thresh <- 0.8; MULT <- 500
  barcodes_sing <- names(which(my_beta[barcodes,cell_type] > sing_thresh))
  Y <- myRCTD@originalSpatialRNA@counts[gene, barcodes_sing] /
    myRCTD@originalSpatialRNA@nUMI[barcodes_sing] * MULT
  plot_puck_continuous(myRCTD@originalSpatialRNA, barcodes_sing, Y, ylimit = c(0,ymax), title = gene)
}
