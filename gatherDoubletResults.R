library(RCTD)
library(Matrix)
library(dplyr)
library(ggplot2)
require(reshape2)
source('Plotting/figure_utils.R')
iv <- init_RCTD(gene_list_reg = F, get_proportions = F)
puck = iv$puck
iv <- init_RCTD(load_info_renorm = T) #initial variables
resultsdir <- paste0(iv$slideseqdir,"/results")
cell_type_names = iv$cell_type_info[[2]]
results <- list()
for (fold_index in 1:iv$n_puck_folds) {
  results <- append(results, readRDS(paste0(iv$slideseqdir,"/SplitPuckResults/results",fold_index,".RDS")))
}
barcodes <- colnames(puck@counts)
N <- length(results)
weights = Matrix(0, nrow = N, ncol = iv$cell_type_info[[3]])
weights_doublet = Matrix(0, nrow = N, ncol = 2)
rownames(weights) = barcodes; rownames(weights_doublet) = barcodes
colnames(weights) = iv$cell_type_info[[2]]; colnames(weights_doublet) = c('first_type', 'second_type')
empty_cell_types = factor(character(N),levels = cell_type_names)
spot_levels <- c("reject", "singlet", "doublet_certain", "doublet_certain_class", "doublet_uncertain")
doublet_mode <- T
if(doublet_mode) {
  spot_levels <- c("reject", "singlet", "doublet_certain", "doublet_uncertain")
  results_df <- data.frame(spot_class = factor(character(N),levels=spot_levels),
                           first_type = empty_cell_types, second_type = empty_cell_types,
                           first_class = logical(N), second_class = logical(N),
                           min_score = numeric(N), singlet_score = numeric(N),
                           conv_all = logical(N), conv_doublet = logical(N))
  for(i in 1:N) {
    if(i %% 100 == 0)
      print(i)
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
  }
} else {
  results_df <- data.frame(spot_class = factor(character(N),levels=spot_levels),
                           first_type = empty_cell_types, second_type = empty_cell_types,
                           other_type = empty_cell_types, min_score = numeric(N), second_score = numeric(N),
                           conv_all = logical(N), conv_doublet = logical(N))
  for(i in 1:N) {
    weights_doublet[i,] = results[[i]]$doublet_weights
    weights[i,] = results[[i]]$all_weights
    results_df[i, "spot_class"] = results[[i]]$spot_class
    results_df[i, "first_type"] = results[[i]]$first_type
    results_df[i, "second_type"] = results[[i]]$second_type
    results_df[i, "other_type"] = results[[i]]$other_type
    results_df[i, "min_score"] = results[[i]]$min_score
    results_df[i, "second_score"] = results[[i]]$second_score
    results_df[i, "conv_all"] = results[[i]]$conv_all
    results_df[i, "conv_doublet"] = results[[i]]$conv_doublet
  }
}
rownames(results_df) = barcodes
marker_data = get_marker_data(iv$cell_type_info[[2]], NULL, iv$cell_type_info[[1]], iv$gene_list, marker_provided = TRUE)
norm_weights = sweep(weights, 1, rowSums(weights), '/')

#slideseq here
marker_scores_df <- get_marker_scores(marker_data, puck, iv$cell_type_info)
norm_marker_scores <- sweep(marker_scores_df, 1, rowSums(marker_scores_df), '/')

#make the plots
plot_meta_genes(iv$cell_type_info,marker_scores_df, resultsdir)
plot_norm_meta_genes(iv$cell_type_info,norm_marker_scores, resultsdir)
plot_weights(iv$cell_type_info, puck, resultsdir, norm_weights)
plot_weights_unthreshold(iv$cell_type_info, puck, resultsdir, norm_weights)
plot_confidence_rate(puck, resultsdir, norm_weights)
plot_weight_distribution(iv$cell_type_info, puck, resultsdir, norm_weights)
get_corr(iv$cell_type_info, norm_marker_scores, norm_weights)
plot_weights_doublet(iv$cell_type_info, puck, resultsdir, weights_doublet)
plot_cond_occur(cell_type_info, resultsdir, norm_weights)
plot_occur_unthreshold(cell_type_info, resultsdir, norm_weights)

#doublets :)
doublets <- results_df[results_df$spot_class == "doublet_certain",]
plot_doublets(doublets, resultsdir, iv$cell_type_info)
plot_doublets_type(doublets, resultsdir, iv$cell_type_info)
doub_occur <- table(doublets$second_type, doublets$first_type)
plot_heat_map(doub_occur/sum(doub_occur)*20, normalize = F)
plot_coloc(results_df, puck, resultsdir)

#decomposing doublets
get_decomposed_data <- function(results_df, iv, puck, weights_doublet) {
  doublets <- results_df[results_df$spot_class == "doublet_certain",]
  first_DGE <- Matrix(0, nrow = dim(doublets)[1], ncol = length(iv$gene_list))
  second_DGE <- Matrix(0, nrow = dim(doublets)[1], ncol = length(iv$gene_list))
  rownames(first_DGE) = rownames(doublets); rownames(second_DGE) = rownames(doublets)
  colnames(first_DGE) = iv$gene_list; colnames(second_DGE) = iv$gene_list
  for(ind in 1:dim(doublets)[1]) {
    print(ind)
    barcode = rownames(doublets)[ind]
    doub_res <- decompose_doublet(puck@counts[iv$gene_list,barcode], weights_doublet[barcode,], iv$gene_list, cell_type_info, results_df[barcode,"first_type"],results_df[barcode,"second_type"])
    first_DGE[barcode,] <- doub_res$expect_1; second_DGE[barcode,] <- doub_res$expect_2
  }
  singlet_id <- results_df$spot_class == "singlet"
  norm1 <- sweep(first_DGE, 1, weights_doublet[rownames(doublets),"first_type"] * puck@nUMI[rownames(first_DGE)], '/')
  norm2 <- sweep(second_DGE, 1, weights_doublet[rownames(doublets),"second_type"] * puck@nUMI[rownames(second_DGE)], '/')
  norm_sing <- sweep(t(puck@counts[iv$gene_list, singlet_id]),1,puck@nUMI[singlet_id],'/')
  all_DGE <- rbind(norm1, norm2, norm_sing)
  cell_type_labels <- unlist(list(doublets$first_type, doublets$second_type, results_df[singlet_id, "first_type"]))
  coords <- rbind(puck@coords[rownames(doublets),c('x','y')], puck@coords[rownames(doublets),c('x','y')], puck@coords[singlet_id,c('x','y')])
  nUMI <- c(weights_doublet[rownames(doublets),"first_type"] *puck@nUMI[rownames(first_DGE)], weights_doublet[rownames(doublets),"second_type"]*puck@nUMI[rownames(second_DGE)], puck@nUMI[singlet_id])
  rownames(coords) = 1:dim(coords)[1]; names(nUMI) = 1:dim(coords)[1]
  rownames(all_DGE) = 1:dim(coords)[1]
  puck_d <- Slideseq(coords, t(all_DGE), nUMI)
  puck_d@cell_labels <- cell_type_labels
  names(puck_d@cell_labels) = 1:dim(coords)[1]
  puck_d@cell_type_names <- cell_type_info[[2]]
  return(puck_d)
}

puck_d <- get_decomposed_data(results_df, iv, puck, weights_doublet)

plot_puck_gene(puck_d, "Sparcl1", cell_type = "Purkinje", min_val = 0, max_val = .03)
plot_puck_gene(puck_d, "Sparcl1", cell_type = "Bergmann", min_val = 0, max_val = .03)
mean(puck_d@counts["Sparcl1", puck_d@cell_labels == "Purkinje"])
mean(puck_d@counts["Sparcl1", puck_d@cell_labels == "Bergmann"])
my_barc <- colnames(puck_d@counts)[puck_d@cell_labels == "Purkinje"]
mean(puck_d@counts["Sparcl1",my_barc])
my_barc <- colnames(puck_d@counts)[puck_d@cell_labels == "Bergmann"]
mean(puck_d@counts["Sparcl1",my_barc])

plot_puck_continuous(puck_d, colnames(puck_d@counts)[puck_d@cell_labels == "Purkinje"], puck_d@counts["Sparcl1",], ylimit = c(0,0.03))
plot_puck_continuous(puck_d, puck_d@cell_labels == "Bergmann", puck_d@counts["Sparcl1",], ylimit = c(0,0.03))
abs(val2[puck_d@cell_labels == "Purkinje"] - val2[puck_d@cell_labels == "Purkinje"])

my_barc <- barcodes[results_df$first_type == "Bergmann" & (results_df$spot_class != "reject")]
norm_gene <- puck@counts["Sparcl1",] / puck@nUMI
plot_puck_continuous(puck, my_barc, norm_gene, ylimit = c(0,0.03))
mean(pmin(norm_gene[my_barc],0.03))
mean(norm_gene[my_barc])
#otherwise you use all of it




#figure 1 and figure 2 here
metadir <- file.path(iv$slideseqdir,"MetaData")
meta_data <- readRDS(file.path(metadir,"meta_data.RDS"))
meta_df <- meta_data$meta_df
UMI_tot <- meta_data$UMI_tot; UMI_list <- meta_data$UMI_list
class_df <- get_class_df(cell_type_names)
DropViz <- F
if(DropViz) {
  common_cell_types = c("Astrocytes", "Bergmann", "Endothelial", "Fibroblast", "Golgi", "Granule", "MLI1", "MLI2", "Oligodendrocytes", "Polydendrocytes", "Purkinje", "UBCs")
} else {
  common_cell_types <- iv$cell_type_info[[2]]
}
resultsdir = file.path(iv$slideseqdir, "results")
square_results <- plot_doublet_identification(meta_data, common_cell_types, resultsdir, meta_df, results_df, class_df, UMI_list)
plot_heat_map(as.matrix(square_results), file_loc = file.path(resultsdir,'doublet_avg_accuracy.png'), save.file = T, normalize = F)
write.csv(as.matrix(square_results), file.path(resultsdir,'doublet_avg_accuracy.csv'))
square_results <- plot_doublet_identification_certain(meta_data, common_cell_types, resultsdir, meta_df, results_df, class_df, UMI_list)
plot_heat_map(as.matrix(square_results), file_loc = file.path(resultsdir,'doublet_avg_accuracy_certain.png'), save.file = T, normalize = F)
write.csv(as.matrix(square_results), file.path(resultsdir,'doublet_avg_accuracy_certain.csv'))
#next make the confusion matrix
true_types = unlist(list(meta_df[meta_df$first_UMI == 0,"second_type"], meta_df[meta_df$first_UMI == UMI_tot,"first_type"]))
pred_types = unlist(list(results_df[meta_df$first_UMI == 0, "first_type"], results_df[meta_df$first_UMI == UMI_tot, "first_type"]))
conf_mat <- caret::confusionMatrix(pred_types,factor(true_types,levels = iv$cell_type_info[[2]]))
plot_heat_map(conf_mat$table, file_loc = file.path(resultsdir,'confusion_matrix.png'), save.file = T)
write.csv(conf_mat$table, file.path(resultsdir,'confusion_matrix.csv'))
if(F) {
  #plot doublet weights
  singlets = meta_df$first_UMI == 0 | meta_df$first_UMI == UMI_tot
  hist(weights_doublet[singlets,"first_type"],breaks = 30)
  p1 <- hist(weights_doublet[singlets,"first_type"], breaks = 20)
  p2 <- hist(weights_doublet[!singlets,"first_type"], breaks = 20)
  plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1))  # first histogram
  plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T)

  #plot all weights
  p1 <- hist(apply(weights[singlets,],1,max), breaks = 20)
  p2 <- hist(apply(weights[!singlets,],1,max), breaks = 20)
  plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1))  # first histogram
  plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T)
}

#plot different conditions
N_UMI_cond <- (length(UMI_list) + 1)/2
spot_class_df <- as.data.frame(Matrix(0,nrow = N_UMI_cond, ncol = length(spot_levels)))
rownames(spot_class_df) <- UMI_list[1:N_UMI_cond]; colnames(spot_class_df) <- spot_levels
for(nUMI in UMI_list[1:N_UMI_cond]) {
  vec = table(results_df[barcodes[meta_df$first_UMI == nUMI],]$spot_class)
  spot_class_df[as.character(nUMI), names(vec)] = as.numeric(vec)
}
spot_class_df <- sweep(spot_class_df, 1, rowSums(spot_class_df),'/')
spot_class_df[,"nUMI"] <- UMI_list[1:N_UMI_cond]
df <- melt(spot_class_df,  id.vars = 'nUMI', variable.name = 'series')
pdf(file.path(resultsdir, "doublet_categories.pdf"))
ggplot(df, aes(nUMI,value)) + geom_line(aes(colour = series)) + ggplot2::ylim(c(0,1))
dev.off()

#plot doublet proportions
N_UMI_cond <- (length(UMI_list) + 1)/2
spot_class_df <- as.data.frame(Matrix(0,nrow = N_UMI_cond, ncol = 2))
rownames(spot_class_df) <- UMI_list[1:N_UMI_cond]; colnames(spot_class_df) <- c('singlets','doublets')
for(nUMI in UMI_list[1:N_UMI_cond]) {
  vec = table(results_df[barcodes[meta_df$first_UMI == nUMI],]$spot_class)
  spot_class_df[as.character(nUMI), 'singlets'] = vec["singlet"]
  spot_class_df[as.character(nUMI), 'doublets'] = sum(vec) - vec["singlet"] - vec["reject"]
}
spot_class_df <- sweep(spot_class_df, 1, rowSums(spot_class_df),'/')
spot_class_df[,"nUMI"] <- UMI_list[1:N_UMI_cond]
df <- melt(spot_class_df,  id.vars = 'nUMI', variable.name = 'series')
pdf(file.path(resultsdir, "doublet_detection.pdf"))
ggplot(df, aes(nUMI,value)) + geom_line(aes(colour = series)) + ggplot2::ylim(c(0,1))
dev.off()

