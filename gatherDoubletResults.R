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

marker_data = get_marker_data(cell_type_info[[2]], NULL, cell_type_info[[1]], iv$gene_list, marker_provided = TRUE)
norm_weights = sweep(weights, 1, rowSums(weights), '/')
marker_scores_df <- get_marker_scores(cell_type_info, marker_data, puck)
norm_marker_scores <- sweep(marker_scores_df, 1, rowSums(marker_scores_df), '/')

#make the plots
plot_meta_genes(cell_type_info,marker_scores_df, resultsdir)
plot_norm_meta_genes(cell_type_info,norm_marker_scores, resultsdir)
plot_weights(iv$cell_type_info, puck, resultsdir, norm_weights)
plot_weights_unthreshold(iv$cell_type_info, puck, resultsdir, norm_weights)
plot_confidence_rate(puck, resultsdir, norm_weights)

#only for ground truth
metadir <- file.path(iv$slideseqdir,"MetaData")
meta_data <- readRDS(file.path(metadir,"meta_data.RDS"))
meta_df <- meta_data$meta_df
UMI_tot <- meta_data$UMI_tot; UMI_list <- meta_data$UMI_list
UMI_list <- c(0,250,500,750,1000) #for now
class_df <- get_class_df(cell_type_names)
common_cell_types = c("Astrocytes", "Bergmann", "Endothelial", "Fibroblast", "Golgi", "Granule", "MLI1", "MLI2", "Oligodendrocytes", "Polydendrocytes", "Purkinje", "UBCs")
resultsdir = file.path(iv$slideseqdir, "results")
square_results <- plot_doublet_identification(meta_data, common_cell_types, resultsdir, meta_df, results_df, class_df, UMI_list)
square_results <- plot_doublet_identification_certain(meta_data, common_cell_types, resultsdir, meta_df, results_df, class_df, UMI_list)
plot_heat_map(as.matrix(square_results), file_loc = file.path(resultsdir,'doublet_avg_accuracy.png'), save.file = T, normalize = F)
write.csv(as.matrix(square_results), file.path(resultsdir,'doublet_avg_accuracy.csv'))
#next make the confusion matrix
true_types = unlist(list(meta_df[meta_df$first_UMI == 0,"second_type"], meta_df[meta_df$first_UMI == UMI_tot,"first_type"]))
pred_types = unlist(list(results_df[meta_df$first_UMI == 0, "first_type"], results_df[meta_df$first_UMI == UMI_tot, "first_type"]))
conf_mat <- caret::confusionMatrix(pred_types,factor(true_types,levels = iv$cell_type_info[[2]]))
plot_heat_map(conf_mat$table, file_loc = file.path(resultsdir,'confusion_matrix.png'), save.file = T)
write.csv(conf_mat$table, file.path(resultsdir,'confusion_matrix.csv'))
#plot doublet weights
singlets = meta_df$first_UMI == 0 | meta_df$first_UMI == UMI_tot
hist(weights_doublet[singlets,"first_type"],breaks = 30)
p1 <- hist(weights_doublet[singlets,"first_type"], breaks = 20)
p2 <- hist(weights_doublet[!singlets,"first_type"], breaks = 20)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T)

#plot all weights
do.call(pmax, weights[singlets,])
apply(weights[singlets,],1,max)
p1 <- hist(apply(weights[singlets,],1,max), breaks = 20)
p2 <- hist(apply(weights[!singlets,],1,max), breaks = 20)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T)

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
ggplot(df, aes(nUMI,value)) + geom_line(aes(colour = series)) + ggplot2::ylim(c(0,1))


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
ggplot(df, aes(nUMI,value)) + geom_line(aes(colour = series)) + ggplot2::ylim(c(0,1))


