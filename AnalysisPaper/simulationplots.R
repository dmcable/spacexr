#making plots for learning of simulated scRNA doublets
metadir <- file.path(iv$slideseqdir,"MetaData")
meta_data <- readRDS(file.path(metadir,"meta_data.RDS"))
meta_df <- meta_data$meta_df
UMI_tot <- meta_data$UMI_tot; UMI_list <- meta_data$UMI_list
class_df <- get_class_df(cell_type_names, use_classes = T)

resultsdir = file.path(iv$slideseqdir, "results")
square_results <- plot_doublet_identification(meta_data, common_cell_types, resultsdir, meta_df, results_df, class_df, UMI_list)
plot_heat_map(as.matrix(square_results), file_loc = file.path(resultsdir,'doublet_avg_accuracy.png'), save.file = T, normalize = F)
write.csv(as.matrix(square_results), file.path(resultsdir,'doublet_avg_accuracy.csv'))
square_results <- plot_doublet_identification_certain(meta_data, common_cell_types, resultsdir, meta_df, results_df, class_df, UMI_list)
plot_heat_map(as.matrix(square_results), file_loc = file.path(resultsdir,'doublet_avg_accuracy_certain.png'), save.file = T, normalize = F)
write.csv(as.matrix(square_results), file.path(resultsdir,'doublet_avg_accuracy_certain.csv'))
square_results <- as.matrix(square_results)
#next make the confusion matrix
true_types = unlist(list(meta_df[meta_df$first_UMI == 0,"second_type"], meta_df[meta_df$first_UMI == UMI_tot,"first_type"]))
pred_types = unlist(list(results_df[meta_df$first_UMI == 0, "first_type"], results_df[meta_df$first_UMI == UMI_tot, "first_type"]))
conf_mat <- caret::confusionMatrix(pred_types,factor(true_types,levels = iv$cell_type_info[[2]]))
plot_heat_map(conf_mat$table, file_loc = file.path(resultsdir,'confusion_matrix.png'), save.file = T)
write.csv(conf_mat$table, file.path(resultsdir,'confusion_matrix.csv'))

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
df$std <- sqrt(df$value*(1 - df$value)/(table(meta_df$first_UMI)*2)[1])
df$std[df$nUMI == 500] <- df$std[df$nUMI == 500]*sqrt(2)
pdf(file.path(resultsdir, "doublet_detection.pdf"))
ggplot(df, aes(nUMI,value)) + geom_line(aes(colour = series)) + ggplot2::ylim(c(0,1))
dev.off()
doublet_df <- df
save(iv, meta_data, common_cell_types, meta_df, results_df, class_df, UMI_list, square_results, doublet_df,file = "Plotting/Results/doublet_ref.RData")
