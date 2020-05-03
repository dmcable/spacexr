#used to find genes within astrocytes dependent on environment
#platform effect analysis
true_proportions <- proportions * 0
true_proportions[common_cell_types] = 1
true_proportions = true_proportions / sum(true_proportions)
gene_means <- as.matrix(cell_type_info_unnorm[[1]][iv$gene_list,]) %*% true_proportions
bulk_vec = rowSums(puck@counts);
total_UMI <- sum(puck@nUMI)
true_platform_effect = log(bulk_vec[iv$gene_list] / total_UMI,2) -log(gene_means,2)
pdf(file.path(resultsdir,'platform_histogram.pdf'))
hist(true_platform_effect, breaks =30,xlab = "log2(Platform Effect)", main = "Measured Platform effects between dropviz and 10x")
dev.off()
#gene_means_pred <- as.matrix(iv$cell_type_info[[1]][iv$gene_list,]) %*% true_proportions
#res_platform_effect = log(bulk_vec[iv$gene_list] / total_UMI,2) -log(gene_means_pred,2)
#hist(res_platform_effect, breaks =30,xlab = "log2(Platform Effect)", main = "Measured Platform effects between dropviz and 10x")
gene_means_unnorm_pred <- as.matrix(cell_type_info_unnorm[[1]][iv$gene_list,]) %*% proportions / sum(proportions)
pred_platform_effect = log(bulk_vec[iv$gene_list] / total_UMI,2) -log(gene_means_unnorm_pred,2)
hist(pred_platform_effect, breaks =30,xlab = "log2(Platform Effect)", main = "Measured Platform effects between dropviz and 10x")
R2 = cor(true_platform_effect, pred_platform_effect)^2
df <- data.frame(estimated_platform_effect = pred_platform_effect,true_platform_effect=true_platform_effect)
pdf(file.path(resultsdir,'platform_estimation.pdf'))
ggplot(df,aes(x=estimated_platform_effect,y=true_platform_effect)) + geom_point(alpha = 0.2) + geom_line(aes(x=estimated_platform_effect,y=estimated_platform_effect)) + ggplot2::ggtitle(paste('R2=',R2))
dev.off()
saveRDS(df, file.path(resultsdir,'platform_estimation_df.RDS'))
sig_pe <- c(sd(pred_platform_effect)*log(2), sd(true_platform_effect)*log(2))
names(sig_pe) <- c('sigma_pe_est','sigma_pe_true')
write.csv(sig_pe, file.path(resultsdir,'sig_pe.csv'))
