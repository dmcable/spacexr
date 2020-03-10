library(RCTD)
library(Matrix)
library(dplyr)
library(ggplot2)
require(reshape2)
source('Plotting/figure_utils.R')
iv <- init_RCTD(load_info_renorm = T) #initial variables
cell_type_names = iv$cell_type_info[[2]]
results <- list()
for (fold_index in 1:iv$n_puck_folds) {
  results <- append(results, readRDS(paste0(iv$slideseqdir,"/DecomposeResults/results",fold_index,".RDS")))
}
metadir <- file.path(iv$slideseqdir,"MetaData")
meta_data <- readRDS(file.path(metadir,"meta_data.RDS"))
meta_df <- meta_data$meta_df
barcodes <- rownames(meta_df)
first_beads <- readRDS(file.path(metadir,"firstbeads.RDS"))
second_beads <- readRDS(file.path(metadir,"secbeads.RDS"))
N <- length(results)
expect1 = Matrix(0, nrow = N, ncol = length(iv$gene_list))
expect2 = Matrix(0, nrow = N, ncol = length(iv$gene_list))
var = Matrix(0, nrow = N, ncol = length(iv$gene_list))
weights_doublet = Matrix(0, nrow = N, ncol = 2)
rownames(expect1) = barcodes; rownames(weights_doublet) = barcodes
rownames(expect2) = barcodes; rownames(var) = barcodes
colnames(expect1) = iv$gene_list; colnames(expect2) = iv$gene_list
colnames(var) = iv$gene_list; colnames(weights_doublet) = c('first_type', 'second_type')
for(i in 1:N) {
  weights_doublet[i,] = results[[i]]$weights
  expect1[i,] = results[[i]]$decompose_results$expect_1
  expect2[i,] = results[[i]]$decompose_results$expect_2
  var[i,] = results[[i]]$decompose_results$variance
}

#weight recovery plot
type1 = "Astrocytes"; type2 = "Bergmann"
cur_df <- meta_df[meta_df$first_type == type1 & meta_df$second_type == type2,]
plot_df <- aggregate(weights_doublet[rownames(cur_df),1], list(cur_df$first_UMI/meta_data$UMI_tot), mean)
colnames(plot_df) = c('proportion','type1_avg')
plot_df[,"st_dev"] = aggregate(weights_doublet[rownames(cur_df),1], list(cur_df$first_UMI/meta_data$UMI_tot), sd)$x*1.96
p<- ggplot2::ggplot(plot_df, ggplot2::aes(x=proportion, y=type1_avg, colour = "type1_avg")) +
  ggplot2::geom_line() +
  ggplot2::geom_point()+
  ggplot2::geom_line(ggplot2::aes(y=proportion,colour = "proportion")) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin=type1_avg-st_dev, ymax=type1_avg+st_dev), width=.05,
                         position=ggplot2::position_dodge(0.05))

#RMSE
RMSE <- sqrt(mean((weights_doublet[rownames(cur_df),1] - cur_df$first_UMI/meta_data$UMI_tot)^2))

#
v = as.vector(var[rownames(cur_df),])
x = as.vector(expect1[rownames(cur_df),])
y = as.vector(t(first_beads[iv$gene_list,rownames(cur_df)]))
z = x + as.vector(expect2[rownames(cur_df),]) #total
cor(x[z>0]/z[z>0],y[z>0]/z[z>0])


pred_por = x[z>0]/z[z>0]
true_por = y[z>0]/z[z>0]
my_var = v[z > 0]/z[z > 0]^2
#true vs expected mean squared error.
num_bins = 10
mean_res = numeric(num_bins)
mean_pred = numeric(num_bins)
exp_mse = numeric(num_bins)
true_mse = numeric(num_bins)
null_mse = numeric(num_bins)
for(ind in 1:num_bins) {
  cur_ind = pred_por > (ind-1)/(num_bins) & pred_por < (ind)/(num_bins)
  mean_res[ind] = mean(true_por[cur_ind]); mean_pred[ind] = mean(pred_por[cur_ind])
  true_mse[ind] = mean((true_por[cur_ind] - pred_por[cur_ind])^2)
  exp_mse[ind] = mean(my_var[cur_ind])
  null_mse[ind] = mean((true_por[cur_ind] - 0.5)^2)
}
err_df <- data.frame(mean_pred, exp_mse, true_mse)
err_df <- data.frame(mean_pred, sqrt(exp_mse), sqrt(true_mse), sqrt(null_mse))
colnames(err_df) = c('mean_pred', 'exp_mse', 'true_mse', 'null_mse')
p<- ggplot2::ggplot(err_df, ggplot2::aes(x=mean_pred, y=exp_mse, colour = "exp_mse")) +
  ggplot2::geom_line() +
  ggplot2::geom_point()+
  ggplot2::geom_line(ggplot2::aes(y=true_mse,colour = "true_mse")) +
  ggplot2::geom_line(ggplot2::aes(y=null_mse,colour = "null_mse")) + ggplot2::ylim(c(0,1))
p
plot(mean_pred, mean_res)
#overall predictions distribution
hist(pred_por)
overall_mse <- mean((true_por - pred_por)^2)
null_mse <- mean((true_por - 0.5)^2)
samp_mse <- mean(my_var)
R2 <- (null_mse - overall_mse)/(null_mse - samp_mse)
#DE genes
beads = expect1 + expect2
epsilon = 1e-9
de_score <- log(iv$cell_type_info[[1]][,type1] + epsilon) - log(iv$cell_type_info[[1]][,type2] + epsilon)
names(de_score) <- iv$gene_list
de_score <- de_score[colSums(beads[rownames(cur_df_mid), ] > 0) >= 15] # highly expressed
N_genes <- 10
de_genes <- c(head(de_score[order(de_score)],N_genes), tail(de_score[order(de_score)],N_genes))
de_gene_names <- names(de_genes)
cur_df_mid <- meta_df[meta_df$first_type == type1 & meta_df$second_type == type2 & meta_df$first_UMI == meta_data$UMI_tot/2,]
e_por <- colMeans(expect1[rownames(cur_df_mid), de_gene_names] / (beads[rownames(cur_df_mid), de_gene_names]), na.rm=T)
t_por <- colMeans((t(first_beads[de_gene_names, rownames(cur_df_mid)]) / (beads[rownames(cur_df_mid), de_gene_names])), na.rm=T)
plot(e_por, t_por)
plot(log(e_por), log(t_por))
