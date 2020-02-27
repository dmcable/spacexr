library(RCTD)
library(Matrix)
source('Plotting/figure_utils.R')
config <- config::get()
refdir <- file.path("Data/Reference",config$reffolder)
dir.create(file.path(refdir, "results"), showWarnings = FALSE) #folder to save results
resultsdir = file.path(refdir, "results")
reference <- readRDS(paste(refdir,config$reffile,sep="/"))
cell_type_info <- get_cell_type_info(reference@assays$RNA@counts, reference@meta.data$liger_ident_coarse, reference@meta.data$nUMI)
cell_type_means = cell_type_info[[1]]; cell_type_names = cell_type_info[[2]]
n_cell_types = cell_type_info[[3]]
slideseqdir <- file.path("Data/Slideseq",config$slideseqfolder)
resultsdir = file.path(slideseqdir,"results")
if(!dir.exists(resultsdir))
  dir.create(resultsdir)
bulkdir <- paste(slideseqdir,"results/Bulk",sep="/")
if(!dir.exists(bulkdir))
  dir.create(bulkdir)
print(paste("prepareBulkData: number of cells in reference:", dim(reference@assays$RNA@counts)[2]))
print(paste("prepareBulkData: number of genes in reference:", dim(reference@assays$RNA@counts)[1]))
print("prepareBulkData: amount of each cell type in the reference:")
cell_counts = table(reference@meta.data$liger_ident_coarse)
print(cell_counts)
puck <- seurat.to.slideseq(reference, cell_type_info)
test_reference <- readRDS("Data/Reference/DropVizCerAnnotated/scRefSubsampled1000.RDS")
puck_test <- seurat.to.slideseq(test_reference, cell_type_info)
gene_list = get_de_genes(cell_type_means, puck_test, fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg)
print(paste("callBeads: number of genes used for regression:", length(gene_list)))
#UMI_min = 500
UMI_min = 100
puck = restrict_counts(puck, gene_list, UMI_thresh = UMI_min, UMI_max = 10000000)
puck = restrict_puck(puck, colnames(puck@counts))
table(puck@cell_labels)

#Figure 2A: OLS Prediction works on Training data, but not cross-reference
test_results = process_data(puck, gene_list, cell_type_info, proportions = NULL, trust_model = F, constrain = T, OLS = T)
plot_heat_map(test_results[[1]]$table)

test_results = process_data(puck_test, gene_list, cell_type_info, proportions = NULL, trust_model = F, constrain = T, OLS = T)
plot_heat_map(test_results[[1]]$table)

# note: for the next step in evaluating platform effect estimation, we must have a 1 to 1 mapping between cell types in both references,
# so we restrict the classes to exclude the dropviz clusters that are known to correspond to two clusters in 10x: Macrophages_Microglia and Ependymal_Choroid
common_cell_types = c("Astrocytes", "Bergmann", "Endothelial", "Fibroblast", "Golgi", "Granule", "MLI1", "MLI2", "Oligodendrocytes", "Polydendrocytes", "Purkinje", "UBCs")
restrict_test_ref <- create_downsampled_data(test_reference, NULL, cell_type_import = common_cell_types, save.file = F)
puck_test <- seurat.to.slideseq(restrict_test_ref, cell_type_info)
gene_list = get_de_genes(cell_type_means, puck_test, fc_thresh = config$fc_cutoff, expr_thresh = config$gene_cutoff)

bulk_vec = rowSums(puck_test@counts)
nUMI = sum(bulk_vec)
X = as.matrix(cell_type_means[gene_list,] * nUMI)
b = bulk_vec[gene_list]

write.csv(as.matrix(X),file.path(bulkdir,"X_bulk.csv"))
write.csv(as.matrix(b),file.path(bulkdir,"b_bulk.csv"))

#debugging
pred = X %*% (proportions * weight_sum)
sum((log(pred,2) - log(b,2))^2)
pred_true = X %*% (true_proportions)
pred_true <- pred_true
sum((log(pred_true,2) - log(b,2))^2)
diff = (log(pred_true) - log(b))^2
diff = diff[order(diff),]
cell_type_info_viz <- get_cell_type_info(test_reference@assays$RNA@counts, test_reference@meta.data$liger_ident_coarse, test_reference@meta.data$nUMI)
cell_type_means_viz = cell_type_info_viz[[1]]
X_viz = as.matrix(cell_type_means_viz[gene_list,] * nUMI)
pred_know = X_viz %*% (true_proportions[cell_type_info_viz[[2]]])
sum((log(pred_know,2) - log(b,2))^2)

diff = (log(pred_true) - log(b))
diff = diff[order(diff),]
bad_genes = names(tail(diff,100))
write.csv(bad_genes,"genes.csv")

gene_list = get_de_genes(cell_type_means, puck_test, fc_thresh = 2, expr_thresh = config$gene_cutoff)
gene_list = get_de_genes(cell_type_means, puck_test, fc_thresh = 0.1, expr_thresh = .000003)
marker_data = get_marker_data(cell_type_names, reference, cell_type_means, gene_list, marker_provided = TRUE)


true_proportions <- aggregate(puck_test@nUMI, list(puck_test@cell_labels),sum)$x
names(true_proportions) <- common_cell_types
true_proportions <- true_proportions / sum(true_proportions)
true_proportions <- true_proportions[cell_type_names]
true_proportions[is.na(true_proportions)] <- 0
names(true_proportions) <- cell_type_names
gene_means <- as.matrix(cell_type_means[gene_list,]) %*% true_proportions
true_platform_effect = log(bulk_vec[gene_list] / nUMI,2) -log(gene_means,2)
#true_platform_effect = log((bulk_vec[gene_list] + 1) / nUMI / (gene_means+1/nUMI),2)
hist(true_platform_effect, breaks =30,xlab = "log2(Platform Effect)", main = "Measured Platform effects between dropviz and 10x")
avg_platform_effect <- aggregate(true_platform_effect, list(marker_data$cell_type), mean)
plot(avg_platform_effect$Group.1, avg_platform_effect$V1)
rownames(avg_platform_effect) <- avg_platform_effect$Group.1
new_norm <- avg_platform_effect[marker_data$cell_type,]
rownames(new_norm) <- gene_list
true_residuals <- true_platform_effect - new_norm$V1

avg_est_platform_effect <- aggregate(estimated_platform_effect, list(marker_data$cell_type), mean)
plot(avg_est_platform_effect$Group.1, avg_est_platform_effect$V1)
rownames(avg_est_platform_effect) <- avg_est_platform_effect$Group.1
new_norm <- avg_est_platform_effect[marker_data$cell_type,]
rownames(new_norm) <- gene_list
est_residuals <- estimated_platform_effect - new_norm$V1
ordered_true = abs(true_residuals)[order(abs(true_residuals))]
ordered_est = abs(est_residuals)[order(abs(est_residuals))]
plot(1:length(gene_list),ordered_true,type="n")
lines(1:length(gene_list),ordered_true,col="red")
#plot(1:length(gene_list),ordered_est,type="n")
lines(1:length(gene_list),ordered_est)
ordered_true = abs(true_platform_effect)[order(abs(true_platform_effect))]
ordered_est = abs(estimated_platform_effect)[order(abs(estimated_platform_effect))]

cell_type = "Ependymal"
hist(true_platform_effect[marker_data$cell_type == cell_type],breaks = 50)
plot(density(true_platform_effect[marker_data$cell_type == cell_type]))

hist(estimated_platform_effect[marker_data$cell_type == cell_type],breaks = 50)
plot(density(estimated_platform_effect[marker_data$cell_type == cell_type]))

prepareBulkData(bulkdir, cell_type_means, puck_test, gene_list)
#TODO: run the python code here (might have to make an individual script)
proportions <- read.csv(file.path(bulkdir,"weights.csv"))$Weight
weight_sum <- sum(proportions)
names(proportions) = cell_type_names
proportions = proportions / weight_sum
plot(true_proportions, proportions, ylab = "estimated_proportions", main = "Accuracy of predicted cell type composition")
gene_means_est <- as.matrix(cell_type_means[gene_list,]) %*% proportions
estimated_platform_effect <- log(bulk_vec[gene_list] / nUMI / gene_means_est,2)
residual_pe <- true_platform_effect - estimated_platform_effect
#mean should now be zero
hist(residual_pe, breaks =30,xlab = "log2(Platform Effect)", main = "Remaining Residual in Platform Effect Estimates")
print(paste("Mean Platform Effect:", mean(residual_pe)))
print(paste("Std Platform Effect:",  sqrt(var(residual_pe))))

#fit a linear function to predict true platform effect from estimated platform effect
#intercept should be 0
#this will give us the shrinkage
plot(estimated_platform_effect, true_platform_effect)
pe_df <- data.frame(estimated_platform_effect, true_platform_effect)
reg <- lm(true_platform_effect~estimated_platform_effect,data=pe_df)
shrinkage <- reg$coefficients[2]
summary(reg)
#these three things should hopefully be the same estimate of bias:
mean_val <- mean(estimated_platform_effect)
print(paste("Bias Estimate 1 (mean PE):",mean_val))
print(paste("Bias Estimate 2 (true PE of mean PE):",reg$coefficients[1] + shrinkage*mean_val))
print(paste("Bias estimate 3 (weight sum)", log(weight_sum,2)))

with(pe_df,plot(estimated_platform_effect, true_platform_effect))
abline(reg)

final_pe_estimate <- estimated_platform_effect * shrinkage
residual_pe_shr <- true_platform_effect - final_pe_estimate
hist(residual_pe_shr, breaks =30,xlab = "log2(Platform Effect)", main = "Remaining Residual in Platform Effect Estimates with Shrinkage")
print(paste("Mean Platform Effect after Shrinkage:", mean(residual_pe_shr)))
print(paste("Std Platform Effect after Shrinkage:",  sqrt(var(residual_pe_shr))))


#scratch
gene = "Ptgds"
plot(log(unlist(cell_type_means_viz[gene,common_cell_types]) + 1e-9), log(unlist(cell_type_means[gene,common_cell_types]) + 1e-9))

cell_type = "Granule"
bulk_vec = rowSums(puck_test@counts[,puck_test@cell_labels==cell_type])
cur_gene_list = names(which(bulk_vec[gene_list] >= 3))
nUMI = sum(bulk_vec)
X = as.matrix(cell_type_means[cur_gene_list,] * nUMI)
b = bulk_vec[cur_gene_list]
write.csv(as.matrix(X),file.path(bulkdir,"X_bulk.csv"))
write.csv(as.matrix(b),file.path(bulkdir,"b_bulk.csv"))

gene_list = get_de_genes(cell_type_means, puck_test, fc_thresh = 0.75, expr_thresh = .0002)
test_results = process_data(puck_test, gene_list, cell_type_info, proportions, trust_model = F, constrain = F, OLS = F)
plot_heat_map(test_results[[1]]$table)
weights <- t(as.data.frame(matrix(unlist(test_results[[2]]), nrow=n_cell_types)))
rownames(weights) = colnames(puck_test@counts)
colnames(weights) = cell_type_names
norm_weights = sweep(weights, 1, rowSums(weights), '/')
avg_weights <- aggregate.data.frame(as.data.frame(weights), list(puck_test@cell_labels), mean)
rownames(avg_weights) <- common_cell_types
avg_weights <- avg_weights[cell_type_names,2:dim(avg_weights)[2]]
avg_weights[is.na.data.frame(avg_weights)] <- 0
rownames(avg_weights) <- cell_type_names
avg_weights <- t(as.matrix(avg_weights))
avg_weights_norm <- aggregate.data.frame(as.data.frame(norm_weights), list(puck_test@cell_labels), mean)
plot_heat_map(avg_weights)

#make simulated doublets (50/50 mixture, 1000 total UMI)

#scratch one of each
nUMI = test_reference@meta.data$nUMI
firstInd = sample(intersect(which(test_reference@meta.data$liger_ident_coarse == type1) , which(nUMI > 1000)),1)
secondInd = sample(intersect(which(test_reference@meta.data$liger_ident_coarse == type2) , which(nUMI > 1000)),1)
bead1 = as.matrix(sub_sample_cell(iv$gene_list, test_reference@assays$RNA@counts, firstInd, 500))[,1]
bead2 = as.matrix(sub_sample_cell(iv$gene_list, test_reference@assays$RNA@counts, secondInd, 500))[,1]
results = process_bead(cell_type_info_renorm, iv$gene_list, 1000, bead1 + bead2, constrain = F)
print(results)
bead1 = test_reference@assays$RNA@counts[iv$gene_list,firstInd]
results1 = process_bead(cell_type_info_renorm, iv$gene_list, test_reference@meta.data$nUMI[firstInd], bead1, constrain = F)
results2 = process_bead(cell_type_info_renorm, iv$gene_list, 500, bead2, constrain = F)
#scratch quick test
type1 = "Fibroblast"
type2 = "Granule"
cell_type_info_renorm = iv$cell_type_info
cell_type_info_renorm[[1]] = get_norm_ref(puck_test, cell_type_info[[1]], gene_list, proportions)
cell_type_means_renorm <- cell_type_info_renorm[[1]]
weight_recovery_test(test_reference, iv$gene_list, cell_type_means_renorm, type1, type2, conditions = 5, trials = 10)
doublet_accuracy_test(test_reference, gene_list, cell_type_info_renorm, type1, type2, conditions = 10, trials = 20)
bead = bead_mix(test_reference, gene_list, 500, 500, type1, type2)
decompose(cell_type_means_renorm, gene_list, 1000, bead, constrain=F)
decompose_sparse(cell_type_means_renorm, gene_list, 1000, bead, type1, type2, score_mode = T)
decompose_sparse(cell_type_means_renorm, gene_list, 1000, bead, "MLI2", type1, score_mode = T)
decompose_sparse_score(cell_type_means_renorm, gene_list, 1000, bead, type1, type2)
decompose_sparse_score(cell_type_means_renorm, gene_list, 1000, bead, type1, "MLI2")
decompose_sparse(cell_type_means_renorm, gene_list, 1000, bead, "MLI2", type1, score_mode = T)
weights= solveIRWLS.weights(reg_data,bead,nUMI,OLS = FALSE, constrain = F)

#here we plot the quasiliklihood for Y = 1
Q_mat <- get_Q_mat()
bead = bead_mix(test_reference, gene_list, 500,500, type1, type2)
results = process_bead(cell_type_info_renorm, gene_list, 1000, bead, constrain = F)
print(results)
weights <- decompose_sparse(cell_type_means_renorm, gene_list, 1000, bead, type1, type2, score_mode = F, constrain = F, verbose = T)
print(weights)
cell_types = c(type1,type2)
reg_data = cell_type_means_renorm[gene_list,] * 1000
reg_data = data.matrix(reg_data)
reg_data = reg_data[,cell_types]
prediction = reg_data %*% weights

weights <- decompose_sparse(cell_type_means_renorm, gene_list, 1000, bead, type1, "Astrocytes", score_mode = F, constrain = F, verbose = T)
cell_types = c(type1,"MLI2")
reg_data = cell_type_means_renorm[gene_list,] * 1000
reg_data = data.matrix(reg_data)
reg_data = reg_data[,cell_types]
prediction = reg_data %*% weights




Y = 0; n = 1000003
print(J_mat[Y+1,n])
print(get_QL(Y,n*delta))

plot(x,J_mat[1,],type="n",xlim=c(0,10))
lines(x,J_mat[1,])
bead = bead_mix(test_reference, gene_list, 500,500, type1, type2)

#weights <- decompose_sparse(cell_type_means_renorm, gene_list, 1000, bead, type1, type2, score_mode = F, constrain = F)
print(weights)

bead = bead_mix(test_reference, gene_list, 300,700, type1, type2)
UMI_tot = 1000


cell_types = c(type1,type2)
reg_data = cell_type_means_renorm[gene_list,] * UMI_tot
reg_data = data.matrix(reg_data)
reg_data = reg_data[,cell_types]
weights= solveIRWLS.weights(reg_data,bead,nUMI,OLS = FALSE, constrain = F)$weights
J = 0
J_list = numeric(length(gene_list))
names(J_list) = gene_list
prediction = reg_data %*% weights
for(gene in gene_list) {
  Y = bead[gene,]
  J = J + get_QL(Y, prediction[gene,])
  J_list[gene] = get_QL(Y, prediction[gene,])
}
scaled_residual = (prediction - bead)/(sqrt(get_V(prediction)))
my_score = sum(abs(scaled_residual))

#weights <- decompose_sparse(cell_type_means_renorm, gene_list, 1000, bead, type2, "MLI2", score_mode = F, constrain = F)
cell_types = c(type2,"MLI2")
reg_data = cell_type_means_renorm[gene_list,] * UMI_tot
reg_data = data.matrix(reg_data)
reg_data = reg_data[,cell_types]
#weights= solveIRWLS.weights(reg_data,bead,nUMI,OLS = FALSE, constrain = constrain, verbose = verbose)
weights= solveIRWLS.weights(reg_data,bead,nUMI,OLS = FALSE, constrain = F)$weights
prediction_bad = reg_data %*% weights
#print(my_score)
print(J)
J = 0
J_list_bad = numeric(length(gene_list))
names(J_list_bad) = gene_list
for(gene in gene_list) {
  Y = bead[gene,]
  J = J + get_QL(Y, prediction_bad[gene,])
  J_list_bad[gene] = get_QL(Y, prediction_bad[gene,])
}
print(J)
scaled_residual_bad = (prediction_bad - bead)/(sqrt(get_V(prediction_bad)))
my_score = sum(abs(scaled_residual_bad))
#print(my_score)
scaled_residual <-
score_df <- data.frame(J_list,J_list_bad,scaled_residual,scaled_residual_bad)
abs_score_df <- abs(score_df)
abs_score_df <- abs_score_df[order(abs_score_df$J_list),]
abs_score_df[gene_list,"prediction"] <- prediction[gene_list,]
abs_score_df[gene_list,"prediction_bad"] <- prediction_bad[gene_list,]
abs_score_df[gene_list,"Y"] <- bead[gene_list,]
print(colSums(abs_score_df))
head_df <- head(abs_score_df,400)
head_df <- head_df[order(head_df$scaled_residual_bad),]
tail(head_df)
norm_l =colSums(abs_score_df)
norm_l["J_list_bad"] = norm_l["J_list"]
norm_l["scaled_residual_bad"] = norm_l["scaled_residual"]
my_df <- abs_score_df[abs_score_df$Y == 0 & ((abs_score_df$prediction > 0.05) | (abs_score_df$prediction_bad > 0.05)),]
my_df <- abs_score_df[abs_score_df$Y > 0 & ((abs_score_df$prediction < 0.5) | (abs_score_df$prediction_bad < 0.5)),]
print(colSums(my_df)/norm_l)
