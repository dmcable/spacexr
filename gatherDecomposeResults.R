library(RCTD)
library(Matrix)
library(dplyr)
library(ggplot2)
require(reshape2)
source('Plotting/figure_utils.R')
iv <- init_RCTD(load_info_renorm = T) #initial variables
cell_type_names = iv$cell_type_info[[2]]
metadir <- file.path(iv$slideseqdir,"MetaData")
meta_data <- readRDS(file.path(metadir,"meta_data.RDS"))
meta_df <- meta_data$meta_df
barcodes <- rownames(meta_df)
N = length(barcodes)
first_beads <- readRDS(file.path(metadir,"firstbeads.RDS"))
second_beads <- readRDS(file.path(metadir,"secbeads.RDS"))
expect1 = Matrix(0, nrow = N, ncol = length(iv$gene_list))
expect2 = Matrix(0, nrow = N, ncol = length(iv$gene_list))
var = Matrix(0, nrow = N, ncol = length(iv$gene_list))
weights_doublet = Matrix(0, nrow = N, ncol = 2)
index = 1
expect1_list <- list(); expect2_list <- list(); var_list <- list();
for (fold_index in 1:iv$n_puck_folds) {
  print(fold_index)
  results <- readRDS(paste0(iv$slideseqdir,"/DecomposeResults/results",fold_index,".RDS"))
  weights_doublet[index:(index+length(results)-1),] <- do.call(rbind,lapply(results,function(x) x$weights))
  exp_list <- lapply(results,function(x) x$decompose_results$expect_1)
  expect1_list[[fold_index]] <- matrix(unlist(exp_list), ncol = length(iv$gene_list), byrow = TRUE)
  exp_list <- lapply(results,function(x) x$decompose_results$expect_2)
  expect2_list[[fold_index]] <- matrix(unlist(exp_list), ncol = length(iv$gene_list), byrow = TRUE)
  exp_list <- lapply(results,function(x) x$decompose_results$variance)
  var_list[[fold_index]] <- matrix(unlist(exp_list), ncol = length(iv$gene_list), byrow = TRUE)
  index = index + length(results)
}

expect1 <- do.call(rbind, expect1_list)
rm(expect1_list)
expect2 <- do.call(rbind, expect2_list)
rm(expect2_list)
var <- do.call(rbind, var_list)
rm(var_list)
rownames(expect1) = barcodes[1:N];
rownames(expect2) = barcodes[1:N]; rownames(var) = barcodes[1:N]
colnames(expect1) = iv$gene_list; colnames(expect2) = iv$gene_list
colnames(var) = iv$gene_list;
rownames(weights_doublet) = barcodes[1:N]; colnames(weights_doublet) = c('first_type', 'second_type')

beads = expect1 + expect2

#weight recovery plot
DropViz <- F
if(DropViz) {
  common_cell_types = c("Astrocytes", "Bergmann", "Endothelial", "Fibroblast", "Golgi", "Granule", "MLI1", "MLI2", "Oligodendrocytes", "Polydendrocytes", "Purkinje", "UBCs")
} else {
  common_cell_types <- iv$cell_type_info[[2]]
}
type1 = "Astrocytes"; type2 = "Bergmann"
RMSE_dat = Matrix(0, nrow = meta_data$n_cell_types, ncol = meta_data$n_cell_types)
R2_dat = Matrix(0, nrow = meta_data$n_cell_types, ncol = meta_data$n_cell_types)
rownames(RMSE_dat) = common_cell_types[1:meta_data$n_cell_types]
colnames(RMSE_dat) = common_cell_types[1:meta_data$n_cell_types]
rownames(R2_dat) = common_cell_types[1:meta_data$n_cell_types]
colnames(R2_dat) = common_cell_types[1:meta_data$n_cell_types]
plot_results <- as.list(rep(0,meta_data$n_cell_types^2))
dim(plot_results) <- c(meta_data$n_cell_types, meta_data$n_cell_types)
dimnames(plot_results) <- list(common_cell_types[1:meta_data$n_cell_types],common_cell_types[1:meta_data$n_cell_types])
for(ind1 in 1:(meta_data$n_cell_types-1))
  for(ind2 in (ind1+1):(meta_data$n_cell_types)) {
    type1 = common_cell_types[ind1]; type2 = common_cell_types[ind2]
    print(paste(type1,type2))
    plot_results[[type1,type2]]<- get_decompose_plots(meta_df, type1, type2, weights_doublet, meta_data, iv, expect1, expect2, var, first_beads, second_beads, beads)
    R2_dat[type1,type2] <- plot_results[[type1,type2]]$R2
    RMSE_dat[type1,type2] <- plot_results[[type1,type2]]$RMSE
    R2_dat[type2,type1] <- plot_results[[type1,type2]]$R2
    RMSE_dat[type2,type1] <- plot_results[[type1,type2]]$RMSE
  }
decomp_dir = file.path(iv$slideseqdir,"DecomposePlots")
if(!dir.exists(decomp_dir))
  dir.create(decomp_dir)
for(plot_title in c('bias_plot', 'err_plot','hist_plot','de_gene_plot','weight_plot', 'de_ind_plot')) {
  pdf(file.path(decomp_dir,paste0(plot_title,".pdf")))
  for(ind1 in 1:(meta_data$n_cell_types-1))
    for(ind2 in (ind1+1):(meta_data$n_cell_types)) {
      type1 = common_cell_types[ind1]; type2 = common_cell_types[ind2]
      invisible(print(plot_results[[type1,type2]][[plot_title]]))
    }
  dev.off()
}

plot_heat_map(as.matrix(R2_dat), normalize = F, file_loc = file.path(decomp_dir,"geneR2.png"), save.file=T)
plot_heat_map(as.matrix(RMSE_dat), normalize = F, file_loc = file.path(decomp_dir,"RMSE_weights.png"), save.file=T)

