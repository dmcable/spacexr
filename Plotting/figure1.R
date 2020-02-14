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
print(paste("prepareBulkData: number of cells in reference:", dim(reference@assays$RNA@counts)[2]))
print(paste("prepareBulkData: number of genes in reference:", dim(reference@assays$RNA@counts)[1]))
print("prepareBulkData: amount of each cell type in the reference:")
cell_counts = table(reference@meta.data$liger_ident_coarse)
puck <- seurat.to.slideseq(reference, cell_type_info)
gene_list = get_de_genes(cell_type_means, puck, fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg)
print(paste("callBeads: number of genes used for regression:", length(gene_list)))
#UMI_min = 500
UMI_min = 100
puck = restrict_counts(puck, gene_list, UMI_thresh = UMI_min, UMI_max = 10000000)
puck = restrict_puck(puck, colnames(puck@counts))
table(puck@cell_labels)
test_results = process_data(puck, gene_list, cell_type_info, proportions = NULL, trust_model = F, constrain = T)
conf_mat = test_results[[1]]; weights = test_results[[2]]; pred_labels = test_results[[3]]
plot_heat_map(conf_mat$table)

library(ggplot2)
decompose(cell_type_means, gene_list, puck@nUMI[1], puck@counts[, 1], constrain = T)
doublet_df <- test_doublet_beads(puck, gene_list, cell_type_info)
doublet_df$true_types = puck@cell_labels
plot_df = aggregate(doublet_df[,"singlet"], list(doublet_df$true_types), mean)
colnames(plot_df) = c('Cell.Type', 'Singlet.Proportion')
#figure 1b: Singlet identification
gg <- ggplot(plot_df, ggplot2::aes(x=Cell.Type, y=Singlet.Proportion)) + geom_point() + ylim(0,1)
gg

#figure 1c: singlet mis-asignment
doublet_df[doublet_df$singlet,"max_weight"] <- 1
avg_weights <- matrix(0,ncol=n_cell_types, nrow=n_cell_types,byrow=TRUE)
colnames(avg_weights) <- cell_type_names
rownames(avg_weights) <- cell_type_names
#smoke <- as.table(smoke)

for (i in 1:dim(doublet_df)[1]) {
  true_type <- doublet_df[i,"true_types"]
  first_type <- doublet_df[i,"first_label"]
  sec_type <- doublet_df[i,"sec_label"]
  avg_weights[first_type, true_type] <- avg_weights[first_type, true_type] + doublet_df[i,"max_weight"]
  avg_weights[sec_type, true_type] <- avg_weights[sec_type, true_type] + 1 - doublet_df[i,"max_weight"]
}
plot_heat_map(avg_weights)

