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
test_reference <- readRDS("Data/Reference/DropVizCerAnnotated/scRefSubsampled25.RDS")
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
