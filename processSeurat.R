library(readr)
config_data <- config::get(file = "conf/dataset.yml", use_parent = FALSE)
config <- config::get(file = paste0("conf/",config_data$config_mode,".yml"), use_parent = FALSE)
slideseqdir <- file.path("Data/Slideseq",config_data$slideseqfolder)
resultsdir = file.path(slideseqdir,"results")
outdir = 'NMFreg/DataCer'
refdir <- file.path("Data/Reference",config_data$reffolder)
reference <- readRDS(paste(refdir,config_data$reffile,sep="/"))
puck = readRDS(file.path(slideseqdir, config_data$puckrds))
puck <- restrict_counts(puck, rownames(puck@counts), UMI_thresh = 100)
all_ref <- readRDS('Seurat/allen_cortex.rds')
my_ref <- CreateSeuratObject(reference@assays$RNA@counts,project = "SeuratProject",assay = "RNA",min.cells = 0,min.features = 0,names.field = 1,names.delim = "_",meta.data = NULL)
Idents(my_ref) <- reference@meta.data$liger_ident_coarse
my_ref <- NormalizeData(my_ref, verbose = FALSE)
my_ref <- FindVariableFeatures(my_ref, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

target <- CreateSeuratObject(puck@counts, project = "Slideseq",assay = "RNA",min.cells = 0,min.features = 0,names.field = 1,names.delim = "_",meta.data = NULL)
target <- NormalizeData(target, verbose = FALSE)
target <- FindVariableFeatures(target, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

pancreas.anchors <- FindTransferAnchors(reference = my_ref, query = target, dims = 1:30)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = my_ref@active.ident, dims = 1:30)
write_rds(predictions, file.path(resultsdir,'predictions.RDS'))
table(predictions$predicted.id)
weights_seurat <- predictions[,2:(cell_type_info[[3]]+1)]
colnames(weights_seurat) <- cell_type_info[[2]]
resultsdir = "Data/Slideseq/NewCerPuck_190926_08/SeuratResults/"
dir.create(resultsdir)
thresh <- rep(0.25, cell_type_info[[3]]); names(thresh) <- cell_type_info[[2]]
plot_weights_nmf(iv$cell_type_info, puck, resultsdir, weights_seurat, thresh)
plot_weights_unthreshold(iv$cell_type_info, puck, resultsdir, weights_seurat)
plot_cond_occur_nmf(cell_type_info, resultsdir, weights_seurat, thresh)
plot_occur_unthreshold(cell_type_info, resultsdir, weights_seurat)
plot_val <- weights_seurat$Microglia + weights_seurat$MLI1 + weights_seurat$MLI2
names(plot_val) <- colnames(puck@counts)
plot_puck_continuous(puck, colnames(puck@counts), plot_val)
gene_list = get_de_genes(cell_type_info, puck, fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg, MIN_OBS = 3)
RCTD::get_de_genes(cell_type_info, puck, fc_thresh = 2, expr_thresh = config$gene_cutoff_reg, MIN_OBS = 3)
