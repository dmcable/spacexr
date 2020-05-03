library(readr)
library(RCTD)
library(Seurat)
config_data <- config::get(file = "conf/dataset.yml", use_parent = FALSE)
config <- config::get(file = paste0("conf/",config_data$config_mode,".yml"), use_parent = FALSE)
slideseqdir <- file.path("Data/Slideseq",config_data$slideseqfolder)
resultsdir = file.path(slideseqdir,"results")
refdir <- file.path("Data/Reference",config_data$reffolder)
reference <- readRDS(paste(refdir,config_data$reffile,sep="/"))
puck = readRDS(file.path(slideseqdir, config_data$puckrds))
puck <- restrict_counts(puck, rownames(puck@counts), UMI_thresh = 100)
my_ref <- CreateSeuratObject(reference@assays$RNA@counts,project = "SeuratProject",assay = "RNA",min.cells = 0,min.features = 0,names.field = 1,names.delim = "_",meta.data = NULL)
Idents(my_ref) <- reference@meta.data$liger_ident_coarse
my_ref <- NormalizeData(my_ref, verbose = FALSE)
my_ref <- FindVariableFeatures(my_ref, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

target <- CreateSeuratObject(puck@counts, project = "Slideseq",assay = "RNA",min.cells = 0,min.features = 0,names.field = 1,names.delim = "_",meta.data = NULL)
target <- NormalizeData(target, verbose = FALSE)
target <- FindVariableFeatures(target, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
resultsdir = "Data/Slideseq/NewCerPuck_190926_08/SeuratResults/"
target <- ScaleData(target)
target <- RunPCA(target, assay = "RNA", verbose = FALSE)
target <- RunUMAP(target, dims = 1:30)
target <- FindNeighbors(target, dims = 1:30)
target <- FindClusters(target, resolution = 0.3, verbose = FALSE)
puck@cell_labels <- target@active.ident
n_clusters = length(levels(puck@cell_labels))
pdf(file.path(resultsdir,'cluster_images.pdf'))
prop_mat <- Matrix(0,nrow=8,ncol = 19)
colnames(prop_mat) <- iv$cell_type_info[[2]]
for(i in 0:(n_clusters - 1)) {
  prop_mat[i+1,] <- as.vector(table(results_df[results_df$spot_class!="reject" & puck@cell_labels==i,"first_type"]) + table(results_df[results_df$spot_class=="doublet_certain" & puck@cell_labels==i,"second_type"]))
  #invisible(print(plot_puck_wrapper(puck, puck@nUMI, cell_type = i, minUMI = 0, maxUMI = 200000, min_val = 0, max_val = 2000)))
}
dev.off()
puck_seurat <- puck
saveRDS(puck_seurat, file.path("Data/Slideseq/NewCerPuck_190926_08/SeuratResults",'results_df.RDS'))


#Finding DE genes / marker genes
cur_cell_types <- c("Bergmann","Granule","Purkinje","MLI1","Oligodendrocytes")
de_genes <- get_de_genes(cell_type_info_restr, puck, fc_thresh = 3, expr_thresh = .001, MIN_OBS = 3)
