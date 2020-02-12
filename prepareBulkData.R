#this script takes the single cell reference and puck and prepares the bulk data files so that you can run fitBulkWeights.py
#working directory must be the Cell Demixing folder.
library(RCTD)
library(Matrix)
print("prepareBulkData: begin")
config <- config::get()
refdir <- file.path("Data/Reference",config$reffolder)
dir.create(file.path(refdir, "results"), showWarnings = FALSE) #folder to save results
slideseqdir <- file.path("Data/Slideseq",config$slideseqfolder)
resultsdir = file.path(slideseqdir,"results")
if(!dir.exists(resultsdir))
  dir.create(resultsdir)
bulkdir <- paste(slideseqdir,"results/Bulk",sep="/")
if(!dir.exists(bulkdir))
  dir.create(bulkdir)
puck = read.slideseq(slideseqdir, count_file = config$puckfile)
puck = restrict_counts(puck, rownames(puck@counts), UMI_thresh = config$UMI_min)
reference <- readRDS(paste(refdir,config$reffile,sep="/"))
print(paste("prepareBulkData: number of cells in reference:", dim(reference@assays$RNA@counts)[2]))
print(paste("prepareBulkData: number of genes in reference:", dim(reference@assays$RNA@counts)[1]))
cell_type_info <- get_cell_type_info(reference@assays$RNA@counts, reference@meta.data$liger_ident_coarse, reference@meta.data$nUMI)
cell_type_means = cell_type_info[[1]]; cell_type_names = cell_type_info[[2]]
print("prepareBulkData: amount of each cell type in the reference:")
cell_counts = table(reference@meta.data$liger_ident_coarse)
print(cell_counts)
CELL_MIN = 25 # need at least this for each cell type
if(min(cell_counts) < CELL_MIN)
  stop(paste0("prepareBulkData error: need a minimum of ",CELL_MIN, " cells for each cell type in the reference"))
print(paste("prepareBulkData: number of cell types used:", cell_type_info[[3]]))
bulk_vec = rowSums(puck@counts)
print(paste("prepareBulkData: number of spots used in test data passing UMI threshold:", dim(puck@counts)[2]))
gene_list = get_de_genes(cell_type_means, puck, fc_thresh = config$fc_cutoff, expr_thresh = config$gene_cutoff)
print(paste("prepareBulkData: number of genes used for Platform Effect Estimation:", length(gene_list)))
nUMI = sum(bulk_vec)
X = cell_type_means[gene_list,] * nUMI
b = bulk_vec[gene_list]
write.csv(as.matrix(X),file.path(bulkdir,"X_bulk.csv"))
write.csv(as.matrix(b),file.path(bulkdir,"b_bulk.csv"))
print("prepareBulkData: end")
