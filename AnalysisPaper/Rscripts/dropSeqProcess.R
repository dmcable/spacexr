library(DropSeq.util)
refdir = "Data/Reference/DropVizCerAnnotated"
dge.path <- file.path(refdir,"F_GRCm38.81.P60Cerebellum_ALT.raw.dge.txt.gz")
class_path <- file.path(refdir,"F_GRCm38.81.P60Cerebellum_ALT.subcluster.assign.RDS")
dge <- loadSparseDge(dge.path)
cluster <- readRDS(class_path)
common_barcodes = intersect(names(cluster), colnames(dge))
raw.data = dge[,common_barcodes]
cluster = cluster[common_barcodes]
meta_data = as.data.frame(cluster)
meta_data$nUMI = colSums(raw.data)
cell_dict_file <- file.path(refdir,"cell_type_dict.csv")
true_type_names <- remap_celltypes(cell_dict_file, meta_data$cluster)
meta_data$liger_ident_coarse = true_type_names
reference = Seurat::CreateSeuratObject(raw.data, meta.data = meta_data)
saveRDS(reference, paste(refdir,"SCRef.RDS",sep="/"))
dref <- create_downsampled_data(reference, refdir)
dref <- create_downsampled_data(dref, refdir, n_samples = 300)
create_downsampled_data(dref, refdir, n_samples = 25)
