#input: refdir, the directory of the reference
#e.g. dgeToSeurat('Data/Reference/KidneyHumanReference')
dgeToSeurat <- function(refdir) {
  raw.data = t(read.csv(file.path(refdir,"dge.csv")))
  meta_data = read.csv(file.path(refdir,"meta_data.csv"))
  rownames(meta_data) = meta_data$barcode
  meta_data$barcode = NULL
  common_barcodes = intersect(colnames(raw.data), rownames(meta_data))
  raw.data = raw.data[,common_barcodes]
  meta_data = meta_data[common_barcodes, ]
  cell_dict_file <- paste(refdir,"cell_type_dict.csv",sep="/")
  true_type_names <- remap_celltypes(cell_dict_file, meta_data$cluster)
  meta_data$liger_ident_coarse = true_type_names
  reference = Seurat::CreateSeuratObject(raw.data, meta.data = meta_data)
  saveRDS(reference, paste(refdir,"SCRef.RDS",sep="/"))
}
