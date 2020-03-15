library(readr)
library(Matrix)
library(RCTD)
dge <- read_csv('Data/Reference/Liver/GSM3714747_chow1_filtered_gene_bc_matrices.csv')
gene_list <- dge$X1
dgeMat <- as(as(dge[,-1],"matrix"),"dgCMatrix")
rownames(dgeMat) <- gene_list
colnames(dgeMat) <- lapply(colnames(dgeMat), function (x) paste0("chow1_",x))

dge <- read_csv('Data/Reference/Liver/GSM3714748_chow2_filtered_gene_bc_matrices.csv')
dgeMat2 <- as(as(dge[,-1],"matrix"),"dgCMatrix")
rownames(dgeMat2) <- dge$X1
colnames(dgeMat2) <- lapply(colnames(dgeMat2), function (x) paste0("chow2_",x))

dge <- read_csv('Data/Reference/Liver/GSM3714749_chow3_filtered_gene_bc_matrices.csv')
dgeMat3 <- as(as(dge[,-1],"matrix"),"dgCMatrix")
rownames(dgeMat3) <- dge$X1
colnames(dgeMat3) <- lapply(colnames(dgeMat3), function (x) paste0("chow3_",x))

counts <- cbind(dgeMat,dgeMat2,dgeMat3)
barcodes <- colnames(counts)
nUMI = meta_data[barcodes,"nUMI"]
cluster = meta_data[barcodes,"res.0.07"]
cell_type_names = c('Endo', 'Macrophage', 'T_Cell', 'B_Cell', 'Dendtritic', 'Cholangiocyte', 'Hepatocyte', 'Dividing', 'Plasma_B','HSC')

meta_data = as.data.frame(cluster)
meta_data$nUMI = nUMI
refdir = "Data/Reference/Liver"
cell_dict_file <- file.path(refdir,"cell_type_dict.csv")
true_type_names <- remap_celltypes(cell_dict_file, meta_data$cluster)
meta_data$liger_ident_coarse = true_type_names
rownames(meta_data) = barcodes
reference = Seurat::CreateSeuratObject(counts, meta.data = meta_data)
saveRDS(reference, paste(refdir,"SCRef.RDS",sep="/"))
dref <- create_downsampled_data(reference, refdir)
dref <- create_downsampled_data(dref, refdir, n_samples = 300)
create_downsampled_data(dref, refdir, n_samples = 25)
