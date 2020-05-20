#simulates doublets for testing on scRNA and snRNA cerebellum reference.
library(RCTD)
library(Matrix)
DropViz <- T
if(DropViz) {
  test_reference <- readRDS("Data/Reference/DropVizCerAnnotated/scRefSubsampled1000.RDS")
  common_cell_types = c("Astrocytes", "Bergmann", "Endothelial", "Fibroblast", "Golgi", "Granule", "MLI1", "MLI2", "Oligodendrocytes", "Polydendrocytes", "Purkinje", "UBCs")
  restrict_test_ref <- create_downsampled_data(test_reference, NULL, cell_type_import = common_cell_types, save.file = F)
} else {
  restrict_test_ref <- readRDS("Data/Reference/10xCer/scRefSubsampled1000.RDS")
}
iv <- init_RCTD(gene_list_reg = F, test_reference = restrict_test_ref)
if(!DropViz)
  common_cell_types <- iv$cell_type_info[[2]]
metadir <- file.path(iv$slideseqdir,"MetaData")
if(!dir.exists(metadir))
  dir.create(metadir)

#next, generate the fake dataset
n_cell_types = length(common_cell_types)
trials = 30 # 30
n_conditions = 13 # 10
N_samples = (n_cell_types * trials * n_conditions * (n_cell_types - 1))/2
first_UMI = numeric(N_samples); first_type = character(N_samples); second_type = character(N_samples)
UMI_tot = 1000; UMI_step = UMI_tot / (n_conditions-1)
UMI_tot = round(UMI_step * (n_conditions-1)); UMI1_vec = round(0:(n_conditions-1)*UMI_step)
beads = Matrix(0, nrow = length(iv$gene_list), ncol = N_samples)
rownames(beads) = iv$gene_list; colnames(beads) = 1:N_samples
firstbeads = Matrix(0, nrow = length(iv$gene_list), ncol = N_samples)
rownames(firstbeads) = iv$gene_list; colnames(firstbeads) = 1:N_samples
secbeads = Matrix(0, nrow = length(iv$gene_list), ncol = N_samples)
rownames(secbeads) = iv$gene_list; colnames(secbeads) = 1:N_samples
index = 1
nUMI = restrict_test_ref@meta.data$nUMI
for(i in 1:(n_cell_types-1)) {
  print(paste("Progress",i))
  for(j in (i+1):n_cell_types) {
    print(paste("ProgressSecond",j))
    type1 = common_cell_types[i]; type2 = common_cell_types[j]
    for (condition in 1:n_conditions) {
      UMI1 = UMI1_vec[condition]; UMI2 = UMI_tot - UMI1
      for(t in 1:trials) {
        first_UMI[index] = UMI1; first_type[index] = type1; second_type[index] = type2
        firstInd = sample(intersect(which(restrict_test_ref@meta.data$liger_ident_coarse == type1) , which(nUMI > 1000)),1)
        secondInd = sample(intersect(which(restrict_test_ref@meta.data$liger_ident_coarse == type2) , which(nUMI > 1000)),1)
        firstbeads[,index] = as.vector(sub_sample_cell(iv$gene_list, restrict_test_ref@assays$RNA@counts, firstInd, UMI1))
        secbeads[,index] = as.vector(sub_sample_cell(iv$gene_list, restrict_test_ref@assays$RNA@counts, secondInd, UMI2))
        beads[,index] = firstbeads[,index] + secbeads[,index]
        index = index + 1
      }
    }
  }
}
puck <- Slideseq(NULL, beads, nUMI = rep(UMI_tot,N_samples))
saveRDS(puck, file.path(iv$slideseqdir, iv$puckrds))
meta_df <- data.frame(first_type, second_type, first_UMI, row.names = colnames(beads))
meta_data <- list(n_cell_types = n_cell_types, trials = trials, n_conditions = n_conditions, N_samples = N_samples,
                  meta_df = meta_df, UMI_tot = UMI_tot, UMI_list = UMI1_vec)
saveRDS(meta_data,file.path(metadir,"meta_data.RDS"))
saveRDS(firstbeads, file.path(metadir,"firstbeads.RDS"))
saveRDS(secbeads, file.path(metadir,"secbeads.RDS"))
