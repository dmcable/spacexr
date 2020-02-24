library(RCTD)
library(Matrix)
test_reference <- readRDS("Data/Reference/DropVizCerAnnotated/scRefSubsampled1000.RDS")
common_cell_types = c("Astrocytes", "Bergmann", "Endothelial", "Fibroblast", "Golgi", "Granule", "MLI1", "MLI2", "Oligodendrocytes", "Polydendrocytes", "Purkinje", "UBCs")
restrict_test_ref <- create_downsampled_data(test_reference, NULL, cell_type_import = common_cell_types, save.file = F)
iv <- init_RCTD(test_reference = restrict_test_ref, get_proportions = T)
metadir <- file.path(iv$slideseqdir,"MetaData")
#generate the celltypeinfo-renorm
cell_type_info_renorm = iv$cell_type_info
cell_type_info_renorm[[1]] = get_norm_ref(iv$puck, iv$cell_type_info[[1]], iv$gene_list, iv$proportions)
saveRDS(cell_type_info_renorm, file.path(metadir, "cell_type_info_renorm.RDS"))

#next, generate the fake dataset
n_cell_types = 6 # length(common_cell_types)
trials = 10 # 30
n_conditions = 5 # 10
N_samples = (n_cell_types * trials * n_conditions * (n_cell_types - 1))/2
first_UMI = numeric(N_samples); first_type = character(N_samples); second_type = character(N_samples)
UMI_tot = 1000; UMI_step = round(UMI_tot / (n_conditions-1))
UMI_tot = UMI_step * (n_conditions-1); UMI1_vec = 0:(n_conditions-1)*UMI_step
beads = Matrix(0, nrow = length(iv$gene_list), ncol = N_samples)
rownames(beads) = iv$gene_list; colnames(beads) = 1:N_samples
index = 1
for(i in 1:(n_cell_types-1))
  for(j in (i+1):n_cell_types) {
    type1 = common_cell_types[i]; type2 = common_cell_types[j]
    for (condition in 1:n_conditions) {
      UMI1 = UMI1_vec[condition]; UMI2 = UMI_tot - UMI1
      for(t in 1:trials) {
        first_UMI[index] = UMI1; first_type[index] = type1; second_type[index] = type2
        beads[,index] = bead_mix(restrict_test_ref, iv$gene_list, UMI1, UMI2, type1, type2)
        index = index + 1
      }
    }
  }
puck <- Slideseq(NULL, beads, nUMI = rep(UMI_tot,N_samples))
split_puck(puck, iv$slideseqdir, iv$config$n_puck_folds)
meta_df <- data.frame(first_type, second_type, first_UMI, row.names = colnames(beads))
meta_data <- list(n_cell_types = n_cell_types, trials = trials, n_conditions = n_conditions, N_samples = N_samples,
                  meta_df = meta_df, UMI_tot = UMI_tot, UMI_list = UMI1_vec)
saveRDS(meta_data,file.path(metadir,"meta_data.RDS"))
