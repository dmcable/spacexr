library(RCTD)
library(Matrix)
#scratch
iv <- init_RCTD(gene_list_reg = F, get_proportions = F) #initial variables
split_puck(iv$puck, iv$slideseqdir, iv$n_puck_folds)
cell_type_info_renorm = iv$cell_type_info
gene_list = get_de_genes(iv$cell_type_info, iv$puck, fc_thresh = iv$config$fc_cutoff_reg, expr_thresh = iv$config$gene_cutoff_reg, MIN_OBS = 0)
cell_type_info_renorm[[1]] <- cell_type_info_renorm[[1]][gene_list,]
saveRDS(cell_type_info_renorm, file.path(metadir, "cell_type_info_renorm.RDS"))
