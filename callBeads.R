library(RCTD)
library(Matrix)

iv <- init_RCTD(get_proportions = T) #initial variables
test_results = process_data(iv$puck, iv$gene_list, iv$cell_type_info, iv$proportions, trust_model = T, constrain = F)
saveRDS(test_results, file = paste(iv$resultsdir,"test_results.RDS",sep="/"))
conf_mat = test_results[[1]]; weights = test_results[[2]]; pred_labels = test_results[[3]]
puck@cell_labels = factor(names(test_results[[3]]),iv$cell_type_info[[2]])
names(iv$puck@cell_labels) = colnames(iv$puck@counts)
saveRDS(iv$puck,file.path(iv$resultsdir,"puck_labeled.RDS"))
plot_cell_types_ind(iv$puck, iv$resultsdir)
plot_cell_types(iv$puck, colnames(iv$puck@counts), iv$resultsdir)
print("callBeads: occurences of cell types:")
print(diag(conf_mat$table))
print("callBeads: end")
