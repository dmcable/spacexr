library(RCTD)
library(Matrix)
library(dplyr)
library(ggplot2)
iv <- init_RCTD(load_info_renorm = T) #initial variables
puck = iv$puck
resultsdir <- paste0(iv$SpatialRNAdir,"/results")
results <- gather_results(puck,iv)
norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/')

#make the plots
plot_weights(iv$cell_type_info, puck, resultsdir, norm_weights)
plot_weights_unthreshold(iv$cell_type_info, puck, resultsdir, norm_weights)
plot_weights_doublet(iv$cell_type_info, puck, resultsdir, results$weights_doublet, results$results_df)
plot_cond_occur(iv$cell_type_info, resultsdir, norm_weights)
plot_all_cell_types(results$results_df, puck@coords, iv$cell_type_info, resultsdir)

#doublets
doublets <- results$results_df[results$results_df$spot_class == "doublet_certain",]
plot_doublets(puck, doublets, resultsdir, iv$cell_type_info)
plot_doublets_type(puck, doublets, resultsdir, iv$cell_type_info)
doub_occur <- table(doublets$second_type, doublets$first_type)
plot_doub_occur_stack(doub_occur, resultsdir, iv)

#get decomposed puck
puck_d <- get_decomposed_data(results$results_df, iv$gene_list, puck, results$weights_doublet, iv$cell_type_info)
