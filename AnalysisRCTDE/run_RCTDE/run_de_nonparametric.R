myRCTD@config$max_cores <- 4
myRCTDde <- run.de.nonparam(myRCTD, df = 15, gene_threshold = .001) 
myRCTDde <- run.de.nonparam(myRCTD, df = 15) #, gene_threshold = .001) 
myRCTDde@de_results$gene_fits$error_mat
myRCTDde@de_results$gene_fits$con_mat
myRCTDde@de_results$gene_fits$all_vals[gene,,]
res2$res$alpha2
myRCTDde <- add_res_genes(myRCTDde, datadir = './', plot_genes = F, param_position = 2:15, 
                          fdr = 0.1, log_fc_thresh = 0.01)
saveRDS(myRCTDde,file.path(resultsdir,'myRCTDde.rds'))
