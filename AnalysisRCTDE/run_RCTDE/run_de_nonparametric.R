library(Matrix)
library(spacexr)
library(doParallel)
library(ggplot2)
library(devtools)
library(mgcv)
load_all()

pwd = getwd()
datadir <- paste0(pwd,'/data/tumor','/')
resultsdir <- paste0(pwd,'/results/ResultsTumorNonparam','/')

# Load RCTD which contains cropped puck, also load full puck
myRCTD = readRDS(paste0(datadir,'RCTD_merged.rds')) # 21902 genes time start 2:30pm 640 genes done at 2:53, estimated finish time 13 hours from start...
cropped_puck = myRCTD@spatialRNA
myRCTD@originalSpatialRNA = cropped_puck

myRCTD@config$max_cores <- 4
myRCTDde <- run.de.nonparam(myRCTD, df = 15) #, gene_threshold = .001)
myRCTDde@de_results$gene_fits$error_mat
myRCTDde@de_results$gene_fits$con_mat
myRCTDde@de_results$gene_fits$all_vals[gene,,]
myRCTDde <- add_res_genes(myRCTDde, datadir = './', plot_genes = F, param_position = 2:15,
                          fdr = 0.1, log_fc_thresh = 0.01)
saveRDS(myRCTDde,file.path(resultsdir,'myRCTDde.rds'))
