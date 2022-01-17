library(Matrix)
library(spacexr)
library(doParallel)
library(ggplot2)
#source('~/Documents/MIT/Research/Rafalab/Projects/slideseq/Cell Demixing/ContentStructure/RCTD/PersonalFiles/private_utils.R')
#source('~/Documents/MIT/Research/Rafalab/Projects/slideseq/Cell Demixing/ContentStructure/RCTD/de.R')
#source('~/Documents/MIT/Research/Rafalab/Projects/slideseq/Cell Demixing/ContentStructure/RCTD/R/prob_model.R')
#load in RCTD obj
id <- '08'
puck_no <- paste0('190926_', id)
datadir <- paste0('/Users/dcable/Documents/MIT/Research/Rafalab/Projects/slideseq/Cell Demixing/ContentStructure/RCTD/data/SpatialRNA/CerebellumReplicates/Puck_', '190926_11')
resultsdir <- paste0('/Users/dcable/Documents/MIT/Research/Rafalab/Projects/slideseq/Cell Demixing/ContentStructure/RCTD/data/SpatialRNA/CerebellumReplicates/Puck_', puck_no)
myRCTD<- readRDS(file.path(datadir,'myRCTD_cer_reps.rds'))
load(file.path(datadir,"regions.RData"))
#nodular = substr(nodular_08, start=1,stop=nchar(nodular_08)-3)
#anterior = substr(anterior_08, start=1,stop=nchar(anterior_08)-3)
explanatory.variable <- c(rep(0,length(nodular_08)), rep(1,length(anterior_08)))#FILL IN
names(explanatory.variable) <- c(nodular_08, anterior_08)
#Check cell types
cell_types = c('Astrocytes','Bergmann','Granule','Purkinje','Oligodendrocytes')
CSIDE:::aggregate_cell_types(myRCTD, names(explanatory.variable))
#de
puck <- readRDS(file.path(resultsdir, 'puckCropped.rds'))
myRCTD@originalSpatialRNA <- CSIDE:::coerce_old(puck)
colnames(myRCTD@originalSpatialRNA@counts) <- unlist(lapply(colnames(myRCTD@originalSpatialRNA@counts), function(x) paste0(x,'_',id)))
rownames(myRCTD@originalSpatialRNA@coords) <- unlist(lapply(rownames(myRCTD@originalSpatialRNA@coords), function(x) paste0(x,'_',id)))
names(myRCTD@originalSpatialRNA@nUMI) <- unlist(lapply(names(myRCTD@originalSpatialRNA@nUMI), function(x) paste0(x,'_',id)))
myRCTD@config$max_cores <- 4
SHUFFLE_VAR <- F
if(SHUFFLE_VAR) {
  explanatory.variable <- sample(explanatory.variable)
  names(explanatory.variable) <- c(nodular_08, anterior_08)
}
#if time:
#system.time({myRCTD <- run.de.single(myRCTD, explanatory.variable, cell_types = cell_types)})
myRCTD <- run.de.single(myRCTD, explanatory.variable, cell_types = cell_types)
if(SHUFFLE_VAR) {
  saveRDS(myRCTD,file.path(resultsdir,'myRCTDde_shuffle.rds'))
  myRCTD@de_results$res_gene_list
  table(myRCTD@de_results$gene_fits$con_mat[,'Granule'])
  Z_score <- myRCTD@de_results$gene_fits$mean_val[,'Granule'] / myRCTD@de_results$gene_fits$I_mat[,6]
  qqnorm(Z_score)
  abline(0,1)
} else {
  saveRDS(myRCTD,file.path(resultsdir,'myRCTDde.rds'))
}
#plots
cell_types_present = c('Astrocytes','Bergmann','Granule','Purkinje','MLI1','MLI2','Oligodendrocytes')
make_all_de_plots(myRCTD, resultsdir, cell_types_present = cell_types_present)
