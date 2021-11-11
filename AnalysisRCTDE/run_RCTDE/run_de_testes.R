library(DEGLAM)
library(Matrix)
library(ggplot2)
library(devtools)
source('~/Documents/MIT/Research/Rafalab/Projects/slideseq/Cell Demixing/ContentStructure/DEGLAM/analysis/helper_functions/testes_helper.R')
load_all()
datadir <- '../RCTD/data/SpatialRNA/Testes'
puck <- readRDS(file.path(datadir,'puck.rds'))
### END PRELUDE
myRCTD <- readRDS(file.path(datadir,'myRCTD_testes.rds'))
myRCTD@originalSpatialRNA <- puck

codes_1 <- c()
stadir <- file.path(datadir,'Stage/I-III')
for(file in list.files(stadir)) {
  labels <- read.csv(file.path(stadir,file))
  tubule_codes <- as.character(labels$barcode)
  codes_1 <- c(codes_1,tubule_codes)
  print(length(codes_1))
}

codes_2 <- c()
stadir <- file.path(datadir,'Stage/IV-VI')
for(file in list.files(stadir)) {
  labels <- read.csv(file.path(stadir,file))
  tubule_codes <- as.character(labels$barcode)
  codes_2 <- c(codes_2,tubule_codes)
  print(length(codes_2))
}

codes_3 <- c()
stadir <- file.path(datadir,'Stage/VII-VIII')
for(file in list.files(stadir)) {
  labels <- read.csv(file.path(stadir,file))
  tubule_codes <- as.character(labels$barcode)
  codes_3 <- c(codes_3,tubule_codes)
  print(length(codes_3))
}

codes_4 <- c()
stadir <- file.path(datadir,'Stage/IX-XII')
for(file in list.files(stadir)) {
  labels <- read.csv(file.path(stadir,file))
  tubule_codes <- as.character(labels$barcode)
  codes_4 <- c(codes_4,tubule_codes)
  print(length(codes_4))
}

barcodes_1 <- intersect(codes_1,names(myRCTD@spatialRNA@nUMI))
barcodes_2 <- intersect(codes_2,names(myRCTD@spatialRNA@nUMI))
barcodes_3 <- intersect(codes_3,names(myRCTD@spatialRNA@nUMI))
barcodes_4 <- intersect(codes_4,names(myRCTD@spatialRNA@nUMI))
region_list <- list(barcodes_1, barcodes_2, barcodes_3, barcodes_4)
saveRDS(region_list, file.path(datadir, 'region_list.rds'))
region_list <- readRDS(file.path(datadir, 'region_list.rds'))
all_barc <- intersect(Reduce(union,region_list), colnames(myRCTD@spatialRNA@counts))
aggregate_cell_types(myRCTD,all_barc,doublet_mode = F)
myRCTD@config$max_cores <- 4
cur_cell_types = c('1','2','4','5','6')
#region_list_sm <- list(sample(barcodes_1,100), sample(barcodes_2,100), sample(barcodes_3,100), sample(barcodes_4,100))
#myRCTDde <- run.de.regions(myRCTD, region_list_sm, datadir = datadir, gene_threshold = 0.01, cell_types = cur_cell_types,
#                           doublet_mode = F, cell_type_threshold = 5)
myRCTDde <- run.de.regions(myRCTD, region_list, datadir = datadir, cell_types = cur_cell_types,
                           doublet_mode = F)
myRCTDde@internal_vars_de$delta <- 0; myRCTDde@internal_vars_de$test_mode <- 'multi'
myRCTDde@reference <- gen_small_reference()
myRCTDde <- add_res_genes(myRCTDde, datadir = datadir, plot_genes = F, p_thresh = 1)#, p_thresh = 1)
saveRDS(myRCTDde, file.path(datadir,'myRCTDde_updated2.rds'))
myRCTDde <- readRDS(file.path(datadir,'myRCTDde_updated2.rds'))
make_all_de_plots(myRCTDde, datadir)
