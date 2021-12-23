library(Matrix)
library(spacexr)
library(doParallel)
library(ggplot2)

datadir <- '../RCTD/data/SpatialRNA/Puck_210605/21'
reference <- readRDS('../RCTD/data/Reference/DropVizHC/scRefSubsampled1000.RDS')
nUMI <- reference@meta.data$nUMI
names(nUMI) <- rownames(reference@meta.data)
cell_type_names <- reference@meta.data$liger_ident_coarse
names(cell_type_names) <- rownames(reference@meta.data)
reference <- Reference(reference@assays$RNA@counts,cell_type_names, nUMI)
puck <- readRDS(file.path(datadir,'puckCropped.rds'))

myRCTD <- create.RCTD(puck, reference, max_cores = 2)
myRCTD <- fitBulk(myRCTD)
myRCTD <- choose_sigma_c(myRCTD)
myRCTD <- fitPixels(myRCTD)
saveRDS(myRCTD,file.path(datadir,'myRCTD_j20.rds'))
create_RCTD_plots(myRCTD, datadir)
