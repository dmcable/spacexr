library(Matrix)
library(DEGLAM)
library(doParallel)
library(ggplot2)

datadir <- '../RCTD/data/SpatialRNA/Puck_210605/21'
myRCTD<- readRDS(file.path(datadir,'myRCTD_j20.rds'))
load(file = file.path(datadir,'exvar_thresh.RData')) # exvar, explanatory_variable
aggregate_cell_types(myRCTD, names(exvar))
cell_types <- c("Astrocyte","CA1","CA3", "Denate", "Interneuron" ,  "Microglia_Macrophages", "Oligodendrocyte" )
myRCTD@config$max_cores <- 2
myRCTD <- run.de.single(myRCTD,exvar, datadir = datadir, cell_types = cell_types) #testing gene_threshold = .001
saveRDS(myRCTD, file.path(datadir, 'myRCTDde_thresh.rds'))   
#make_all_de_plots(myRCTD, datadir)
