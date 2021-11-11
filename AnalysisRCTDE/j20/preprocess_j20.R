library(Matrix)
library(DEGLAM)
library(doParallel)
library(ggplot2)
library(geometry)
source('~/Documents/MIT/Research/Rafalab/Projects/slideseq/Cell Demixing/ContentStructure/DEGLAM/analysis/personal_utils.R')
#datadir <- '../RCTD/data/SpatialRNA/Puck_210605/24'
datadir <- '../RCTD/data/SpatialRNA/Puck_210605/James_Results'
puck_no <- 'Puck_210605_24'
puck <- read.Slideseq.mtx(datadir, puck_no)
puck <- restrict_counts(puck, rownames(puck@counts), UMI_thresh = 100)
puck = restrict_puck(puck, colnames(puck@counts))

datadir <- '../RCTD/data/SpatialRNA/Puck_210605/24'
puck_old <- readRDS(file.path(datadir,'puck_old.rds'))
puckCropped_old <- readRDS(file.path(datadir,'puckCropped_old.rds'))
# crops puck according to convex hull of puckCropped_old
puckCropped <- crop.puck.convex(puck, puck_old, puckCropped_old, DIST.THRESH = 30, DIST.THRESH.BAD = 30)
dim(puckCropped@counts)
dim(puckCropped@counts)[2] / dim(puckCropped_old@counts)[2]
saveRDS(puckCropped,file.path(datadir,'puckCropped.rds'))
saveRDS(puck,file.path(datadir,'puck.rds'))
puck <- readRDS(file.path(datadir,'puck.rds'))
plot_puck_continuous(puck,rownames(puck@coords), puck@nUMI, ylimit = c(0,1500))
plot_puck_continuous(puck,rownames(puck@coords), puck@counts['C1ql2',], ylimit = c(0,1)) #Dentate
plot_puck_continuous(puck,rownames(puck@coords), puck@counts['Trf',], ylimit = c(0,1)) #Oligo

library(gatepoints)
X <- puck@coords
pal = colorRampPalette(c("blue", "red"))
dforder = findInterval(puck@nUMI, sort(puck@nUMI))
plot(X, pch=20, col=pal(nrow(X))[dforder], cex = 0.5)
selectedPoints <- fhs(X, mark = TRUE,  pch=20, col=pal(nrow(X))[dforder], cex = 0.5)
arr <- as.integer(names(puck@nUMI) %in% selectedPoints); names(arr) <- names(puck@nUMI)
RCTD::plot_puck_continuous(puck, names(puck@nUMI), arr, ylimit = c(-1,2))
puckCropped = restrict_puck(puck,names(which(arr == 1)))
plot_puck_continuous(puckCropped,rownames(puckCropped@coords), puckCropped@counts['Trf',], ylimit = c(0,1)) #Oligo
saveRDS(puckCropped,file.path(datadir,'puckCropped.rds'))

reference <- readRDS('../RCTD/data/Reference/DropVizHC/scRefSubsampled1000.RDS')
nUMI <- reference@meta.data$nUMI
names(nUMI) <- rownames(reference@meta.data)
cell_type_names <- reference@meta.data$liger_ident_coarse
names(cell_type_names) <- rownames(reference@meta.data)
reference <- Reference(reference@assays$RNA@counts,cell_type_names, nUMI)
