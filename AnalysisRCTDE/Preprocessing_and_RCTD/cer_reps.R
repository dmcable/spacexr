library(Matrix)
library(DEGLAM)
library(doParallel)
datadir <- '../RCTD/data/SpatialRNA/Puck_190926_11/'
puck <- read.SpatialRNA(datadir)
puck <- restrict_counts(puck, rownames(puck@counts), UMI_thresh = 100)
puck = restrict_puck(puck, colnames(puck@counts))
saveRDS(puck,file.path(datadir,'puck.rds'))
plot_puck_continuous(puck,rownames(puck@coords), puck@nUMI, ylimit = c(0,1500))
plot_puck_continuous(puck,rownames(puck@coords), puck@counts['Pcp4',], ylimit = c(0,10))
library(gatepoints)
X <- puck@coords
pal = colorRampPalette(c("blue", "red"))
dforder = findInterval(puck@nUMI, sort(puck@nUMI))
plot(X, pch=20, col=pal(nrow(X))[dforder], cex = 0.5)
selectedPoints <- fhs(X, mark = TRUE,  pch=20, col=pal(nrow(X))[dforder], cex = 0.5)
arr <- as.integer(names(puck@nUMI) %in% selectedPoints); names(arr) <- names(puck@nUMI)
RCTD::plot_puck_continuous(puck, names(puck@nUMI), arr, ylimit = c(-1,2))
puckCropped = restrict_puck(puck,names(which(arr == 1)))
saveRDS(puckCropped,file.path(datadir,'puckCropped.rds'))
### REPEAT THE ABOVE FOR PUCKS 08, 09, and 11 to create puckCropped
#### end generic puck cropper

### start merge puck
puck <- readRDS('data/SpatialRNA/NewCerPuck_190926_08/puckCropped.RDS')
puck <- restrict_puck(puck, colnames(puck@counts))
newnames = paste(rownames(puck@coords), "08",sep="_")
puck <- rename_puck(puck,newnames)
rename_puck <- function(puck, newnames) {
  rownames(puck@coords) = newnames
  colnames(puck@counts) = newnames
  names(puck@nUMI) = newnames
  return(puck)
}
newnames = paste(rownames(puckCropped@coords), "11",sep="_")
puckCropped <- rename_puck(puckCropped,newnames)
reference <- readRDS('data/Reference/10xCer/scRefSubsampled1000.RDS')
merge_puck <- function(puck1, puck2) {
  m_puck <- puck1
  gene_list <- intersect(rownames(m_puck@counts), rownames(puck2@counts))
  name_overlap <- intersect(rownames(m_puck@coords), rownames(puck2@coords))
  if(length(name_overlap) > 0)
    stop('overlapping pixel barcodes')
  m_puck@counts <- cbind(m_puck@counts[gene_list,], puck2@counts[gene_list,])
  m_puck@coords <- rbind(m_puck@coords[,1:2], puck2@coords[,1:2])
  m_puck@nUMI <- c(m_puck@nUMI, puck2@nUMI)
  return(m_puck)
}
m_puck <- merge_puck(puck, puckCropped)
puck3 <- readRDS('data/SpatialRNA/Puck_190926_09/puckCropped.rds')
newnames = paste(rownames(puck3@coords), "09",sep="_")
puck3 <- rename_puck(puck3,newnames)
m_puck <- merge_puck(m_puck, puck3)
### END MERGE PUCK

### RUN RCTD
myRCTD <- create.RCTD(m_puck, reference, test_mode = FALSE)
myRCTD <- fitBulk(myRCTD)
myRCTD <- choose_sigma_c(myRCTD)
myRCTD <- fitPixels(myRCTD, doublet_mode = 'doublet')
saveRDS(myRCTD,'myRCTD_cer_reps.rds')
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
M = 4500
puck <- myRCTD@spatialRNA
puck@coords[substrRight(rownames(puck@coords),2) == "09",'x'] = puck@coords[substrRight(rownames(puck@coords),2) == "09",'x'] + M
puck@coords[substrRight(rownames(puck@coords),2) == "11",'y'] = puck@coords[substrRight(rownames(puck@coords),2) == "11",'y'] + M
myRCTD@spatialRNA <- puck
### END RUN RCTD

myRCTD <- readRDS(file.path(datadir,'myRCTD_cer_reps.rds'))
plot_weights_doublet(myRCTD@cell_type_info$renorm[[2]], puck, datadir, myRCTD@results$weights_doublet, myRCTD@results$results_df)
### DEFINE REGIONS FOR DE
plot_puck_continuous(puck,rownames(puck@coords), puck@nUMI, ylimit = c(0,1500))
X <- puck@coords
pal = colorRampPalette(c("blue", "red"))
dforder = findInterval(puck@nUMI, sort(puck@nUMI))
plot(X, pch=20, col=pal(nrow(X))[dforder], cex = 0.5)
selectedPoints <- fhs(X, mark = TRUE,  pch=20, col=pal(nrow(X))[dforder], cex = 0.5)
arr <- as.integer(names(puck@nUMI) %in% selectedPoints); names(arr) <- names(puck@nUMI)
RCTD::plot_puck_continuous(puck, names(puck@nUMI), arr, ylimit = c(-1,2))
nodular_08 <- names(which(arr == 1))
anterior_08 <- names(which(arr == 1))
nodular_09 <- names(which(arr == 1))
anterior_09 <- names(which(arr == 1))
nodular_11 <- names(which(arr == 1))
anterior_11 <- names(which(arr == 1))
## REPEAT THE ABOVE TO DEFINE EACH OF THE 6 REGIONS ACROSS 3 PUCKS
save(nodular_08, nodular_09, nodular_11, anterior_08, anterior_09, anterior_11, file = file.path(datadir,'regions.RData'))