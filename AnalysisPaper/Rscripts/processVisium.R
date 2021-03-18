#this script converts a DGE and coords matrix to a puck object saved as RDS
library(RCTD)
library(Matrix)
print("prepareCropper: begin")
slideseqdir <- "Data/Slideseq/Visium"
counts <- readMM(file.path(slideseqdir, 'matrix.mtx'))
barcodes <- read.table(file = file.path(slideseqdir, 'barcodes.tsv'), sep = '\t', header = F)
genes <- read.table(file = file.path(slideseqdir, 'features.tsv'), sep = '\t', header = F)
rownames(counts) = genes$V2
colnames(counts) = barcodes$V1
coords <- readr::read_csv(file = paste(slideseqdir,"BeadLocationsForR.csv",sep="/"))
colnames(coords)[2] = 'x' #renaming xcoord -> x
colnames(coords)[3] = 'y' #renaming ycoord -> y
coords = tibble::column_to_rownames(coords, var = "barcodes")
coords$barcodes <- NULL
puck = Slideseq(coords, as(as(counts,"matrix"),"dgCMatrix"))
puck <- restrict_counts(puck, rownames(puck@counts), UMI_thresh = 10, UMI_max = 200000)
puck = restrict_puck(puck, colnames(puck@counts))
plot_puck_wrapper(puck, puck@nUMI, max_val = 100000, maxUMI = 200000, minUMI=00000)
plot_puck_wrapper(puck, puck@counts['C1ql2',], max_val = 10, maxUMI = 200000) #Dentate
plot_puck_wrapper(puck, puck@counts['Trf',], max_val = 20, maxUMI = 200000) #Oligo
plot_puck_wrapper(puck, puck@counts['Nxph3',], max_val = 20, maxUMI = 200000) #Entorihinal
plot_puck_wrapper(puck, puck@counts['Fibcd1',], max_val = 20, maxUMI = 200000) #CA3
plot_puck_wrapper(puck, puck@counts['Pvrl3',], max_val = 20, maxUMI = 200000) #CA3
coordMat = puck@coords
coordMat[names(puck@nUMI), "nUMI"] = puck@counts['Trf',]
write.table(coordMat, file = file.path(slideseqdir,"cropper_input.txt"),row.names=FALSE,col.names= FALSE, sep="\t")
to_keep = read.table(file = file.path(slideseqdir,"cropper_output.txt"), sep="\t")
rownames(to_keep) = rownames(puck@coords)
puck_restr <- restrict_puck(puck, rownames(puck@coords[as.logical(to_keep[,1]),]))
plot_puck_wrapper(puck_restr, puck_restr@nUMI, max_val = 200000, maxUMI = 200000)
plot_puck_wrapper(puck_restr, puck_restr@counts['Trf',], max_val = 20, maxUMI = 200000) #Oligo
saveRDS(puck_restr, file.path(slideseqdir, "puckCropped.RDS"))
