#this script converts a DGE and coords matrix to a puck object saved as RDS
library(RCTD)
library(Matrix)
print("prepareCropper: begin")
config_data <- config::get(file = "conf/dataset.yml", use_parent = FALSE)
slideseqdir <- file.path("Data/Slideseq",config_data$slideseqfolder)
puck = read.slideseq(slideseqdir, count_file = config_data$countfile)
puck <- restrict_counts(puck, rownames(puck@counts), UMI_thresh = 50)
puck = restrict_puck(puck, colnames(puck@counts))
coordMat = puck@coords
coordMat[names(puck@nUMI), "nUMI"] = puck@nUMI
plot_puck_wrapper(puck, puck@counts['Trf',], max_val = 3, maxUMI = 200000) #Oligo
plot_puck_wrapper(puck, puck@nUMI, max_val = 2000)
write.table(coordMat, file = file.path(slideseqdir,"cropper_input.txt"),row.names=FALSE,col.names= FALSE, sep="\t")
to_keep = read.table(file = file.path(slideseqdir,"cropper_output.txt"), sep="\t")
rownames(to_keep) = rownames(puck@coords)
puck_restr <- restrict_puck(puck, rownames(puck@coords[as.logical(to_keep[,1]),]))
plot_puck_wrapper(puck_restr, puck_restr@nUMI, max_val = 2000)
plot_puck_wrapper(puck_restr, puck_restr@counts['Trf',], max_val = 3, maxUMI = 200000) #Oligo
saveRDS(puck_restr, file.path(slideseqdir, "puckCropped.RDS"))
