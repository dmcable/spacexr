#prepares for NMFreg
library(readr)
config_data <- config::get(file = "conf/dataset.yml", use_parent = FALSE)
config <- config::get(file = paste0("conf/",config_data$config_mode,".yml"), use_parent = FALSE)
slideseqdir <- file.path("Data/Slideseq",config_data$slideseqfolder)
resultsdir = file.path(slideseqdir,"results")
outdir = 'NMFreg/DataCer'
refdir <- file.path("Data/Reference",config_data$reffolder)
reference <- readRDS(paste(refdir,config_data$reffile,sep="/"))
puck = readRDS(file.path(slideseqdir, config_data$puckrds))
counts <- as.data.frame(t(puck@counts))
counts2 <- counts
counts2 <- cbind(rownames(counts2), counts2)
colnames(counts2)[1] = "barcode"
write_csv(counts2, file.path(outdir,"counts.csv"))
coords <- as.data.frame(puck@coords)
coords <- cbind(rownames(coords), coords)
colnames(coords)[1] = "barcode"
coords[,"nUMI"] <- NULL
write_csv(coords, file.path(outdir,"coords.csv"))
refcounts <- as.data.frame(t(reference@assays$RNA@counts))
counts2 <- refcounts
counts2 <- cbind(rownames(counts2), counts2)
colnames(counts2)[1] = "barcode"
write_csv(counts2, file.path(outdir,"dge_hvgs.csv"))
#remember to delete the "barcode" by hand
cluster <- as.data.frame(as.integer(reference@meta.data$liger_ident_coarse))
cluster <- cbind(rownames(reference@meta.data), cluster)
colnames(cluster) = c('barcode', 'cluster')
write_csv(cluster, file.path(outdir,"cell_cluster_outcomes.csv"))
