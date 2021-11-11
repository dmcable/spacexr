#puck 24
library(RCTD)
library(Matrix)
counts <- readr::read_csv('data/Reference/Testes/Testis_Ref_DGE.csv')
counts = tibble::column_to_rownames(counts, var = "X1")
metaData <- readr::read_csv('data/Reference/Testes/Major_Cluster_Ident.csv')
metaData = tibble::column_to_rownames(metaData, var = "barcode")
myCounts <- as(counts,"matrix")
metaData <- as.data.frame(metaData)
metaData$X1 <- NULL
metaData$liger_ident_coarse <- factor(metaData$cluster)
metaData$nUMI <- rowSums(counts)

refdir <- 'data/Reference/Testes'
datadir <- 'data/SpatialRNA/Testes'
myCounts <- t(myCounts)
reference = Seurat::CreateSeuratObject(myCounts, meta.data = metaData)

saveRDS(reference, paste(refdir,"SCRef.RDS",sep="/"))
dref <- create_downsampled_data(reference, refdir)
reference <- dref

puck <- read.SpatialRNA('data/SpatialRNA/Testes/')
saveRDS(puck, file = 'data/SpatialRNA/Testes/puck.rds')
puck <- readRDS('data/SpatialRNA/Testes/puck.rds')

myRCTD <- create.RCTD(puck, reference, test_mode = FALSE) # need to select some more genes
myRCTD <- run.RCTD(myRCTD)
saveRDS(myRCTD,'data/SpatialRNA/Testes/myRCTD_testes.rds')
plot_weights_doublet(myRCTD@cell_type_info$info[[2]], myRCTD@spatialRNA, datadir, myRCTD@results$weights_doublet,
                     myRCTD@results$results_df)
plot_all_cell_types(myRCTD@results$results_df, myRCTD@spatialRNA@coords, myRCTD@cell_type_info$info[[2]], datadir)
plot_cond_occur(myRCTD@cell_type_info$info[[2]], datadir, norm_weights, spatialRNA)