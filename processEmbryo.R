library(readr)
library(Matrix)
library(RCTD)
meta_data <- read_csv('Data/Reference/Embryo/GSE119945_cell_annotate.csv')
meta_data <- as.data.frame(meta_data)
restr_ind <- meta_data$day %in% c(11.5,12.5) & !meta_data$detected_doublet & !meta_data$potential_doublet_cluster & !is.na(meta_data$Total_mRNAs)
restr_ind[is.na(restr_ind)] <- F
restr <- meta_data[restr_ind,]


restr <- restr[index_keep,]

library('monocle')
dge <- readRDS('Data/Reference/Embryo/cds_cleaned.RDS')
index_keep = c(); n_samples = 10000
for(i in 1:38) {
  new_index = which(dge$Cluster == i & dge$day %in% c(11.5,12.5))
  new_samples = min(n_samples, length(new_index))
  index_keep = c(index_keep, sample(new_index,new_samples,replace=FALSE))
}
counts <- readRDS(paste0('Data/Reference/Embryo/counts',ind,'.RDS'))
saveRDS(paste0('Data/Reference/Embryo/counts',ind,'.RDS'))
for(ind in 1:7) {
  start = 1 + 40000*ind
  end = min(length(index_keep), 40000*(ind+1))
  #saveRDS(dge@assayData$exprs[,index_keep[start:end]],paste0('Data/Reference/Embryo/counts',ind,'.RDS'))
  counts_new <- readRDS(paste0('Data/Reference/Embryo/counts',ind,'.RDS'))
  counts <- cbind(counts,counts_new)
}
saveRDS(dge@assayData$exprs[,index_keep],'Data/Reference/Embryo/counts.RDS')
nUMI <- colSums(counts)
cluster = dge$Cluster[index_keep]
saveRDS(cluster, paste0('Data/Reference/Embryo/cluster.RDS'))

cell_type_names = c('Endo', 'Macrophage', 'T_Cell', 'B_Cell', 'Dendtritic', 'Cholangiocyte', 'Hepatocyte', 'Dividing', 'Plasma_B','HSC')

meta_data = as.data.frame(cluster)
meta_data$nUMI = nUMI
rownames(meta_data) = colnames(counts)
refdir = "Data/Reference/Embryo"
cell_dict_file <- file.path(refdir,"cell_type_dict.csv")
true_type_names <- remap_celltypes(cell_dict_file, meta_data$cluster)
meta_data$liger_ident_coarse = droplevels(true_type_names)
reference = Seurat::CreateSeuratObject(counts, meta.data = meta_data)
saveRDS(reference, paste(refdir,"SCRef.RDS",sep="/"))


#Get cell type info
config_data <- config::get(file = "conf/dataset.yml", use_parent = FALSE)
print(paste("init_RCRD: using config mode:",config_data$config_mode))
config <- config::get(file = paste0("conf/",config_data$config_mode,".yml"), use_parent = FALSE)
slideseqdir <- file.path("Data/Slideseq",config_data$slideseqfolder)
resultsdir = file.path(slideseqdir,"results")
if(!dir.exists(resultsdir))
  dir.create(resultsdir)
bulkdir <- paste(slideseqdir,"results/Bulk",sep="/")
if(!dir.exists(bulkdir))
  dir.create(bulkdir)
refdir <- file.path("Data/Reference",config_data$reffolder)
reference <- readRDS(paste(refdir,config_data$reffile,sep="/"))
print(paste("init_RCTD: number of cells in reference:", dim(reference@assays$RNA@counts)[2]))
print(paste("init_RCTD: number of genes in reference:", dim(reference@assays$RNA@counts)[1]))
cell_counts = table(reference@meta.data$liger_ident_coarse)
print(cell_counts)
CELL_MIN = 25 # need at least this for each cell type
if(min(cell_counts) < CELL_MIN)
  stop(paste0("init_RCTD error: need a minimum of ",CELL_MIN, " cells for each cell type in the reference"))
cell_type_info <- get_cell_type_info(reference@assays$RNA@counts, reference@meta.data$liger_ident_coarse, reference@meta.data$nUMI)
saveRDS(cell_type_info,file.path(refdir,'MetaData/cell_type_info.RDS'))
cell_type_info <- readRDS(file.path(refdir,'MetaData/cell_type_info.RDS'))
puck = readRDS(file.path(slideseqdir, config_data$puckrds))
puck = restrict_counts(puck, rownames(puck@counts), UMI_thresh = 100, UMI_max = 200000)
gene_ann <- read.csv(file.path(refdir,'GSE119945_gene_annotate.csv'))
rownames(gene_ann) = gene_ann$gene_id

dup <- which(table(gene_ann$gene_short_name) > 1)
for (gene in names(dup)) {
  genes <- rownames(gene_ann)[gene_ann$gene_short_name == gene]
  cell_type_info[[1]][genes[1],] <- colSums(cell_type_info[[1]][genes,])
  cell_type_info[[1]] <- cell_type_info[[1]][!(rownames(cell_type_info[[1]]) %in% genes[2:length(genes)]),]
}

rownames(cell_type_info[[1]]) = gene_ann[rownames(cell_type_info[[1]]),]$gene_short_name
gene_list = get_de_genes(cell_type_info, puck, fc_thresh = config$fc_cutoff, expr_thresh = config$gene_cutoff, MIN_OBS = 3)
gene_list = get_de_genes(cell_type_info, puck, fc_thresh = 0.01, expr_thresh = .000125, MIN_OBS = 3)
gene_list_reg = get_de_genes(cell_type_info, puck, fc_thresh = 0.01, expr_thresh = .00015, MIN_OBS = 3)


