#script to call the beads on a puck as single cell types
library(config) #
library(Seurat) #
library(purrr) #
source("slideseq.R")
source("utils.R")
source("processRef.R")
source("slideseq_sim.R")
source("topicmodel.R")
config <- config::get()
refdir <- file.path("Data/Reference",config$reffolder)
dir.create(file.path(refdir, "results"), showWarnings = FALSE) #folder to save results
slideseqdir <- file.path("Data/Slideseq",config$slideseqfolder)
bulkdir <- paste(slideseqdir,"results/Bulk",sep="/")
if(!dir.exists(bulkdir))
  dir.create(bulkdir)
puck = read.slideseq(slideseqdir, count_file = config$puckfile)
reference <- readRDS(paste(refdir,config$reffile,sep="/"))
cell_type_info <- get_cell_type_info(reference@raw.data, reference@meta.data$liger_ident_coarse)
cell_type_means = cell_type_info[[1]]; cell_type_names = cell_type_info[[2]]
n_cell_types = cell_type_info[[3]]

proportions <- read.csv(file.path(bulkdir,"weights.csv"))$Weight
names(proportions) = cell_type_names
#now we switch to a different gene list
gene_list = get_gene_list(cell_type_means, puck, cutoff_val = config$gene_cutoff_reg)
print(paste("callBeads: number of genes used for regression:", length(gene_list)))
puck = restrict_counts(puck, gene_list, UMI_thresh = config$UMI_min)
puck = restrict_puck(puck, colnames(puck@counts))

test_results = process_data(puck, gene_list, cell_type_info, proportions, trust_model = T, constrain = T)
resultsdir = file.path(slideseqdir,"results")
saveRDS(test_results, file = paste(resultsdir,"test_results.RDS",sep="/"))
conf_mat = test_results[[1]]; weights = test_results[[2]]; pred_labels = test_results[[3]]
puck@cell_labels = factor(names(test_results[[3]]),cell_type_names)
names(puck@cell_labels) = colnames(puck@counts)
saveRDS(puck,file.path(resultsdir,"puck_labeled.RDS"))
plot_cell_types_ind(puck, resultsdir)
plot_cell_types(puck, colnames(puck@counts), resultsdir)
