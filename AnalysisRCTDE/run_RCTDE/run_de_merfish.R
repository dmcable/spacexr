library(Matrix)
library(doParallel)
library(ggplot2)
library(dplyr)
library(RCTD)
library(GSA)
library(xlsx)
library(readr)
# Data is obtained from the following article.
# https://science.sciencemag.org/content/362/6416/eaau5324/tab-figures-data
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248 merfish data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113576 reference data
# normalization data was obtained by emailing the authors of the article
# Directory setup
pwd = getwd()
datadir <- paste0(pwd,'/data/moffitt','/')
resultsdir <- paste0(pwd,'/results/ResultsMerfish','/')

# Load in spatialRNA data and create spatialRNA object
# colnames are Cell_ID, Animal_ID, Animal_sex, Behavior, Bregma, Centroid_X, Centroid_Y, Neuron_cluster_ID, All_the_other_genes
# 161 genes in the columns here
merfish_data = as.data.frame(read_csv(paste0(datadir,'Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv')))
merfish_barcodes = merfish_data[,1]
normalization_data = as.data.frame(read_csv(paste0(datadir,'volume_and_batch_correction.csv')))

# Filter coords, counts, and normalization_data so they only contain data where
# animal ID = 1 and bregma = .01.
# This represents data from one slice of tissue
# nUMI is calculated later from counts
normalization_data <- normalization_data[merfish_data$Animal_ID==1 & merfish_data$Bregma==.01,]
merfish_data <- merfish_data[merfish_data$Animal_ID==1 & merfish_data$Bregma==.01,]
rownames(merfish_data) <- merfish_data[,1]; rownames(normalization_data) <- merfish_data[,1]
merfish_data[,1] <- NULL
normalization_data = select(normalization_data, abs_volume, batch_correction)

# Create SpatialRNA object using coords, counts, and nUMI
coords = select(merfish_data,Centroid_X,Centroid_Y)
counts = select(merfish_data, -Animal_ID, -Animal_sex, -Behavior, -Bregma, -Centroid_X, -Centroid_Y, -Neuron_cluster_ID, -Cell_class)
counts = counts[,1:140] # Throw out genes at and after Adcyap1

# Reverse the normalization done on expression values to get their raw integer counts.
# Counts is barcodes x gene, original counts were  (counts / volume) * batch_correction value.
counts = sweep(counts, 1,normalization_data[,2], '/') # reverse batch correction
counts = sweep(counts,1, normalization_data[,1],'*') # reverse division by cell volume

# Verify most of the values in counts are now basically integers
close_integer = function(val){
	if(val==0){return(FALSE)} # don't care really about the 0's
	eps = .05
	remainder = val %% 1
	if (remainder + eps > 1 | remainder - eps < 0)
		return(TRUE)
	print(val)
	return(FALSE)
}
non_zero_count = sum(counts > 0)
close_integers_matrix = apply(counts, MARGIN=c(1,2) ,FUN=close_integer)
total_close_integers = sum(close_integers_matrix)
close_integer_ratio = total_close_integers / non_zero_count # should be close to 1

# Make the counts values integers
counts = as.matrix(round(counts))
mode(counts) = "integer"
# counts = data.frame(counts) # If counts needs to be in dataframe form.

counts = t(counts)
counts <- Matrix(counts)
nUMI <- colSums(counts)

# Clear massive variables not going to be used again. Anything >1GB
# rm(merfish_data)

puck = SpatialRNA(coords,counts,nUMI)

# Examine spatialRNA object (optional)
print(dim(puck@counts)) # 140 genes x 6111 cells
hist(log(puck@nUMI,2))
puck_barcodes <- colnames(puck@counts)
plot_puck_continuous(puck, puck_barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))),
                     title ='plot of nUMI')
saveRDS(puck, file.path(resultsdir,'puck.rds'))
# Load in reference data and then create Reference Object
# 27998x31299 genes x barcodes
reference_matrix_data = readMM(paste0(datadir,'GSE113576_matrix.mtx'))
reference_barcodes_data = read.table(file=paste0(datadir,'GSE113576_barcodes.tsv'))
reference_genes_data = read.table(file=paste0(datadir,'GSE113576_genes.tsv'))
cell_type_data = xlsx::read.xlsx(file=paste0(datadir,'aau5324_Moffitt_Table-S1.xlsx'),sheetIndex =1 ,header=TRUE)
colnames(cell_type_data) = c("barcodes","sex","replicate_number","cell_type","non_neuronal_cluster","neuronal_cluster")
cell_type_data = cell_type_data[-1,]; # Set actual column names since file is loaded oddly

# Create Reference Object using counts, and cell_types
counts = Matrix(reference_matrix_data)
mode(counts) = "integer" # convert matrix data from numeric to integers
reference_barcodes = reference_barcodes_data[,1]
colnames(counts) = reference_barcodes
reference_genes = reference_genes_data[,2]
rownames(counts) = reference_genes

# Create nUMI before filtering out any data from counts
nUMI = colSums(counts)
names(nUMI) <- colnames(counts)

# There are duplicated genes in counts. Remove them.
duplicate_rows_boolean = duplicated(rownames(counts))
duplicate_row_indices = which(duplicate_rows_boolean)
duplicate_gene_names = unique(rownames(counts)[duplicate_row_indices])
# Can't convert to dataframe and filter out rows because conversion creates unique rownames.
counts = counts[!(rownames(counts) %in% duplicate_gene_names), ]

cell_types = cell_type_data$cell_type
names(cell_types) = cell_type_data$barcodes
cell_types <- as.factor(cell_types)

# Clear massive variables no longer being used, anything >1GB
# rm(reference_matrix_data)

# 27877 genes x 31299 cells
reference <- Reference(counts, cell_types, nUMI)
saveRDS(reference, file.path(resultsdir,'reference.rds'))
# Create and run RCTD
cell_type_names <- levels(cell_types)[2:13] # remove ambiguous, unstable
myRCTD <- create.RCTD(puck, reference, max_cores = 4, gene_cutoff = 0.000, fc_cutoff = 0, gene_cutoff_reg = 0,
                      fc_cutoff_reg = 0, cell_type_names = cell_type_names)
saveRDS(myRCTD, file.path(resultsdir,'myRCTD_pre.rds'))
# Save/Load RCTD object
# saveRDS(myRCTD,file.path(resultsdir,'preRCTD.rds'))
myRCTD = readRDS(paste0(resultsdir,'myRCTD_pre.rds')) # 53 DE genes found from 140

myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
# Save/Load postRCTD object
create_RCTD_plots(myRCTD, resultsdir)
saveRDS(myRCTD,file.path(resultsdir,'postRCTD.rds'))
myRCTD = readRDS(paste0(resultsdir,'postRCTD.rds'))

# Create explanatory.variable
cell_types = myRCTD@cell_type_info$info[[2]]
barcodes <- colnames(myRCTD@spatialRNA@counts)
explanatory.variable = c(rep(0,length(barcodes)))
names(explanatory.variable) = barcodes
x_coords = myRCTD@spatialRNA@coords[,1]
midline_x = median(x_coords)
for(i in 1:length(barcodes)) {
  if(i%%50==0) { print(i) } # Track progress
  barcode = barcodes[i]
  curr_x = myRCTD@spatialRNA@coords[barcode,1]
  distance = abs(curr_x - midline_x)
  explanatory.variable[i]=distance
}
print(max(abs(explanatory.variable))) #896.1267
explanatory.variable = explanatory.variable / max(explanatory.variable)
hist(explanatory.variable) # Observe the distribution of explanatory.variable
plot_puck_continuous(myRCTD@spatialRNA, colnames(myRCTD@spatialRNA@counts), explanatory.variable, # ylimit = c(0,round(quantile(myRCTD@spatialRNA@nUMI,0.9))),
                     title ='plot of explanatory variable')


# run DE
myRCTD <- run.de.single(myRCTD, explanatory.variable, datadir = resultsdir, cell_type_threshold = 75) # outputs plots into pwd() + "/de_plots/"
make_all_de_plots(myRCTD, resultsdir)
# Save/Load and plot resulting RCTDde
saveRDS(myRCTD,file.path(resultsdir,'myRCTDde.rds'))


#### RUN QUADRATIC

### try quadratic
X2 <- myRCTD@internal_vars_de$X2
X2 <- cbind(X2, X2[,2]^2)
myRCTDQ <- find_de_genes(myRCTD, myRCTD@internal_vars_de$X1, X2, myRCTD@internal_vars_de$all_barc, cell_types, resultsdirQ,
                         cell_types_present = cell_types_present)
gene_fitsQ <- myRCTDQ@de_results$gene_fits
resultsdirQ <- paste0(pwd,'/results/ResultsMerfish/Quadratic/')
make_all_de_plots(myRCTDQ, resultsdirQ)
#cell_type <- 'Excitatory'
#ct_ind <- which(cell_types == cell_type)
#gene_list <- get_gene_list_type_wrapper(myRCTDQ, cell_type, cell_types_present)
#param_position <- 3
#gene <- 'Galr1'
#gene <- 'Ebf3'
#mu <- gene_fitsQ$all_vals[gene,param_position,ct_ind]
#sd <- gene_fitsQ$I_mat[gene,3*(ct_ind-1) + param_position]
#mu <- gene_fitsQ$all_vals[gene_list,param_position,ct_ind]
#sd <- gene_fitsQ$I_mat[gene_list,3*(ct_ind-1) + param_position]
#Z <- mu / sd
#hist(Z)
myRCTDQ <- add_res_genes(myRCTDQ, datadir = resultsdirQ, plot_genes = T, param_position = param_position)
myRCTDQ@de_results$res_gene_list$Excitatory
saveRDS(myRCTDQ,file.path(resultsdir,'myRCTDQ.rds'))
myRCTDQ <- readRDS(file.path(resultsdir,'myRCTDQ.rds'))


