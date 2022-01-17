library(Matrix)
library(spacexr)
library(doParallel)
library(ggplot2)
library(xlsx)
library(dplyr)
library(fields)
library(stringr)
library(GSA)

# Load in spatialRNA data and Reference data
pwd = getwd()
datadir <- paste0(pwd,'/data/tumor','/')
resultsdir <- paste0(pwd,'/results/ResultsTumor','/')

# Load RCTD which contains cropped puck, also load full puck
myRCTD = readRDS(paste0(datadir,'RCTD_merged.rds')) # 21902 genes time start 2:30pm 640 genes done at 2:53, estimated finish time 13 hours from start...
cropped_puck = myRCTD@spatialRNA
myRCTD@originalSpatialRNA = cropped_puck
full_puck = readRDS(paste0(datadir,'puck.rds'))

## Examine SpatialRNA object (optional)
print(dim(cropped_puck@counts))
hist(log(cropped_puck@nUMI,2))
print(head(cropped_puck@coords))
barcodes <- colnames(cropped_puck@counts)
plot_puck_continuous(cropped_puck, barcodes, cropped_puck@nUMI, ylimit = c(0,round(quantile(cropped_puck@nUMI,0.9))),
                     title ='plot of nUMI cropped_puck')
plot_puck_continuous(full_puck, colnames(full_puck@counts) , full_puck@nUMI, ylimit = c(0,round(quantile(full_puck@nUMI,0.9))),
                     title ='plot of nUMI full_puck') + coord_fixed()
types <- myRCTD@results$results_df[,'first_type']
names(types) <- rownames(myRCTD@results$results_df)
plot_class(full_puck, colnames(full_puck@counts), types) + coord_fixed()
plot_class(full_puck, names(types[types %in% c('tumor II', 'CAF')]), types) + coord_fixed()
# Create explanatory.variable
cell_types = myRCTD@cell_type_info$info[[2]]
target_type = cell_types[6] # "monocyte/DC"
doublet_df = myRCTD@results$results_df
weights_doublet=myRCTD@results$weights_doublet
coords = cropped_puck@coords
emt_list <- c('Col7a1', 'Dpysl3','Inhba','Lrp1','Col6a2','Pmepa1','Mgp','Timp3','Ecm1')
Y <- colSums(cropped_puck@counts[emt_list,])
cancer_sing <- rownames(myRCTD@results$results_df)[which(myRCTD@results$results_df$spot_class == 'singlet' & myRCTD@results$results_df$first_type == 'CAF' )]
MAX_EXPR <- 7
my_barc <- rownames(myRCTD@results$results_df)[which((myRCTD@results$results_df$first_type == 'monocyte/DC' & myRCTD@results$results_df$spot_class != 'reject') | (myRCTD@results$results_df$second_type == 'monocyte/DC' & myRCTD@results$results_df$spot_class == 'doublet_certain' ))]
plot_puck_continuous(cropped_puck, cancer_sing, Y, ylimit = c(0,MAX_EXPR),
                     title ='plot of EMT Signature') + geom_point(data = myRCTD@spatialRNA@coords[my_barc,], size = 0.1, color = 'red')+
  scale_color_gradientn(colors = pals::brewer.blues(20)[2:20],name = "EMT Signature Expression", labels = c(0,MAX_EXPR),breaks = c(0,MAX_EXPR), limits = c(0,MAX_EXPR))
plot_puck_continuous(cropped_puck, cancer_sing, Y, ylimit = c(0,MAX_EXPR),
                     title ='plot of EMT Signature') +
  scale_color_gradientn(colors = pals::brewer.blues(20)[2:20],name = "EMT Signature Expression", labels = c(0,MAX_EXPR),breaks = c(0,MAX_EXPR), limits = c(0,MAX_EXPR))



# Get a list of barcodes for cells of target_type
# Filter so we have cells in the cropped puck, that are "singlets" or "certain doublets" with first or second type being the target type
target_df = dplyr::filter(doublet_df, (rownames(doublet_df) %in% barcodes) &
                            ((first_type == target_type  & (spot_class != 'reject')) |
                               ((second_type == target_type) & (spot_class == 'doublet_certain'))
                             )
                          )
target_barcodes = rownames(target_df)

# Names are the barcodes, value is a score computed using euclidean distance from the cells of target_type
all_barcodes = barcodes # The cropped puck subset, use rownames(doublet_df) for all barcodes
explanatory.variable = c(rep(0,length(all_barcodes)))
names(explanatory.variable) = all_barcodes

# Calculate proximity score by summing the scores across all cells of target type for each cell in cropped_puck
# Individual scores between a cell and any target cell is calculated as n_i*exp(-d_i/c)
# n_i is the weighted nUMI of the target cell; weighted by the proportion that the pixel is the target cell type. singlets are weighted as 1.0
# d_i is the distance between the current cell and target cell
# c is a rate constant around 30-50 microns

# Create a distance table between all pairs of cells. rdist is so fast there's no need to save this.
# fields::rdist treats rows as coordinates and computes all distances between placing them in a distance matrix.
dist_matrix = fields::rdist(as.matrix(cropped_puck@coords))
rownames(dist_matrix) = rownames(cropped_puck@coords)
colnames(dist_matrix) = rownames(cropped_puck@coords)
# Precompute the exponent component of the proximity score for all pairs of cells
c = 30/.65
exponent_mat = exp(-dist_matrix/c)

# Precompute the weighted nUMI values for all target cells
weighted_nUMIs = c(rep(0,length(target_barcodes)))
names(weighted_nUMIs) = target_barcodes
for(i in 1:length(weighted_nUMIs)) {
  barcode = target_barcodes[i]
  nUMI = cropped_puck@nUMI[barcode]

  spot_class = doublet_df[barcode,"spot_class"]
  first_type = doublet_df[barcode,"first_type"]
  second_type = doublet_df[barcode,"second_type"]

  weight = 0.0;
  if(spot_class == "singlet") {
    weight = if (first_type == target_type) 1.0 else 0.0;
  } else {
    weight = if (first_type == target_type) weights_doublet[barcode,1] else weights_doublet[barcode,2];
  }
  weighted_nUMI = nUMI * weight
  weighted_nUMIs[i] = weighted_nUMI
}

# Use the precomputed components above to compute explanatory.variable
for(i in 1:length(all_barcodes)) {
  barcode = all_barcodes[i]

  exp_dists = exponent_mat[barcode,target_barcodes]
  proximity_score = weighted_nUMIs %*% exp_dists
  explanatory.variable[i]=proximity_score
}

# Normalize explanatory.variable so the values span from 0 to 1.
normalize_ev = function(ev) {
  # Threshold values over the specific 85% percentile to be 1.
  percentile = quantile(explanatory.variable,.85)
  ev = ev - min(ev)
  ev = ev / percentile
  ev[ev>1] = 1
  return(ev)
}
explanatory.variable = normalize_ev(explanatory.variable)
hist(explanatory.variable)
plot_puck_continuous(cropped_puck, colnames(cropped_puck@counts) , explanatory.variable, # ylimit = c(0,round(quantile(explanatory.variable,0.9))),
                     title ='plot of explanatory.variable cropped_puck')

saveRDS(explanatory.variable, file.path(resultsdir, paste0('exvar',target_type,'.rds')))

# run DE. Won't run on all genes, only the highly expressed ones in both spatialRNA and reference
# CSIDE automatically checks for what cell types to use but since ev is distance from monocytes, it doesn't make sense to measure DE of monocytes so we have a custom list
cell_types = c("CAF","hepatocyte 2","vascular smooth mc")
cell_types_present = c("CAF","hepatocyte 2","vascular smooth mc", 'monocyte/DC')
myRCTD@config$max_cores <- 4
myRCTD <- run.de.single(myRCTD, explanatory.variable, cell_types=cell_types, datadir = resultsdir, cell_type_threshold = 100, cell_types_present = cell_types_present) # outputs plots into pwd() + "/de_plots/"

saveRDS(myRCTD, paste0(resultsdir,'myRCTDde.rds'))

