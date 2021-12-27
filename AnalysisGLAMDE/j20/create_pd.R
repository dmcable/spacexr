library(Matrix)
library(imager)
library(raster)
library(rgdal)
source('~/Documents/MIT/Research/Rafalab/Projects/slideseq/Cell Demixing/ContentStructure/DEGLAM/analysis/helper_functions/alzheimers_helper.R')

IM_DIR <- '~/Documents/MIT/Research/Rafalab/Projects/slideseq/Cell Demixing/ContentStructure/RCTD/data/Images/Plaque/210715/Processed/'
for(IMAGE_NUMBER in 8:1) {
  print(IMAGE_NUMBER)
  plaque_img <- as.cimg(raster(paste0(IM_DIR,'C2-AVG_slide ',IMAGE_NUMBER,'.tif')))
  pd_res <- get_plaque_density(plaque_img, a_len = 25)
  saveRDS(pd_res, file = file.path(IM_DIR, paste0('Plaque_Density/',IMAGE_NUMBER,'.rds')))
}

