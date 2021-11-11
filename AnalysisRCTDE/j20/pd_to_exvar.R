library(Matrix)
library(ggplot2)
library(STutility)
library(shiny)
library(imager)
library(raster)
library(rgdal)
library(DEGLAM)
###CONSTANTS
IM_DIR <- '../RCTD/data/Images/Plaque/210715/Processed/'
ALIGN_DIR <- '../RCTD/data/Images/Plaque/210715/AlignmentResults/'
####START CREATE PUCK IMAGES
datadir <- '../RCTD/data/SpatialRNA/Puck_210605/24'
puck_no <- 'Puck_210605_24'
puck_im_name <- file.path(datadir,paste0(puck_no,"_hippo_princ.tif"))
myRCTD<- readRDS(file.path(datadir,'myRCTD_j20.rds'))
puck <- readRDS(file.path(datadir,'puckCropped.rds'))
use_thresh <- T
source('analysis/helper_functions/alzheimers_helper.R')
# define fns
###Only for the first
explanatory_variable <- as.data.frame(cbind(puck@nUMI[colnames(puck@counts)]*0,puck@nUMI[colnames(puck@counts)]*0))
rownames(explanatory_variable) <- colnames(puck@counts)
pd_res_old <- list()
plot_puck_continuous(puck,rownames(puck@coords), puck@nUMI, ylimit = c(0,1500))

xy_shift <- readRDS(file = file.path(datadir, 'xy_shift.rds'))
X_SHIFT <- xy_shift[1]; Y_SHIFT <- xy_shift[2]
### transform slideseq to plaque (repeat for both IMAGE_NUMBERS)
IMAGE_NUMBER = 8
plaque_img <- as.cimg(raster(paste0(IM_DIR,'C2-AVG_slide ',IMAGE_NUMBER,'.tif')))
aligned <- readRDS(file.path(ALIGN_DIR,paste0('aligned_',IMAGE_NUMBER)))
dapi_img <- as.cimg(raster(paste0(IM_DIR,'C1-AVG_slide ',IMAGE_NUMBER,'.tif')))
new_puck <- puck
new_puck@coords <- get_new_coords(aligned, puck, plaque_img, X_SHIFT = X_SHIFT, Y_SHIFT = Y_SHIFT, only_neurons = T)
plot_puck_continuous(new_puck, rownames(new_puck@coords), new_puck@nUMI, ylimit = c(0,1500))
verify_puck_transform(new_puck, dapi_img)
#pd_res <- get_plaque_density(plaque_img, a_len = 25)
if(use_thresh) {
  pd_res <- readRDS(file.path(IM_DIR, paste0('Plaque_Density/',IMAGE_NUMBER,'_thresh.rds')))
} else {
  pd_res <- readRDS(file.path(IM_DIR, paste0('Plaque_Density/',IMAGE_NUMBER,'.rds')))
}
plot(pmin(pd_res$p_img %>% as.cimg, pd_res$quant_cutoff))

new_puck@coords <- get_new_coords(aligned, puck, plaque_img, X_SHIFT = X_SHIFT, Y_SHIFT = Y_SHIFT, only_neurons = F)
modulus <- IMAGE_NUMBER %% 2
if(modulus == 0)
  modulus = 2
for(name_ind in rownames(new_puck@coords))
  explanatory_variable[name_ind, modulus] <- pd_res$p_img[round(new_puck@coords[name_ind,'x']),round(new_puck@coords[name_ind,'y'])]

explanatory.variable <- explanatory_variable[,modulus]; names(explanatory.variable) <- rownames(explanatory_variable)
hist(pmin(explanatory.variable,pd_res$quant_cutoff))
mean(pmin(explanatory.variable,pd_res$quant_cutoff) > pd_res$quant_cutoff/2)
plot_puck_continuous(puck,rownames(puck@coords), explanatory.variable, ylimit = c(0,pd_res$quant_cutoff))
#plot(pmax(pmin(wplaque_img - 100,100),0)) #compare

pd_res_old[[modulus]] <- pd_res$p_img

### GO BACK UP
if(use_thresh) {
  saveRDS(pd_res_old, file = file.path(datadir,'pd_res_old_thresh.rds'))
} else {
  saveRDS(pd_res_old, file = file.path(datadir,'pd_res_old.rds'))
}
explanatory_variable_norm <- explanatory_variable
for(ind in 1:2) {
  q_level <- quantile(explanatory_variable[,ind],.90)
  explanatory_variable_norm[,ind] <- pmin(explanatory_variable_norm[,ind], q_level) / q_level
}
exvar <- rowMeans(explanatory_variable_norm)
hist(exvar)
mean(exvar > 0.5)
plot_puck_continuous(puck,rownames(puck@coords), exvar, ylimit = c(0,1))
if(use_thresh) {
  save(exvar, explanatory_variable, file = file.path(datadir,'exvar_thresh.RData'))
} else {
  save(exvar, explanatory_variable, file = file.path(datadir,'exvar.RData'))
}

