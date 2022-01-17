library(Matrix)
library(ggplot2)
library(STutility)
library(shiny)
library(imager)
library(raster)
library(rgdal)
library(spacexr)
###CONSTANTS
IM_DIR <- '../RCTD/data/Images/Plaque/210715/Processed/'
ALIGN_DIR <- '../RCTD/data/Images/Plaque/210715/AlignmentResults/'
####START CREATE PUCK IMAGES
datadir <- '../RCTD/data/SpatialRNA/Puck_210605/24'
puck_no <- 'Puck_210605_24'
puck_im_name <- file.path(datadir,paste0(puck_no,"_hippo_princ.tif"))
myRCTD<- readRDS(file.path(datadir,'myRCTD_j20.rds'))
puck <- readRDS(file.path(datadir,'puckCropped.rds'))
create_RCTD_plots(myRCTD,datadir)
use_thresh <- T
source('analysis/helper_functions/alzheimers_helper.R')
# define fns

plot_puck_continuous(puck,rownames(puck@coords), puck@nUMI, ylimit = c(0,1500))
X_SHIFT <- 800; Y_SHIFT <- 0 # 1:(800, 0), 2: (300,1000), 3:(750,750), 4:(0,400)
saveRDS(c(X_SHIFT, Y_SHIFT), file = file.path(datadir, 'xy_shift.rds'))
create_slideseq_im(cell_type_list = c('CA1','CA3','Dentate'), puck_no, datadir, puck, myRCTD,
                   X_SHIFT = X_SHIFT, Y_SHIFT = Y_SHIFT, oligo = F);
create_slideseq_im(cell_type_list = c('Oligodendrocytes'), puck_no, datadir, puck, myRCTD,
                   X_SHIFT = X_SHIFT, Y_SHIFT = Y_SHIFT, oligo = T);
#create dapi (repeat for both)
IMAGE_NUMBER = 5
thresh_im_name <- paste0(IM_DIR,'C1-AVG_slide ',IMAGE_NUMBER,'_thresh.png')
create_thresh_image()
####END CREATE PUCK IMAGES
use_cropped = F
aligned_res <- align_images_now(thresh_im_name, puck_im_name, use_cropped = use_cropped, pre_cropped = T)
aligned <- aligned_res$aligned; se <- aligned_res$se
#se <- MaskImages(object = se)
#image_processing.R > ManualAlignImages.Staffli -> custom.edge.detector > grid.from.staffli
#2D_point_patterns.R > grid.from.staffli > scatter_HE
#verification
#plot(my_im2_orig)
#plot(transform_plaque(aligned, my_im2_orig, extra_transform = F))
verify_alignment(aligned, se, aligned_res$my_im2_orig, use_cropped = use_cropped)
saveRDS(aligned,file.path(ALIGN_DIR,paste0('aligned_',IMAGE_NUMBER)))
saveRDS(se,file.path(ALIGN_DIR,paste0('se_',IMAGE_NUMBER)))
aligned <- readRDS(file.path(ALIGN_DIR,paste0('aligned_',IMAGE_NUMBER)))
se <- readRDS(file.path(ALIGN_DIR,paste0('se_',IMAGE_NUMBER)))
plaque_img <- as.cimg(raster(paste0(IM_DIR,'C2-AVG_slide ',IMAGE_NUMBER,'.tif')))
dapi_img <- as.cimg(raster(paste0(IM_DIR,'C1-AVG_slide ',IMAGE_NUMBER,'.tif')))
plot(pmin(dapi_img,1000))
#hist(log(plaque_img,2), breaks = 100, xlim=c(6,10))
#plaque_img <- pmax(plaque_img - 100,0)
plot(pmin(plaque_img,1000))

#plot(pmax(pmin(wplaque_img,200),1000))
#plot(pmax(pmin(wplaque_img - 100,5),0))
#plot(pmin(pmax(plaque_img-100,0),500))
