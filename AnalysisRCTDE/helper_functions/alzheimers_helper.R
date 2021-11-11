align_images_now <- function(thresh_im_name, puck_im_name, use_cropped = F, pre_cropped = F, thresh = 150, thresh2 = 150) {
  samples <- list.files(pattern = "^Hippo[1-6].tsv.gz", path = system.file("extdata/counts", package = "STutility"), full.names = TRUE)
  spotfiles <- list.files(pattern = "^alignment_table_Hippo[1-6].tsv.gz", path = system.file("extdata/counts", package = "STutility"), full.names = TRUE)
  imgs <- list.files(pattern = "^Hippo[1-6].jpg", path = system.file("extdata/counts", package = "STutility"), full.names = TRUE)
  imgs[6] <- thresh_im_name
  imgs[1] <- puck_im_name
  imgs[2] <- file.path(datadir,paste0(puck_no,"_hippo_oligo.tif"))
  infoTable <- data.frame(samples, imgs, spotfiles, stringsAsFactors = F)
  
  se <- InputFromTable(infotable = infoTable, 
                       min.gene.count = 100, 
                       min.gene.spots = 5,
                       min.spot.count = 500,
                       platform =  "Visium")
  se <- LoadImages(se, time.resolve = FALSE, verbose = TRUE)
  ImagePlot(se, method = "raster", type = "raw")
  
  my_im <- se@tools$Staffli['raw'][[1]] %>% as.cimg # this is type = 'masked' or type = 'raw'
  sum_special <- sum(my_im)
  my_im2 <- se@tools$Staffli['raw'][[6]] %>% as.cimg # this is type = 'masked' or type = 'raw'
  my_im2_orig <- my_im2
  if(use_cropped) {
    if(pre_cropped)
      my_im2 <- readRDS(paste0(IM_DIR,'cropped/DAPI_',IMAGE_NUMBER,'.rds'))
    else
      my_im2 <- crop_dapi_im(my_im2)
  }
  se@tools$Staffli['raw'][[6]] <- my_im2
  my_im[my_im < thresh] <- 0
  plot(my_im)
  my_im2[my_im2 < thresh2] <- 0
  plot(my_im2)
  custom.edge.detector <- function(my_im) {
    if(sum(my_im) == sum_special)
      return(my_im > thresh)
    return(my_im > thresh2)
  }
  aligned <- ManualAlignImages(se, type = 'raw', custom.edge.detector = custom.edge.detector,edges = F)
  return(list(aligned = aligned, se = se, my_im2_orig = my_im2_orig))
}

crop_dapi_im <- function(my_im2) {
  library(gatepoints)
  my_df <- as.data.frame(grayscale(my_im2))
  X <- my_df[,c(1,2)]
  pal = colorRampPalette(c("blue", "red"))
  dforder = findInterval(my_df$value, sort(my_df$value))
  plot(X, pch=20, col=pal(nrow(X))[dforder], cex = 0.5)
  selectedPoints <- fhs(X, mark = TRUE,  pch=20, col=pal(nrow(X))[dforder], cex = 0.5)
  arr <- as.integer(rownames(X) %in% selectedPoints); names(arr) <- rownames(X)
  my_im2[!arr] <- 0
  plot(my_im2)
  saveRDS(my_im2, paste0(IM_DIR,'cropped/DAPI_',IMAGE_NUMBER,'.rds'))
  return(my_im2)
}

verify_alignment <- function(algined, se, my_im2_orig, use_cropped = F) {
  m1 <- aligned@tools$Staffli['processed']$`1`
  #m2 <- aligned@tools$Staffli['processed']$`6`
  if(use_cropped)
    m2 <- transform_plaque(aligned, my_im2_orig, extra_transform = F) #
  else
    m2 <- aligned@tools$Staffli['processed']$`6`
  m3 <- se@tools$Staffli['raw']$`2`
  s2 <- (as.cimg(m2))
  s1 <- (as.cimg(m1))
  s3 <- as.cimg(m3)
  B(s1) <- B(s2)
  G(s1) <- G(s3)
  R(s1) <- R(s1)
  plot(s1)
}

create_thresh_image <- function() {
  str_name<-paste0(IM_DIR,'C1-AVG_slide ',IMAGE_NUMBER,'.tif')
  imported_raster=raster(str_name)
  my_mat <- as.matrix(imported_raster)
  blur_im <- pmin(my_mat,1000) %>% as.cimg #7
  plot(blur_im)
  save.image(blur_im, thresh_im_name)
}

create_slideseq_im <- function(cell_type_list = c('CA1','CA3','Dentate'), puck_no, datadir, puck, myRCTD,  X_SHIFT = 0, Y_SHIFT = 0, oligo = F) {
  decay_scale = 3 #30
  #pic <- matrix(0,4401,4401)
  pic <- matrix(0,4401,4401)
  MAX_UMI = 3000;
  ind_list <- (1:length(puck@nUMI))[names(puck@nUMI) %in% rownames(myRCTD@results$results_df)[which(myRCTD@results$results_df$first_type %in% c('Oligodendrocyte') & myRCTD@results$results_df$spot_class != 'reject')]]
  for(ind in ind_list) {#1:length(puck@nUMI)) {
    print(ind)
    cur_UMI = min(puck@nUMI[ind],MAX_UMI)
    x = puck@coords[ind,'x'] - X_SHIFT
    y = puck@coords[ind,'y'] - Y_SHIFT
    for(dx in ((-4*decay_scale) : (4*decay_scale)))
      for(dy in ((-4*decay_scale) : (4*decay_scale)))
        pic[x + dx, y + dy] = pic[x + dx, y + dy] + cur_UMI*exp(-(dx^2 + dy^2)/(2*decay_scale^2))
  }
  pic_rev <- pmin(pic, 1)
  image(pic_rev, useRaster=TRUE, axes=FALSE)
  
  l <- list(x = 1:nrow(pic_rev), y = 1:ncol(pic_rev), z = pic_rev)
  x <- image2Grid(l)
  if(!oligo) {
    writeGDAL(x, file.path(datadir,paste0(puck_no,"_hippo_princ.tif")))
  } else
    writeGDAL(x, file.path(datadir,paste0(puck_no,"_hippo_oligo.tif")))
}

transform_plaque <- function(aligned, plaque_img, extra_transform = T) {
  #micro_img <- as.cimg(raster('../../RCTD/data/Images/Plaque/210715/Processed/C1-AVG_slide 2_thresh.png'))
  #micro_img <- as.cimg(raster('../../RCTD/data/Images/Plaque/210715/Processed/C1-AVG_slide 2.tif'))
  #plot(micro_img)
  #apply transformation to plaqeu
  #plot(micro_img)
  #wimg <- imwarp(micro_img, map.shift)
  #wimg[wimg > 2] <- 2
  #plot(wimg)
  my_t_form <- aligned@tools$Staffli@transformations[[6]]
  my_t_form[c(1,2),3] <- my_t_form[c(1,2),3]*dim(plaque_img)[1]/400
  map.shift <- function(x,y) {
    prod <- (solve(my_t_form) %*% rbind(x,y,1))
    return(list(x = prod[1,], y = prod[2,]))
  }
  #t_img <- plaque_img
  #t_img[t_img < 128] <- 0
  #plot(pmin(t_img,1000))
  #plaque_img <- pmin(plaque_img,500)
  #plot(plaque_img)
  #plot(imager::mirror(plaque_img, 'x'))
  if(extra_transform)
    plaque_img <- imrotate(imager::mirror(plaque_img, 'x'),-90)
  plot(pmin(plaque_img,1000))
  wplaque_img2 <- imwarp(plaque_img, map.shift)
}

get_new_coords <- function(aligned, puck, plaque_img, X_SHIFT = 0, Y_SHIFT = 0, only_neurons = T) {
  my_t_form <- aligned@tools$Staffli@transformations[[6]]
  my_t_form[c(1,2),3] <- my_t_form[c(1,2),3]*dim(plaque_img)[1]/400
  if(only_neurons)
    ind_list <- (1:length(puck@nUMI))[names(puck@nUMI) %in% rownames(myRCTD@results$results_df)[which(myRCTD@results$results_df$first_type %in% c('CA1','CA3','Denate') & myRCTD@results$results_df$spot_class != 'reject')]]
  else
    ind_list <- rownames(puck@coords)
  new_coords <- puck@coords[ind_list,]
  new_coords$x <- (new_coords$x - X_SHIFT)
  new_coords$y <- (new_coords$y - Y_SHIFT)
  new_coords <- new_coords/4401*dim(plaque_img)[1]
  new_coords2 <- new_coords
  new_coords2$y <- new_coords$x
  new_coords2$x <- dim(plaque_img)[1] - new_coords$y
  new_coords <- new_coords2
  
  map.shift.inv_ss <- function(x) {
    prod <- (my_t_form %*% c(x[2],x[1],1))
    return(c(prod[1,],prod[2,]))
  }
  
  rev_2 <- apply(new_coords,1, map.shift.inv_ss)
  rev_2 <- t(matrix(rev_2,nrow = 2))
  colnames(rev_2) <- c('x','y')
  rownames(rev_2) <- rownames(new_coords)
  rev_2 <- as.data.frame(rev_2)
  return(rev_2)
}

verify_puck_transform <- function(new_puck, dapi_img) {
  my_m <- matrix(0, nrow = dim(dapi_img)[1], ncol = dim(dapi_img)[2])
  for(ind in rownames(new_puck@coords))
    for(dx in -5:5)
      for(dy in -5:5)
        my_m[new_puck@coords[ind,'x'] + dx, new_puck@coords[ind, 'y'] + dx] <- 1
      
  s1 <- raster(my_m*1000) %>% as.cimg
  s2 <- pmin(dapi_img,1000)
  s1 <- s1/max(s1)
  s2 <- s2 / max(s2)
  s3 <- add.color(s1)
  B(s3) <- s2
  G(s3) <- 0
  R(s3) <- s1
  plot(s3)
}


get_plaque_density <- function(plaque_img, a_len = 30) {
  plaque_img <- pmax(plaque_img - 100,0)
  my_im_t <- plaque_img
  my_im_t[my_im_t < 300] <- 0
  plot(pmin(my_im_t,1000))
  #plot(my_im_t > 300)
  D <- round(a_len*6.65)*2+1
  unpadded <- as.matrix(my_im_t)
  padded <- matrix(0, nrow = dim(my_im_t)[1] + 2*(D-1)/2, ncol = dim(my_im_t)[2] + 2*(D-1)/2)
  padded[((D+1)/2):((D-1)/2+dim(my_im_t)[1]), ((D+1)/2): (dim(my_im_t)[2] + (D-1)/2)] <- unpadded
  rast <- raster(padded)
  
  cent <- (D+1) / 2
  f = matrix(0, nrow=D, ncol = D)
  for(i in 1:D)
    for(j in 1:D)
      f[i,j] <- exp(-max((sqrt((i - cent)^2 + (j - cent)^2) / a_len) - 1,0))
  print(f[cent,cent] / f[cent,D])
  f <- f / sum(f)
  plot(f %>% as.cimg)
  system.time({rm <- raster::focal(rast, w=f, pad = F)})
  padded_rm <- raster(as.matrix(rm)[((D+1)/2):((D-1)/2+dim(my_im_t)[1]), ((D+1)/2): (dim(my_im_t)[2] + (D-1)/2)])
  plot(pmin(plaque_img,1000))
  plot(pmin(plaque_img %>% as.matrix %>% raster %>% as.cimg, 1000))
  quant_thresh <- .950
  quant_cutoff <- quantile(padded_rm, quant_thresh)
  hist(pmin(padded_rm %>% as.cimg, quant_cutoff))
  plot(pmin(t(padded_rm) %>% as.cimg, quant_cutoff))
  return(list(p_img = t(padded_rm), quant_cutoff = quant_cutoff))
}



