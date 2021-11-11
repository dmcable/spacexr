library(ijtiff)
my_dir <- '/Users/dcable/Documents/MIT/Research/Rafalab/Projects/slideseq/Cell Demixing/ContentStructure/RCTD/data/Images'
im_ind <- 5
mtif <- read_tif(file.path(my_dir,paste0('HCR/Composite',im_ind,'.tif')))
atif <- read_tif(file.path(my_dir,paste0('HCR/c2_i',im_ind,'.tif')))
if(im_ind == 1) {
  purk_thresh <- 200
  berg_thresh <- 40000
  min_x_n <- 5300
  min_y_n <- 2875
  max_y_n <- 3650
  max_x_n <- 6400
  min_y_l <- 3800
  min_x_l <- 1600
  max_x_l <- 4550
  max_y_l <- 5100
  pindex <- 2
  bindex <- 4
}
if(im_ind == 2) {
  purk_thresh <- 25000
  berg_thresh <- 36000
  max_x_n <- 2500
  min_y_n <- 5376
  min_x_n <- 1700
  max_y_n <- 6700
  max_x_l <- 1320
  min_y_l <- 3700
  min_x_l <- 350
  max_y_l <- 4800
  pindex <- 2
  bindex <- 4
}
if(im_ind == 3) {
  purk_thresh <- 25000
  berg_thresh <- 41000
  min_x_n <- 5400
  min_y_n <- 3468
  max_x_n <- 6500
  max_y_n <- 4400
  max_x_l <- 4362
  min_y_l <- 3990 
  min_x_l <- 2850
  max_y_l <- 5200
  pindex <- 2
  bindex <- 4
}
if(im_ind == 4) {
  purk_thresh <- 25000
  berg_thresh <- 14500
  max_x_n <- 3282
  max_y_n <- 1494
  min_x_l <- 4356
  max_y_l <- 1974
  pindex <- 2
  bindex <- 3
  min_x_n <- 2150
  min_y_n <- 800
  max_x_l <- 5600
  min_y_l <- 950
} 
#if(im_ind == 5) {
#  purk_thresh <- 20000
#  berg_thresh <- 47000
#  min_x_n <- 3570
#  min_y_n <- 4650
#  max_x_l <- 2574
#  min_y_l <- 4038
#  pindex <- 3
#  bindex <- 4
#}
purk_ind <- which(mtif[,,pindex,1] >= purk_thresh)
berg_ind <- which(mtif[,,bindex,1] >= berg_thresh)
both_ind <- intersect(purk_ind, berg_ind)
berg_ind <- setdiff(berg_ind, both_ind)
purk_ind <- setdiff(purk_ind, both_ind)
y_im <- mtif[,,3,1]
y_im[] <- rep(1:dim(mtif)[1],dim(mtif)[2])
x_im <- t(y_im)
x_im[] <- rep(1:dim(mtif)[2],dim(mtif)[1])
x_im <- t(x_im)

y_index_l <- which(y_im < max_y_l & y_im > min_y_l)
y_index_n <- which(y_im < max_y_n & y_im > min_y_n)
l_index <- intersect(y_index_l, which(x_im < max_x_l & x_im > min_x_l))
n_index <- intersect(y_index_n, which(x_im < max_x_n & x_im > min_x_n))

l_purk_ind <- intersect(l_index,purk_ind)
n_purk_ind <- intersect(n_index,purk_ind)
length(n_purk_ind) / length(purk_ind)
length(l_purk_ind) / length(purk_ind)
plot(x_im[l_purk_ind],-y_im[l_purk_ind])
plot(x_im[n_purk_ind],-y_im[n_purk_ind])
bg_aldoc_im <- quantile(atif, 0.5)
m_pn <- mean(atif[n_purk_ind]) - bg_aldoc_im
m_pl <- mean(atif[l_purk_ind]) - bg_aldoc_im
log_fc_p <- log(m_pn) - log(m_pl)

l_berg_ind <- intersect(l_index,berg_ind)
n_berg_ind <- intersect(n_index,berg_ind)
length(n_berg_ind) / length(berg_ind)
length(l_berg_ind) / length(berg_ind)
m_bn <- mean(atif[n_berg_ind]) - bg_aldoc_im
m_bl <- mean(atif[l_berg_ind]) - bg_aldoc_im
log_fc_b <- log(m_bn) - log(m_bl)

print(log_fc_p)
print(log_fc_b)


