puck <- readRDS(paste0(iv$slideseqdir,"/SplitPuck/puck",fold_index,".RDS"))
puck <- iv$puck
cell_type_info_renorm <- iv$cell_type_info
ind = 5
bead = puck@counts[,ind]
ind2 = 73
bead2 = puck@counts[,ind2]
results = process_bead(cell_type_info_renorm, iv$gene_list, 1000, bead, constrain = T)
decompose_full(cell_type_info_renorm[[1]], iv$gene_list, 1000, bead, constrain = F, verbose = F)$weights

type1 = "Endothelial"
type2 = "Astrocytes"
Q_mat <- get_Q_mat()
bead = bead_mix(test_reference, iv$gene_list, 500,500, type1, type2)
bead = bead_mix(restrict_test_ref, iv$gene_list, UMI1, UMI2, type1, type2)
#bead2 = bead_mix(test_reference, iv$gene_list, 300,700, type1, type2)
results = process_bead(cell_type_info_renorm, iv$gene_list, 1000, bead[,1], constrain = F)
print(results)
UMI_tot = 1000
get_prediction_sparse(cell_type_info_renorm[[1]], iv$gene_list, UMI_tot, 1, type1, type2)

#cell_types = c("Purkinje","MLI2")
reg_data = cell_type_info_renorm[[1]][iv$gene_list,] * 1000
reg_data = data.matrix(reg_data)
reg_data = reg_data[,cell_types]
prediction_bad = reg_data %*% c(0.73,0.27)
get_log_l <- function(gene_list, prediction, bead) {
  sigma = 0.95
  N = 50
  delta = 1e-1*sigma
  log_l = 0
  index = 0
  for(gene in gene_list) {
    S = 0
    for(n in (-N:N)) {
      z = n*delta
      p = 1/sqrt(2*pi*(sigma^2))*exp(-z^2/(2*sigma^2))
      u = exp(z) * max(prediction[gene,],1e-4)
      Y = bead[gene]
      log_x = -u + Y * log(u) - lgamma(Y+1)
      #x = exp(-u)*(u)^Y/factorial(Y)
      x = exp(log_x)
      #x = 1
      S = S + x*delta*p
    }
    log_l = log_l + log(S)
    if(is.nan(log(S)))
      print(log(S))
    index = index+1
  }
  return(log_l)
}


get_likelihood(iv$gene_list, prediction, bead)
get_likelihood(iv$gene_list, prediction_bad, bead)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   

barcodes = (1:75)*10

change_likelihood = numeric(length(barcodes)); names(change_likelihood) = barcodes
for (barcode in barcodes) {
  print(barcode)
  type1 = results_df[barcode, "first_type"]
  type2 = results_df[barcode, "second_type"]
  bead <- puck@counts[,barcode]
  prediction <- get_prediction_sparse(cell_type_info_renorm[[1]], iv$gene_list, UMI_tot, 1, type1, type2)
  log_l <- get_log_l(iv$gene_list, prediction, bead)
  p <- weights_doublet[barcode,"first_type"]
  prediction_doub <- get_prediction_sparse(cell_type_info_renorm[[1]], iv$gene_list, UMI_tot, p, type1, type2)
  log_l_doub <- get_log_l(iv$gene_list, prediction_doub, bead)
  change_likelihood[barcode] = log_l_doub - log_l
}

singlets = meta_df[barcodes, ]$first_UMI == 0 | meta_df[barcodes, ]$first_UMI == 1000
sing_barcode = barcodes[which(singlets)]; doub_barcode = barcodes[which(!singlets)]
hist(change_likelihood[doub_barcode], breaks = 20, xlim = c(-30,200))
hist(change_likelihood[sing_barcode], breaks = 20, xlim = c(-30,200))
p1 <- hist(change_likelihood[sing_barcode])
p2 <- hist(change_likelihood[doub_barcode])
plot( p1, col=rgb(0,0,1,1/4))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), add=T)
res =  change_likelihood[doub_barcode]
names(res) = doub_barcode
sum(change_likelihood[sing_barcode] < 25) / length(sing_barcode)
sum(change_likelihood[doub_barcode] < 25) / length(doub_barcode)
singlets = meta_df[barcodes, ]$first_UMI == 0 | meta_df[barcodes, ]$first_UMI == 1000
skew_doublets = meta_df[barcodes, ]$first_UMI == 250 | meta_df[barcodes, ]$first_UMI == 750
full_doublets = meta_df[barcodes, ]$first_UMI == 500
skew_barcode = barcodes[which(skew_doublets)]; full_barcode = barcodes[which(full_doublets)]
prop_doublets = numeric(3); names(prop_doublets) = c(0,250,500)
prop_doublets["0"] = sum(change_likelihood[sing_barcode] < 25) / length(sing_barcode)
prop_doublets["250"] = sum(change_likelihood[skew_barcode] < 25) / length(skew_barcode)
prop_doublets["500"] = sum(change_likelihood[full_barcode] < 25) / length(full_barcode)
plot(names(prop_doublets),prop_doublets,ylim=c(0,1),xlab="nUMI",ylab="proportion singlets")
barcodes = (1:750)
#calculates the change in likelihood between doublet and singlet


change_likelihood_q = numeric(length(barcodes)); names(change_likelihood) = barcodes
for (barcode in barcodes) {
  print(barcode)
  type1 = results_df[barcode, "first_type"]
  type2 = results_df[barcode, "second_type"]
  bead <- puck@counts[,barcode]
  prediction <- get_prediction_sparse(cell_type_info_renorm[[1]], iv$gene_list, UMI_tot, 1, type1, type2)
  log_l <- calc_log_l_par(iv$gene_list, prediction, bead)
  p <- weights_doublet[barcode,"first_type"]
  prediction_doub <- get_prediction_sparse(cell_type_info_renorm[[1]], iv$gene_list, UMI_tot, p, type1, type2)
  log_l_doub <- calc_log_l_par(iv$gene_list, prediction_doub, bead)
  change_likelihood_q[barcode] = log_l_doub - log_l
}
hist(-change_likelihood_q[doub_barcode], breaks = 20, xlim = c(-30,200))
hist(-change_likelihood_q[sing_barcode], breaks = 20, xlim = c(-30,200))
sum(-change_likelihood_q[sing_barcode] < 25) / length(sing_barcode)
sum(-change_likelihood_q[doub_barcode] < 25) / length(doub_barcode)

#truth
prediction <- get_prediction_sparse(cell_type_info_renorm[[1]], iv$gene_list, UMI_tot, 0.5, "Astrocytes","Granule")
log_l <- get_log_l(iv$gene_list, prediction, bead)

#make the residual plots

barcodes = rownames(meta_df); M = length(barcodes)
predictions = matrix(0, nrow = length(iv$gene_list), ncol = M)
colnames(predictions) = barcodes; rownames(predictions) = iv$gene_list
for (barcode in barcodes) {
  type1 = meta_df[barcode, "first_type"]
  type2 = meta_df[barcode, "second_type"]
  p = meta_df[barcode, "first_UMI"] / UMI_tot
  predictions[, barcode] <- get_prediction_sparse(cell_type_info_renorm[[1]], iv$gene_list, UMI_tot, p, type1, type2)
}
X = as.vector(predictions)
Y = as.vector(puck@counts)
my_ind = X > 5 & X < 6
mean(X[my_ind])
mean(Y[my_ind])
plot(density(Y[my_ind]))
#look at the dropviz residuals
cell_type_means_renorm = cell_type_info_renorm[[1]]
myref <- create_downsampled_data(test_reference, NULL, n_samples = 250, save.file=F)
counts <- myref@assays$RNA@counts[iv$gene_list,] #counts
res <- apply(myref@meta.data[,c('liger_ident_coarse','nUMI')], 1, function(x) cell_type_means_renorm[,x[1]]*as.numeric(x[2]))
umi = t(replicate(length(iv$gene_list),myref@meta.data$nUMI))
N_UMI = as.vector(umi)
X = as.vector(res)
Y = as.vector(counts)
data_df <- data.frame(X,Y,N_UMI)
index = 2048680
x = floor(index / 3130) + 1
y = index %% 3130
res[y, x]
counts[y,x]
iv$gene_list[y]
myref@meta.data$liger_ident_coarse[x]
cell_type_means_renorm[y,]


bad_list = as.numeric(rownames(data_df[(log(Y + 1e-4) - log(X + 1e-4) > 7) & Y > 2,]))
bad_list_gene = iv$gene_list[bad_list %% length(iv$gene_list)]
bad_list_type = myref@meta.data$liger_ident_coarse[floor(bad_list / length(iv$gene_list)) + 1]
my_ind = X > 20 & X < 22
mean(X[my_ind])
mean(Y[my_ind])
plot(density(log(Y[my_ind])))

my_ind = (X > 2) & (Y > 1)
Y_my = Y[my_ind]; X_my = X[my_ind]
Z = log(Y_my + 1e-4) - log(X_my + 1e-4)
# new epsilon = max(1e-4, 1e-7*nUMI)
plot(density(Z))
hist(Z,breaks=20)
sigma = 0.8 #sqrt(var(Z))
trials = 60
delta = 0.1
v = (1:trials) * delta
obs_cdf = numeric(trials)
norm_cdf = numeric(trials)
new_cdf = numeric(trials)
for (i in 1:trials) {
  x = v[i]
  obs_cdf[i] = sum(Z > sigma*x)/length(Z)
  norm_cdf[i] = pnorm(-x)
  if(x >= 3)
    new_cdf[i] = a/(x-c)
  else
    new_cdf[i] = a/(3-c) + (pnorm(-x) - pnorm(-3))
}

plot( v, log(obs_cdf), type="l", col="red",xlab = 'x', ylab = 'P(Z > x)' )
lines( v, log(norm_cdf), col="green" )
lines( v, log(new_cdf), col="blue" )

a = 4/9*exp(-3^2/2)/sqrt(2*pi)
c = 7/3
a/(3-7/3)^2

x = -1000:1000/10
y = numeric(length(x))
for(i in 1:length(x)) {
  y[i] = ht_pdf(x[i],1)
}
plot(x,log(y),type="l",xlim=c(-10,10))




N_X = 10000; delta = 1e-5
X = (1:N_X)^1.5*delta
Q_mat <- calc_Q_mat(sigma, X, K = 4)
saveRDS(Q_mat, file.path(resultsdir,"Q_mat.RDS"))
index = 200001; S[index,1]
x = X[1]
Y[index]
exp(-x)*(x^k)/1*p[index]
log_R = -x + log(x)


N_Y = 200000; gamma = 1e-4
Y = (-N_Y:N_Y) * gamma
index = 100
x = X[index]
k = 2
S <- 0
for (y in Y) {
  p = ht_pdf(y,sigma)
  pred = exp(y) * x
  log_S = log(p) - pred + k * log(pred) - log(factorial(k))
  S = S + exp(log_S)
}
S = S*gamma
k = 2
#plot(X[(index-20):(index+20)],Q_mat[(k+1),(index-20):(index+20)])
num_der <- (Q_mat[(k+1), index + 1] - Q_mat[(k+1), index])/(X[index + 1] - X[index])
the_der <- 1/x * (-(k+1)*Q_mat[k+2,index] + k*Q_mat[k+1,index])
the_der_next <- 1/X[index+1] * (-(k+1)*Q_mat[k+2,index+1] + k*Q_mat[k+1,index+1])
emp_sec_der <- (the_der_next - the_der)/(X[index + 1] - X[index])
the_sec_der <- 1/(x^2)*((k+1)*(k+2)*Q_mat[k+3,index] - k*(2*(k+1)*Q_mat[k+2,index] - (k-1)*Q_mat[k+1,index]))
print(emp_sec_der)
print(the_sec_der)



trials = 10000; k = 2
X_check = 1:trials*(1e-5)
y = numeric(trials)
for(i in 1:trials) {
  y[i] = calc_Q(X_check[i], k)
}
plot(X_check,y, type="l",col="red", xlim= c(0,0.3), ylim = c(0,1))
lines(X,Q_mat[k+1,], col="green")

gene_list = iv$gene_list
cell_types = c(type1,type2)
reg_data = cell_type_info_renorm[[1]][iv$gene_list,] * 1000
reg_data = data.matrix(reg_data)
reg_data = reg_data[,cell_types]
S = reg_data; B = bead
solution = as.vector(c(0.4,0.6))
prediction = S %*% solution
calc_log_l_par(gene_list, prediction, B)

b = c(0,0)
H = Matrix(0, nrow = 2, ncol = 2)
for (gene in gene_list) {
  b = b - calc_log_p_d1(prediction[gene,], B[gene,]) * S[gene,]
  H = H - calc_log_p_d2(prediction[gene,], B[gene,]) * S[gene,] %*% t(S[gene,])
}
d = c(-0.000001,0.000001)
prediction_d = S %*% (solution+d)
change_est = b %*% d + 0.5*t(d) %*% H %*% d
change_mes = calc_log_l_par(gene_list, prediction_d, B) - calc_log_l_par(gene_list, prediction, B)
print(change_est); print(change_mes)
print(change_est/change_mes)

b1 = c(0,0)
H1 = Matrix(0, nrow = 2, ncol = 2)
d1_vec = calc_log_d1_par(gene_list, prediction, bead)
d2_vec = calc_log_d2_par(gene_list, prediction, bead)
for (gene in gene_list) {
  b1 = b1 + d1_vec[gene] * S[gene,]
  H1 = H1 + d2_vec[gene] * S[gene,] %*% t(S[gene,])
}
b2 = d1_vec %*% S
H2 = t(S) %*% (diag(d2_vec) %*% S)


b <- get_gradient(S, B, gene_list, prediction)
H <- get_hessian(S, B, gene_list, prediction)A<-cbind(diag(dim(S)[2]))
bzero<- (-solution)
A_const = t(rbind(1,A))
b_const <-c(1 - sum(solution),bzero)
quadprog::solve.QP(H,-b,A_const,b_const,meq=1)$solution

#1D scratch
gene_list = iv$gene_list[1:100]
S = reg_data[gene_list,1]; B = bead[gene_list,]
solution = as.vector(c(1))
prediction = S * solution
calc_log_l_par(gene_list, prediction, B)

b = c(0)
H = Matrix(0, nrow = 1, ncol = 1)
for (gene in gene_list) {
  b = b - calc_log_p_d1(prediction[gene], B[gene]) * S[gene]
  H = H - calc_log_p_d2(prediction[gene], B[gene]) * S[gene] %*% t(S[gene])
}
d = c(-0.000001)
prediction_d = S * (solution+d)
change_est = b * d + 0.5* H * (d^2)
change_mes = calc_log_l_par(gene_list, prediction_d, B) - calc_log_l_par(gene_list, prediction, B)
print(change_est); print(change_mes)
print(change_est/change_mes)

x = 1
k = 1
d = 0.001
obs_d = calc_log_p(x + d,k) - calc_log_p(x,k)
calc_d = calc_log_p_d1(x,k) * d + 0.5*(d^2)*calc_log_p_d2(x,k)
obs_dd = calc_log_p_d1(x+d,k) - calc_log_p_d1(x,k)
calc_dd = d*calc_log_p_d2(x,k)
#end 1D scratch

type1 = "Endothelial"
type2 = "Astrocytes"
bead = bead_mix(test_reference, iv$gene_list, 0,1000, type1, type2)
decompose_sparse(cell_type_means_renorm, gene_list, 1000, bead, type1, type2, verbose = T, constrain = F)

