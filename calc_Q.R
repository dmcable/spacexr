library(RCTD)
library(Matrix)
#calculates Q_mat for one value of i = k + 1, taken from args
iv <- init_RCTD(load_info_renorm = T) #initial variables

N_X = 50000; delta = 1e-5;
X_vals = (1:N_X)^1.5*delta
K = 100
resultsdir <- paste0(iv$slideseqdir,"/results")
sigresdir = file.path(resultsdir, "sigma")
sigma <- readRDS(file.path(sigresdir,"sigma.RDS"))
results <- calc_Q_par(K, X_vals, sigma)
Q_mat <-t(as.data.frame(matrix(unlist(results), nrow=length(unlist(results[1])))))
saveRDS(Q_mat, file.path(resultsdir,'Q_mat.RDS'))


#args = commandArgs(trailingOnly=TRUE)
#if (length(args)==0) {
#  stop("calc_Q: At least one argument must be supplied: i = k + 1", call.=FALSE)
#}
#i = as.integer(args[1])

