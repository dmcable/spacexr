library(Matrix)
library(RCTD)

N_X = 20000; delta = 1e-5; sigma = 0.8
X = (1:N_X)^1.5*delta
QQ_mat <- calc_Q_mat(sigma, X, K = 30)
saveRDS(QQ_mat, file.path("Data/Slideseq/Puck_Viz/results","QQ_mat.RDS"))
