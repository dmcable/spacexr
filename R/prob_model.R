#' Sets Precomputed Probabiliites as Global Variable
#'
#' Given a matrix, \code{Q_mat}, of P(y|x), under the Poisson-Lognormal model.
#' Sets this as a global variable for fast computations in the future.
#'
#' @param Q_mat_loc Matrix of precomputed probabiliites, as previously computed by \code{\link{get_Q_mat}}
#' @param X_vals the x-values used for computing the likelihood functions.
#' @export
set_likelihood_vars <- function(Q_mat_loc, X_vals) {
  Q_mat <<- Q_mat_loc
  N_X <<- dim(Q_mat)[2]
  X_vals <<- X_vals
  K_val <<- dim(Q_mat)[1] - 3;
}

set_global_Q_all <- function() {
  Q1 <- readRDS(system.file("extdata", "Qmat/Q_mat_1.rds", package = "spacexr"))
  Q2 <- readRDS(system.file("extdata", "Qmat/Q_mat_2.rds", package = "spacexr"))
  Q_mat_all <<- c(Q1,Q2)
  X_vals <<- readRDS(system.file("extdata", "Qmat/X_vals.rds", package = "spacexr"))
}

ht_pdf <- function(z, sigma) {
  x = z/sigma
  p = ht_pdf_norm(x)
  return(p/sigma)
}

#assumes sigma = 1
ht_pdf_norm <- function(x) {
  a = 4/9*exp(-3^2/2)/sqrt(2*pi); c = 7/3
  C = 1/((a/(3-c) - pnorm(-3))*2 + 1)
  p = numeric(length(x))
  p[abs(x) < 3] = C/sqrt(2*pi)*exp(-(x[abs(x) < 3])^2/2)
  p[abs(x) >= 3] = C*a/(abs(x[abs(x) >= 3])-c)^2
  return(p)
}

get_Q <- function(X_vals, k, sigma, big_params = T) {
  if(big_params) {
    #N_Y = 20000;  gamma = 1e-3
    N_Y = 5000;  gamma = 4e-3
  }
  else {
    N_Y = 5000;  gamma = 4e-3
  }
  N_X = length(X_vals)
  Y = (-N_Y:N_Y) * gamma
  log_p <- log(ht_pdf(Y,sigma))
  log_S <- outer(-exp(Y),X_vals) + replicate(N_X, k*Y + log_p)
  log_S <- (log_S - lgamma(k+1))
  log_S <- sweep(log_S, 2, k*log(X_vals), '+')
  S <- exp(log_S)
  return(colSums(S)*gamma)
}

calc_Q_mat_one <- function(sigma, X_vals, k, batch = 100, big_params = T) {
  N_X = length(X_vals); results = numeric(N_X)
  for(b in 1:ceiling(N_X/batch)) {
    X_ind = (batch*(b-1) + 1):min((batch*b),N_X)
    curr_X = X_vals[X_ind]
    results[X_ind] = get_Q(curr_X, k, sigma, big_params = big_params)
  }
  return(results)
}

calc_Q_par <- function(K, X_vals, sigma, big_params = T) {
  out_file = "logs/calc_Q_log.txt"
  if (file.exists(out_file))
    file.remove(out_file)
  numCores = parallel::detectCores(); MAX_CORES = 8
  if(parallel::detectCores() > MAX_CORES)
    numCores <- MAX_CORES
  cl <- parallel::makeCluster(numCores,outfile="") #makeForkCluster
  doParallel::registerDoParallel(cl)
  environ = c('calc_Q_mat_one','get_Q','ht_pdf','ht_pdf_norm')
  results <- foreach::foreach(i = 1:(K+3), .export = environ) %dopar% {
    cat(paste0("calc_Q: Finished i: ",i,"\n"), file=out_file, append=TRUE)
    k = i-1;
    result = calc_Q_mat_one(sigma, X_vals, k, batch = 100, big_params = big_params)
  }
  parallel::stopCluster(cl)
  return(results)
}

#all values of K
calc_Q_all <- function(x, bead) {
  epsilon = 1e-4; X_max = max(X_vals); delta = 1e-5
  x = pmin(pmax(epsilon, x),X_max - epsilon)
  l = floor((x/delta)^(2/3))
  l <- pmin(l, 900) + floor(pmax(l - 900, 0)/30)
  prop = (X_vals[l+1] - x)/(X_vals[l+1] - X_vals[l])
  v1 <- Q_mat[cbind(bead+1,l+1)]
  k <- Q_mat[cbind(bead+1,l)] - v1
  r1 <- k * prop + v1
  v1 <- Q_mat[cbind(bead+2,l+1)]
  k <- Q_mat[cbind(bead+2,l)] - v1
  r2 <- k * prop + v1
  v1 <- Q_mat[cbind(bead+3,l+1)]
  k <- Q_mat[cbind(bead+3,l)] - v1
  r3 <- k * prop + v1
  return(list(r1,r2,r3))
}

#just one value of k
calc_Q_k <- function(x, bead) {
  bead[bead > K_val] = K_val
  epsilon = 1e-4; X_max = max(X_vals); delta = 1e-5
  x = pmin(pmax(epsilon, x),X_max - epsilon)
  l = floor((x/delta)^(2/3))
  l <- pmin(l, 900) + floor(pmax(l - 900, 0)/30)
  prop = (X_vals[l+1] - x)/(X_vals[l+1] - X_vals[l])
  v1 <- Q_mat[cbind(bead+1,l+1)]
  k <- Q_mat[cbind(bead+1,l)] - v1
  r1 <- k * prop + v1
  return(r1)
}

#negative log likelihood
calc_log_l_vec <- function(lambda, Y, return_vec = FALSE) {
  log_l_vec <- -log(calc_Q_k(lambda,Y))
  if(return_vec)
    return(log_l_vec)
  return(sum(log_l_vec))
}


get_d1_d2 <- function(B, prediction) {
  bead = B; epsilon = 1e-4; X_max = max(X_vals);
  x = pmin(pmax(epsilon, prediction),X_max - epsilon);
  Q_cur <- calc_Q_all(x, bead) #need calc_log_d1 and calc_log_d2
  bead[bead > K_val] = K_val
  Q_k <- Q_cur[[1]]
  Q_k1 <- Q_cur[[2]]; Q_k2 <- Q_cur[[3]]
  Q_d1 = 1/x * (-(bead+1)*Q_k1 + bead*Q_k)
  Q_d2 = 1/(x^2)*((bead+1)*(bead+2)*Q_k2 - bead*(2*(bead+1)*Q_k1 - (bead-1)*Q_k))
  d1_vec = as.vector(Q_d1 / Q_k)
  d2_vec = as.vector(-Q_d1^2/(Q_k^2) + Q_d2/Q_k)
  return(list(d1_vec = d1_vec, d2_vec = d2_vec))
}

get_der_fast <- function(S, B, gene_list, prediction, bulk_mode = F) {
  if(bulk_mode) {
    #d1_vec <- -t(log(prediction) - log(B))
    #d2_vec <- -t(1/prediction)
    d1_vec <- -2*t((log(prediction) - log(B))/prediction)
    d2_vec <- -2*t((1 - log(prediction) + log(B))/prediction^2)
  } else {
    d1_d2 <- get_d1_d2(B, prediction)
    d1_vec <- d1_d2$d1_vec
    d2_vec <- d1_d2$d2_vec
  }
  grad = -d1_vec %*% S;
  hess_c <- -d2_vec %*% S_mat
  hess <- matrix(0,nrow = dim(S)[2], ncol = dim(S)[2])
  counter = 1
  for(i in 1:dim(S)[2]) {
    l <- dim(S)[2] - i
    hess[i,i:dim(S)[2]] <- hess_c[counter:(counter+l)]
    hess[i,i] <- hess[i,i] / 2
    counter <- counter + l + 1
  }
  hess <- hess + t(hess)
  return(list(grad=grad, hess=hess))
}

#return positive semidefinite part
psd <- function(H, epsilon = 1e-3) {
  eig <- eigen(H);
  if(length(H) == 1)
    P <- eig$vectors %*% pmax(eig$values,epsilon) %*% t(eig$vectors)
  else
    P <- eig$vectors %*% diag(pmax(eig$values,epsilon)) %*% t(eig$vectors)
  return(P)
}
