#' Sets Precomputed Probabiliites as Global Variable
#'
#' Given a matrix, \code{Q_mat}, of P(y|x), under the Poisson-Lognormal model.
#' Sets this as a global variable for fast computations in the future.
#'
#' @param Q_mat_loc Matrix of precomputed probabiliites, as previously computed by \code{\link{get_Q_mat}}
#' @export
set_likelihood_vars <- function(Q_mat_loc) {
  Q_mat <<- Q_mat_loc
  N_X <<- dim(Q_mat)[2]
  delta <<- 1e-5
  X_vals <<- (1:N_X)^1.5*delta
  K_val <<- dim(Q_mat)[1] - 3;
  use_Q <<- T
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

ht_pdf_ds <- function(z, sigma) {
  z = abs(z)
  c = 7/3
  p = ht_pdf(z,sigma)
  S = rep(-1/sigma, length(z))
  g_ind = abs(z) > 3*sigma; l_ind = abs(z) <= 3*sigma
  S[g_ind] = S[g_ind] + 2*z[g_ind]/((z[g_ind]/sigma - c)*sigma^2)
  S[l_ind] = S[l_ind] + z[l_ind]^2/(sigma^3)
  return(S*p)
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

get_Q_d <- function(X_vals, k, sigma, big_params = T) {
  if(big_params) {
    #N_Y = 20000;  gamma = 1e-3
    N_Y = 5000;  gamma = 4e-3
  }
  else {
    N_Y = 5000;  gamma = 4e-3
  }
  #N_Y = 200000;  gamma = 1e-4
  N_X = length(X_vals)
  Y = (-N_Y:N_Y) * gamma
  p <- ht_pdf_ds(Y,sigma)
  pos = p > 0
  log_p <- log(p[pos])
  log_S <- outer(-exp(Y[pos]),X_vals) + replicate(N_X, k*Y[pos] + log_p)
  log_S <- (log_S - lgamma(k+1))
  log_S <- sweep(log_S, 2, k*log(X_vals), '+')
  S <- exp(log_S)
  neg = p < 0
  log_p <- log(-p[neg])
  log_S <- outer(-exp(Y[neg]),X_vals) + replicate(N_X, k*Y[neg] + log_p)
  log_S <- (log_S - lgamma(k+1))
  log_S <- sweep(log_S, 2, k*log(X_vals), '+')
  S_neg <- exp(log_S)
  return((colSums(S) - colSums(S_neg))*gamma)
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

calc_Q_d_one <- function(sigma, X_vals, k, batch = 100, big_params = T) {
  N_X = length(X_vals); results = numeric(N_X)
  print(N_X)
  for(b in 1:ceiling(N_X/batch)) {
    X_ind = (batch*(b-1) + 1):min((batch*b),N_X)
    curr_X = X_vals[X_ind]
    results[X_ind] = get_Q_d(curr_X, k, sigma, big_params = big_params)
  }
  return(results)
}

calc_Q_mat <- function(sigma, X_vals, K = 10, big_params = T) {
  N_X = length(X_vals)
  Q_mat <- Matrix(0, nrow= K+1, ncol = N_X)
  batch = 100
  for(i in 1:(K+1)) {
    k = i-1; print(k)
    Q_mat[i, ] = calc_Q_mat_one(sigma, X_vals, k, batch = batch, big_params = big_params)
  }
  return(Q_mat)
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
  environ = c('calc_Q_mat_one')
  results <- foreach::foreach(i = 1:(K+3), .export = environ) %dopar% {
    cat(paste0("calc_Q: Finished i: ",i,"\n"), file=out_file, append=TRUE)
    k = i-1;
    result = calc_Q_mat_one(sigma, X_vals, k, batch = 100, big_params = big_params)
  }
  parallel::stopCluster(cl)
  return(results)
}

#' Calculates Poisson-lognormal Likelihood
#'
#' For a given value of sigma, calculates P(y|x) under the Poisson-Lognormal model
#' for each possible value of y and x.
#'
#' @param iv Initial Variables: meta data obtained from the \code{\link{init_RCTD}} function
#' @param sigma a numeric representing estimated sigma_c from \code{\link{choose_sigma_c}}
#' @return Returns \code{Q_mat}, a matrix containing pre-computed values for the log likelihood.
#' @export
get_Q_mat <- function(iv, sigma, big_params = T) {
  delta <<- 1e-5
  results <- calc_Q_par(iv$config$K_val, (1:iv$config$N_X)^1.5*delta, sigma, big_params = T)
  Q_mat <-t(as.data.frame(matrix(unlist(results), nrow=length(unlist(results[1])))))
  return(Q_mat)
}

#not using Q also uses small params
calc_Q <- function(x, k) {
  epsilon = 1e-4; X_max = max(X_vals); K = K_val+2; delta = 1e-5
  x = pmin(pmax(epsilon, x),X_max - epsilon); k = min(k,K)
  if(use_Q) {
    l = floor((x/delta)^(2/3))
    prop = (X_vals[l+1] - x)/(X_vals[l+1] - X_vals[l])
    return(prop*Q_mat[k+1,l] + (1-prop)*Q_mat[k+1,l+1])
  }
  else {
    return(calc_Q_mat_one(sigmavar, x, k, batch = 100, big_params = F))
  }
}

#all values of K
calc_Q_all <- function(x, bead) {
  epsilon = 1e-4; X_max = max(X_vals); delta = 1e-5
  x = pmin(pmax(epsilon, x),X_max - epsilon)
  l = floor((x/delta)^(2/3))
  prop = (X_vals[l+1] - x)/(X_vals[l+1] - X_vals[l])
  #N <- length(l)
  #query <- matrix(0,nrow <- 6*N,ncol <- 2)
  #query[1:N,1] <- bead + 1; query[(N+1):(2*N),1] <- bead + 2; query[(2*N+1):(3*N),1] <- bead + 3
  #query[(3*N+1):(4*N),1] <- bead + 1; query[(4*N+1):(5*N),1] <- bead + 2; query[(5*N+1):(6*N),1] <- bead + 3
  #query[1:(3*N),2] <- l; query[(3*N+1):(6*N),2] <- l + 1
  #response <- Q_mat[query]

  #v1 <- response[(3*N+1):(6*N)]
  #k <- response[1:(3*N)] - v1
  #r1 <- k * prop + v1

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
  #return(list(r1[1:N],r1[(N+1):(2*N)],r1[(2*N+1):(3*N)]))
}

calc_Q_d1 <- function(x, k) {
  epsilon = 1e-4; X_max = max(X_vals); K = K_val; delta = 1e-5
  x = pmin(pmax(epsilon, x),X_max - epsilon); k = min(k,K)
  return(1/x * (-(k+1)*calc_Q(x, k+1) + k*calc_Q(x, k)))
}

calc_Q_d2 <- function(x, k) {
  epsilon = 1e-4; X_max = max(X_vals); K = K_val; delta = 1e-5
  x = pmin(pmax(epsilon, x),X_max - epsilon); k = min(k,K)
  sec_der <- 1/(x^2)*((k+1)*(k+2)*calc_Q(x,k+2) - k*(2*(k+1)*calc_Q(x,k+1) - (k-1)*calc_Q(x,k)))
  return(sec_der)
}

calc_log_p <- function(x,k) {
  return(log(calc_Q(x,k)))
}

calc_log_p_d1 <- function(x,k) {
  return(calc_Q_d1(x,k) / calc_Q(x,k))
}

calc_log_p_d2 <- function(x,k) {
  p = calc_Q(x,k)
  result = -calc_Q_d1(x,k)^2/(p^2) + calc_Q_d2(x,k)/p
  return(result)
}

#negative log likelihood
calc_log_l <- function(gene_list, prediction, bead) {
  total_score=0
  for(gene in gene_list) {
    total_score = total_score + calc_log_p(prediction[gene], bead[gene])
  }
  return(-total_score)
}

calc_log_l_par <- function(gene_list, prediction, bead) {
  K = K_val
  bead[bead > K] = K
  total_score=0
  for(k in as.numeric(names(table(bead)))) {
    total_score = total_score + sum(calc_log_p(prediction[bead == k], k))
  }
  return(-total_score)
}

calc_log_d1_par <- function(gene_list, prediction, bead) {
  K = K_val
  bead[bead > K] = K
  results=numeric(length(bead)); names(results) = gene_list
  for(k in as.numeric(names(table(bead)))) {
    results[bead == k] = calc_log_p_d1(prediction[bead == k], k)
  }
  return(-results)
}

calc_log_d2_par <- function(gene_list, prediction, bead) {
  K = K_val
  bead[bead > K] = K
  results=numeric(length(bead)); names(results) = gene_list
  for(k in as.numeric(names(table(bead)))) {
    results[bead == k] = calc_log_p_d2(prediction[bead == k], k)
  }
  return(-results)
}

get_gradient <- function(S, B, gene_list, prediction) {
  d1_vec = calc_log_d1_par(gene_list, prediction, B)
  return(d1_vec %*% S)
}

get_hessian <- function(S, B, gene_list, prediction) {
  d2_vec = calc_log_d2_par(gene_list, prediction, B)
  return(t(S) %*% (diag(d2_vec) %*% S))
}

get_der_fast <- function(S, B, gene_list, prediction) {
  if(use_Q) {
    bead = B; epsilon = 1e-4; X_max = max(X_vals); delta = 1e-5
    x = pmin(pmax(epsilon, prediction),X_max - epsilon);
    Q_cur <- calc_Q_all(x, bead) #need calc_log_d1 and calc_log_d2
    bead[bead > K_val] = K_val
    Q_k <- Q_cur[[1]]
    Q_k1 <- Q_cur[[2]]; Q_k2 <- Q_cur[[3]]
    Q_d1 = 1/x * (-(bead+1)*Q_k1 + bead*Q_k)
    Q_d2 = 1/(x^2)*((bead+1)*(bead+2)*Q_k2 - bead*(2*(bead+1)*Q_k1 - (bead-1)*Q_k))
    d1_vec = as.vector(Q_d1 / Q_k)
    d2_vec = as.vector(-Q_d1^2/(Q_k^2) + Q_d2/Q_k)
    grad = -d1_vec %*% S;
    hess_c <- -d2_vec %*% S_mat
    #hess <- cbind(hess_c[1:2],hess_c[2:3])
    hess <- matrix(0,nrow = dim(S)[2], ncol = dim(S)[2])
    counter = 1
    for(i in 1:dim(S)[2]) {
      l <- dim(S)[2] - i
      hess[i,i:dim(S)[2]] <- hess_c[counter:(counter+l)]
      hess[i,i] <- hess[i,i] / 2
      counter <- counter + l + 1
    }
    hess <- hess + t(hess)
  } else {
    grad = get_gradient(S, B, gene_list, prediction)
    hess = get_hessian(S, B, gene_list, prediction)
  }
  return(list(grad=grad, hess=hess))
}

#return positive semidefinite part
psd <- function(H) {
  eig <- eigen(H); epsilon = 1e-3
  if(length(H) == 1)
    P <- eig$vectors %*% pmax(eig$values,epsilon) %*% t(eig$vectors)
  else
    P <- eig$vectors %*% diag(pmax(eig$values,epsilon)) %*% t(eig$vectors)
  return(P)
}
