#' Sets Precomputed Probabiliites as Global Variable
#'
#' Given a matrix, \code{Q_mat}, or log P(y|x), under the Poisson-Lognormal model.
#' Sets this as a global variable for fast computations in the future.
#'
#' @param Q_mat_loc Matrix of precomputed probabiliites, as previously computed by \code{\link{get_Q_mat}}
#' @param X_vals the x-values used for computing the likelihood functions.
#' @param sigma (default NULL). If NULL, computes SQ_mat according to Q_mat_loc.
#' Else, uses precomputed values of SQ_mat stored in SQ_mat_all with index sigma
#' @export
set_likelihood_vars <- function(Q_mat_loc, X_vals, sigma = NULL) {
  Q_mat <<- Q_mat_loc
  N_X <<- dim(Q_mat)[2]
  X_vals <<- X_vals
  K_val <<- dim(Q_mat)[1] - 3;
  if(is.null(sigma))
    SQ_mat <<- solve_sq(Q_mat_loc, X_vals)
  else
    SQ_mat <<- SQ_mat_all[[sigma]]
}

solve_sq <- function(Q_mat, X_vals) {
  # solve for the s coefficients
  n <- dim(Q_mat)[2] - 1
  M <- matrix(0, n-1, n-1)
  del <- diff(X_vals)
  diag(M) <- 2*(del[1:(n-1)] + del[2:n])
  M[cbind(2:(n-1), 1:(n-2))] <- del[2:(n-1)]
  M[cbind(1:(n-2), 2:(n-1))] <- del[2:(n-1)]
  MI <- solve(M)
  fB <- sweep(diff(t(Q_mat)),1,del,'/')
  fBD <- 6*diff(fB)
  SQ_mat <- t(MI %*% fBD)
  ### append the zeros
  SQ_mat <- cbind(0, SQ_mat, 0)
}

set_global_Q_all <- function() {
  message('set_global_Q_all: begin')
  Q_mat_all <<- get_Q_all()
  X_vals <<- readRDS(system.file("extdata", "Qmat/X_vals.rds", package = "spacexr"))
  SQ_mat_all <<- lapply(Q_mat_all, function(x) solve_sq(x, X_vals))
  message('set_global_Q_all: finished')
}

set_likelihood_vars_sigma <- function(sigma) {
  Q_mat_all <- get_Q_all()
  X_vals <- readRDS(system.file("extdata", "Qmat/X_vals.rds", package = "spacexr"))
  set_likelihood_vars(Q_mat_all[[sigma]], X_vals)
}

get_Q_all <- function() {
  Q1 <- readRDS(system.file("extdata", "Qmat/Q_mat_1.rds", package = "spacexr"))
  Q2 <- readRDS(system.file("extdata", "Qmat/Q_mat_2.rds", package = "spacexr"))
  Q3 <- readRDS(system.file("extdata", "Qmat/Q_mat_3.rds", package = "spacexr"))
  Q4 <- readRDS(system.file("extdata", "Qmat/Q_mat_4.rds", package = "spacexr"))
  Q5 <- readRDS(system.file("extdata", "Qmat/Q_mat_5.rds", package = "spacexr"))
  Q_mat_all <- c(Q1, Q2, Q3, Q4, Q5)
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
calc_Q_all <- function(Y, lambda) {
  Y[Y > K_val] <- K_val
  epsilon <- 1e-4; X_max <- max(X_vals); delta <- 1e-6
  lambda <- pmin(pmax(epsilon, lambda),X_max - epsilon)

  l <- floor((lambda/delta)^(1/2))
  m <- pmin(l - 9,40) + pmax(ceiling(sqrt(pmax(l-48.7499,0)*4))-2,0)
  ti1 <- X_vals[m]; ti <- X_vals[m+1]; hi <- ti - ti1
  #Q0 <- cbind(Y+1, m); Q1 <- cbind(Y+1, m+1)
  Q0 <- cbind(Y+1, m); Q1 <- Q0; Q1[,2] <- Q1[,2] + 1
  fti1 <- Q_mat[Q0]; fti <- Q_mat[Q1]
  zi1 <- SQ_mat[Q0]; zi <- SQ_mat[Q1]
  diff1 <- lambda - ti1; diff2 <- ti - lambda
  diff3 <- fti/hi-zi*hi/6; diff4 <- fti1/hi-zi1*hi/6
  zdi <- zi / hi; zdi1 <- zi1 / hi
  # cubic spline interpolation
  d0_vec <- zdi*(diff1)^3/6 + zdi1*(diff2)^3/6 + diff3*diff1 + diff4*diff2
  d1_vec <- zdi*(diff1)^2/2 - zdi1*(diff2)^2/2 + diff3 - diff4
  d2_vec <- zdi*(diff1) + zdi1*(diff2)
  return(list(d0_vec = d0_vec, d1_vec = d1_vec, d2_vec = d2_vec))
}

#negative log likelihood
calc_log_l_vec <- function(lambda, Y, return_vec = FALSE) {
  log_l_vec <- -calc_Q_all(Y,lambda)$d0_vec
  if(return_vec)
    return(log_l_vec)
  return(sum(log_l_vec))
}

# linear interpolation
calc_log_l_vec_fast <- function(lambda, Y) {
  Y[Y > K_val] <- K_val
  epsilon <- 1e-4; X_max <- max(X_vals); delta <- 1e-6
  lambda <- pmin(pmax(epsilon, lambda),X_max - epsilon)
  l <- floor((lambda/delta)^(1/2))
  m <- pmin(l - 9,40) + pmax(ceiling(sqrt(pmax(l-48.7499,0)*4))-2,0)
  Q0 <- cbind(Y+1, m); Q1 <- Q0; Q1[,2] <- Q1[,2] + 1
  fti1 <- Q_mat[Q0]; fti <- Q_mat[Q1]
  prop = (X_vals[m+1] - lambda)/(X_vals[m+1] - X_vals[m])
  r1 <- prop * fti1 + (1 - prop) * fti
  return(-sum(r1))
}



get_d1_d2 <- function(Y, lambda) {
  d_all <- calc_Q_all(Y, lambda)
  return(list(d1_vec = d_all$d1_vec, d2_vec = d_all$d2_vec))
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
