ht_pdf <- function(z, sigma) {
  x = z/sigma
  a = 4/9*exp(-3^2/2)/sqrt(2*pi); c = 7/3
  C = 1/((a/(3-c) - pnorm(-3))*2 + 1)
  p = numeric(length(z))
  p[abs(x) < 3] = C/sqrt(2*pi)*exp(-(x[abs(x) < 3])^2/2)
  p[abs(x) >= 3] = C*a/(abs(x[abs(x) >= 3])-c)^2
  return(p/sigma)
}

get_Q <- function(X_vals, k, sigma) {
  N_Y = 200000;  gamma = 1e-4
  N_X = length(X_vals)
  Y = (-N_Y:N_Y) * gamma
  p <- ht_pdf(Y,sigma)
  log_p <- log(ht_pdf(Y,sigma))
  log_S <- outer(-exp(Y),X_vals) + replicate(N_X, k*Y + log_p)
  log_S <- (log_S - lgamma(k+1))
  log_S <- sweep(log_S, 2, k*log(X_vals), '+')
  S <- exp(log_S)
  return(colSums(S)*gamma)
}

calc_Q_mat <- function(sigma,X_vals, K = 10) {
  N_X = length(X_vals)
  Q_mat <- Matrix(0, nrow= K+1, ncol = N_X)
  batch = 100
  for(i in 1:(K+1)) {
    k = i-1
    print(k)
    for(b in 1:(N_X/batch)) {
      X_ind = (batch*(b-1) + 1):(batch*b)
      curr_X = X_vals[X_ind]
      Q_mat[i, X_ind] = get_Q(curr_X, k, sigma)
    }
  }
  return(Q_mat)
}

calc_Q <- function(x, k) {
  epsilon = 1e-4; X_max = max(X_vals); K = K_val; delta = 1e-5
  x = pmin(pmax(epsilon, x),X_max - epsilon); k = min(k,K)
  l = floor((x/delta)^(2/3))
  prop = (X_vals[l+1] - x)/(X_vals[l+1] - X_vals[l])
  return(prop*Q_mat[k+1,l] + (1-prop)*Q_mat[k+1,l+1])
}

calc_Q_d1 <- function(x, k) {
  epsilon = 1e-4; X_max = max(X_vals); K = K_val - 1; delta = 1e-5
  x = pmin(pmax(epsilon, x),X_max - epsilon); k = min(k,K)
  return(1/x * (-(k+1)*calc_Q(x, k+1) + k*calc_Q(x, k)))
}

calc_Q_d2 <- function(x, k) {
  epsilon = 1e-4; X_max = max(X_vals); K = K_val - 2; delta = 1e-5
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

#return positive semidefinite part
psd <- function(H) {
  eig <- eigen(H); epsilon = 1e-3
  P <- eig$vectors %*% diag(pmax(eig$values,epsilon)) %*% t(eig$vectors)
  return(P)
}
