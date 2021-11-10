
choose_sigma_gene <- function(sigma_init, Y, X1, X2, my_beta, nUMI,test_mode, verbose = F, n.iter = 100, MIN_CHANGE = 0.001, MAX_ITER_SIGMA = 10) {
  sigma_s_best <- sigma_init
  sigma_vals <- names(Q_mat_all)
  for(iter in 1:MAX_ITER_SIGMA) {
    last_sigma <- sigma_s_best
    set_likelihood_vars(Q_mat_all[[as.character(sigma_s_best)]], X_vals)
    res <- estimate_effects_trust(Y,X1,X2,my_beta, nUMI,test_mode, verbose = verbose, n.iter = n.iter, MIN_CHANGE = MIN_CHANGE)
    prediction <- res$prediction
    pred_c <- as.vector(prediction * exp((as.numeric(sigma_s_best) / 100)^2/2))
    res_val <- numeric(length(sigma_vals))
    names(res_val) <- sigma_vals
    for(sigma_s in sigma_vals) {
      set_likelihood_vars(Q_mat_all[[as.character(sigma_s)]], X_vals)
      res_val[as.character(sigma_s)] <- (calc_log_l_vec(pred_c / exp((as.numeric(sigma_s) / 100)^2/2), as.vector(t(Y))))
    }
    sigma_s_best <- names(which.min(res_val))
    if(sigma_s_best == last_sigma) {
      break
    }
  }
  return(list(sigma_s_best = sigma_s_best, res = res))
}


construct_hess_fast <- function(X1,X2,lambda,lambda_k, K, d1_d2, d1_lam) {
  g_1 <- sweep(X1, 1,lambda,'*')
  g_2 <- matrix(0, nrow = dim(X2)[1], K*dim(X2)[2])
  for(k in 1:K)
    g_2[,(1 + dim(X2)[2]*(k-1)):(k*dim(X2)[2])] <- sweep(X2, 1,lambda_k[,k],'*')
  grad <- cbind(g_1, g_2)
  grad_Q <- sweep(grad, 1, d1_d2$d2_vec, '*')
  H1 <- -t(grad_Q) %*% grad
  X1_Q <- sweep(X1,1, lambda*d1_d2$d1_vec, '*')
  H2_11 <- t(X1_Q) %*% X1
  L1 <- dim(X1)[2]; L2 <- dim(X2)[2]
  H2_12 <- matrix(0,nrow = L1, ncol = L2*K)
  for(k in 1:K) {
    X1_Q <- sweep(X1,1, d1_lam[,k], '*')
    H2_12[,(1 + L2*(k-1)):(L2*k)] <- t(X1_Q) %*% X2
  }
  H2 <- matrix(0, nrow = L1 + L2*K, ncol = L1 + L2*K)
  if(L1 > 0) {
    H2[(1:L1),1:L1] <- H2_11
    H2[1:L1,(L1+1):(L1+L2*K)] <- H2_12
    H2[(L1+1):(L1+L2*K),1:L1] <- t(H2_12)
  }
  for(k in 1:K) {
    X2_Q <- sweep(X2,1, d1_lam[,k], '*')
    H2[(L1 + 1 + (k-1)*L2):(L1 + k*L2), (L1 + 1 + (k-1)*L2):(L1 + k*L2)] <- t(X2_Q) %*% X2
  }
  H <- (H1 - H2)
}

solveIRWLS.effects_trust <-function(Y, X1, X2, my_beta, test_mode, verbose = FALSE, n.iter = 200, MIN_CHANGE = .001, PRECISION.THRESHOLD = .01){
  L1 <- dim(X1)[2]; L2 <- dim(X2)[2]
  n_cell_types = dim(my_beta)[2]
  alpha1 <- numeric(dim(X1)[2])
  delta <- 0.1 #trust region radius
  alpha2 <- matrix(0,nrow = dim(X2)[2], ncol = n_cell_types) # initialize it to be the previous cell type means
  alpha2[1,] <- -5
  if(test_mode == 'multi') {
    alpha2[,] <- -5 # multi mode
  }
  pred_decrease_vals <- numeric(n.iter)
  Y[Y > K_val] <- K_val
  K <- dim(my_beta)[2]
  lambda_k <- exp(sweep(X2 %*% (alpha2), 1, X1 %*% (alpha1),'+'))*my_beta #J by K
  lambda <- rowSums(lambda_k)
  prev_ll <- calc_log_l_vec(lambda,Y)
  error_vec <- (1:dim(my_beta)[2]) == 0
  for(itera in 1:n.iter) {
    threshold = 1e-8
    lambda[lambda < threshold] <- threshold
    d1_d2 <- get_d1_d2(Y, lambda)
    grad_1 <- d1_d2$d1_vec %*% sweep(X1, 1,lambda,'*')
    d1_lam <- sweep(lambda_k, 1, d1_d2$d1_vec, '*')
    grad_2 <- t(X2) %*% d1_lam
    H <- construct_hess_fast(X1,X2,lambda,lambda_k, K, d1_d2, d1_lam)
    d_vec_o <- c(grad_1, grad_2)
    D_mat_o <- psd(H)
    norm_factor <- norm(D_mat_o,"2")
    D_mat <- D_mat_o / norm_factor
    d_vec <- d_vec_o / norm_factor
    epsilon <- 1e-7; D_mat <- D_mat + epsilon * diag(length(d_vec))
    A <- cbind(diag(dim(D_mat)[2]), -diag(dim(D_mat)[2]))
    bzero <- rep(-delta,2*dim(D_mat)[2])
    #print(max(abs(bzero)))
    solution <-  quadprog::solve.QP(D_mat,d_vec,A,bzero,meq=0)$solution
    predicted_decrease = -(0.5*t(solution) %*% D_mat_o %*% solution - sum(d_vec_o*solution))

    alpha1_new <- alpha1 + solution[1:L1]
    alpha2_new <- alpha2 + matrix(solution[(L1 + 1):length(solution)], nrow = L2, ncol = K)
    lambda_k_new <- exp(sweep(X2 %*% (alpha2_new), 1, X1 %*% (alpha1_new),'+'))*my_beta #J by K
    error_vec <- is.na(colMeans(lambda_k_new)) | error_vec
    lambda_k_new[is.na(lambda_k_new)] <- 1
    lambda_new <- rowSums(lambda_k_new)
    log_l <- calc_log_l_vec(lambda_new,Y)

    true_decrease = prev_ll - log_l
    pred_decrease_vals[itera] <- predicted_decrease
    if(true_decrease >= 0.1 * predicted_decrease) {
      delta = 1.1*delta
      alpha1 <- alpha1_new
      alpha2 <- alpha2_new
      lambda_k <- lambda_k_new
      lambda <- lambda_new
      prev_ll <- log_l
    } else {
      delta = min(1,0.5*delta)
    }
    if(delta < MIN_CHANGE || (itera > 5 && max(pred_decrease_vals[(itera - 5):itera]) < min(MIN_CHANGE,1e-5))) {
      break
    }
  }

  threshold = 1e-8
  lambda[lambda < threshold] <- threshold
  d1_d2 <- get_d1_d2(Y, lambda)
  d1_lam <- sweep(lambda_k, 1, d1_d2$d1_vec, '*')
  H <- construct_hess_fast(X1,X2,lambda,lambda_k, K, d1_d2, d1_lam)
  eps = 1e-8
  if(min(eigen(H)$values) < eps) {
    I <- solve(H + diag(rep(eps,dim(H)[1])))
  } else {
    I <- solve(H)
  }
  precision <- abs(solve(D_mat) %*% d_vec)
  cell_type_ind <- get_cell_type_ind(X1,X2, dim(my_beta)[2])
  converged_vec <- (1:dim(my_beta)[2]) == 0
  if(itera < n.iter) {
    converged_vec <- !converged_vec
    converged_vec[unique(cell_type_ind[precision > PRECISION.THRESHOLD])] <- FALSE
  }
  converged_vec <- converged_vec & (!error_vec)
  names(converged_vec) <- colnames(my_beta)
  names(error_vec) <- colnames(my_beta)
  return(list(alpha1 = alpha1, alpha2 = alpha2, converged = any(converged_vec), I = I, d = d_vec_o,
              n.iter = itera, log_l = prev_ll, precision = precision, prediction = lambda, 
              converged_vec = converged_vec, error_vec = error_vec))
}

estimate_gene_wrapper <- function(Y,X1,X2,my_beta, nUMI, sigma_init, test_mode, verbose = F, n.iter = 200, MIN_CHANGE = 1e-3, sigma_gene = T, PRECISION.THRESHOLD = 0.01) {
  if(sigma_gene)
    return(choose_sigma_gene(sigma_init, Y, X1, X2, my_beta, nUMI,test_mode, verbose = verbose, n.iter = n.iter, MIN_CHANGE = MIN_CHANGE))
  else {
    res <- estimate_effects_trust(Y,X1,X2,my_beta, nUMI,test_mode, verbose = verbose, n.iter = n.iter, MIN_CHANGE = MIN_CHANGE, PRECISION.THRESHOLD = PRECISION.THRESHOLD)
    return(list(sigma_s_best = -1, res = res))
  }
}

estimate_effects_trust <- function(Y, X1, X2, my_beta, nUMI, test_mode, verbose = F, n.iter = 200, MIN_CHANGE = 1e-3, PRECISION.THRESHOLD = 0.01) {
  my_beta<- sweep(my_beta,1, nUMI, '*')
  solveIRWLS.effects_trust(Y,X1,X2, my_beta, test_mode, verbose = verbose, n.iter = n.iter, MIN_CHANGE = MIN_CHANGE, PRECISION.THRESHOLD = PRECISION.THRESHOLD)
}
