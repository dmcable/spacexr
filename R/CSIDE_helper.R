
choose_sigma_gene <- function(sigma_init, Y, X1, X2, my_beta, nUMI,test_mode, verbose = F, n.iter = 100, MIN_CHANGE = 0.001, MAX_ITER_SIGMA = 10, PRECISION.THRESHOLD = .01) {
  sigma_s_best <- sigma_init
  sigma_vals <- names(Q_mat_all)
  alpha1 <- NULL; alpha2 <- NULL;
  for(iter in 1:MAX_ITER_SIGMA) {
    last_sigma <- sigma_s_best
    set_likelihood_vars(Q_mat_all[[as.character(sigma_s_best)]], X_vals)
    res <- estimate_effects_trust(Y,X1,X2,my_beta, nUMI,test_mode, verbose = verbose, n.iter = n.iter,
                                  MIN_CHANGE = MIN_CHANGE, PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                                  alpha1_init = alpha1, alpha2_init = alpha2)
    alpha1 <- res$alpha1; alpha2 <- res$alpha2
    prediction <- res$prediction
    pred_c <- as.vector(prediction)
    res_val <- numeric(length(sigma_vals))
    names(res_val) <- sigma_vals
    for(sigma_s in sigma_vals) {
      set_likelihood_vars(Q_mat_all[[as.character(sigma_s)]], X_vals)
      res_val[as.character(sigma_s)] <- (calc_log_l_vec(pred_c, as.vector(t(Y))))
    }
    sigma_s_best <- names(which.min(res_val))
    if(sigma_s_best == last_sigma) {
      break
    }
  }
  return(list(sigma_s_best = sigma_s_best, res = res))
}

mysweep <- function(X1,tX2,lambda,lambda_k, K, d1_d2, d1_lam) {
  g_2 <- matrix(0, nrow = K*dim(tX2)[1], dim(tX2)[2])
  if(TRUE) {
  for(k in 1:K) {
    #g_2[,(1 + dim(X2)[2]*(k-1)):(k*dim(X2)[2])] <- sweep(X2, 1,lambda_k[,k],'*')
    mymat <- Rfast::eachrow(tX2, lambda_k[,k], oper = "*")
    g_2[(1 + dim(tX2)[1]*(k-1)):(k*dim(tX2)[1]),] <- mymat

  }
  }
  return(t(g_2))
}

sweep1 <- function(tX1, lambda) {
  #g_1 <- sweep(X1, 1,lambda,'*')
  g_1 <- Rfast::eachrow(tX1, lambda,oper = "*")
  return(t(g_1))
}

sweep2 <- function(tX1, d1_lam, k) {
  #X1_Q <- sweep(X1,1, d1_lam[,k], '*')
  if(dim(tX1)[1] > 0) {
    X1_Q <- t(Rfast::eachrow(tX1, d1_lam[,k],oper = "*"))
  } else {
    X1_Q <- t(tX1)
  }
  return(X1_Q)
}

sweep3 <- function(tX2, d1_lam, k) {
  #X2_Q <- sweep(X2,1, d1_lam[,k], '*')
  X2_Q <- Rfast::eachrow(tX2, d1_lam[,k],oper = "*")
  return(t(X2_Q))
}


mysweept <- function(tX2,tlk, K) {
  #g_2 <- matrix(0, nrow = K*dim(tX2)[1], dim(tX2)[2])
  if(F) {
    for(k in 1:K) {
      #g_2[,(1 + dim(X2)[2]*(k-1)):(k*dim(X2)[2])] <- sweep(X2, 1,lambda_k[,k],'*')
      mymat <- Rfast::eachrow(tX2, lambda_k[,k], oper = "*")
      g_2[(1 + dim(tX2)[1]*(k-1)):(k*dim(tX2)[1]),] <- mymat

    }
  } else {
    g_2 <- tX2[rep(1:dim(tX2)[1],K),] * tlk[rep(1:K, each = dim(tX2)[1]), ]
  }
  return(g_2)
}

sweep1t <- function(tX1, lambda) {
  #g_1 <- sweep(X1, 1,lambda,'*')
  g_1 <- Rfast::eachrow(tX1, lambda,oper = "*")
  return(g_1)
}

sweep2t <- function(tX1, tdl, k) {
  #X1_Q <- sweep(X1,1, d1_lam[,k], '*')
  if(dim(tX1)[1] > 0) {
    X1_Q <- Rfast::eachrow(tX1, tdl[k,],oper = "*")
  } else {
    X1_Q <- tX1
  }
  return(X1_Q)
}

sweep3t <- function(tX2, tdl, k) {
  #X2_Q <- sweep(X2,1, d1_lam[,k], '*')
  X2_Q <- Rfast::eachrow(tX2, tdl[k, ],oper = "*")
  return(X2_Q)
}

sweep3t_all <- function(tX2, tdl, K) {
  #X2_Q <- sweep(X2,1, d1_lam[,k], '*')
  X2_Q <- tX2[rep(1:dim(tX2)[1],K),] * tdl[rep(1:K, each = dim(tX2)[1]), ]
  return(X2_Q)
}

#grad_1 <- d1_d2$d1_vec %*% sweep(X1, 1,lambda,'*')
#d1_lam <- sweep(lambda_k, 1, d1_d2$d1_vec, '*')
#grad_2 <- t(X2) %*% d1_lam

#grad_1 <- d1_d2$d1_vec %*% g_1
#d1_lam <- sweep(lambda_k, 1, d1_d2$d1_vec, '*')
#grad_2 <- t(X2) %*% d1_lam

construct_hess_fast_old <- function(X1,X2,lambda,lambda_k, K, d1_d2, d1_lam) {
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

construct_hess_fast_nont <- function(X1,X2,lambda,lambda_k, K, d1_d2, d1_lam) {
  tX1 <- t(X1); tX2 <- t(X2)
  g_1 <- sweep1(tX1, lambda)
  #g_1 <- sweep(X1, 1,lambda,'*')
  g_2 <- mysweep(X1,tX2,lambda,lambda_k, K, d1_d2, d1_lam)
  #g_2 <- matrix(0, nrow = dim(X2)[1], K*dim(X2)[2])
  #for(k in 1:K)
  #  g_2[,(1 + dim(X2)[2]*(k-1)):(k*dim(X2)[2])] <- sweep(X2, 1,lambda_k[,k],'*')
  grad <- cbind(g_1, g_2)
  #grad_Q <- sweep(grad, 1, d1_d2$d2_vec, '*')
  grad_Q <- t(Rfast::eachrow(t(grad), d1_d2$d2_vec, '*'))
  H1 <- -t(grad_Q) %*% grad
  #X1_Q <- sweep(X1,1, lambda*d1_d2$d1_vec, '*')
  X1_Q <- t(Rfast::eachrow(tX1, lambda*d1_d2$d1_vec, '*'))
  H2_11 <- t(X1_Q) %*% X1
  L1 <- dim(X1)[2]; L2 <- dim(X2)[2]
  H2_12 <- matrix(0,nrow = L1, ncol = L2*K)
  for(k in 1:K) {
    # X1_Q <- sweep(X1,1, d1_lam[,k], '*')
    X1_Q <- sweep2(tX1, d1_lam, k)
    H2_12[,(1 + L2*(k-1)):(L2*k)] <- t(X1_Q) %*% X2
  }
  H2 <- matrix(0, nrow = L1 + L2*K, ncol = L1 + L2*K)
  if(L1 > 0) {
    H2[(1:L1),1:L1] <- H2_11
    H2[1:L1,(L1+1):(L1+L2*K)] <- H2_12
    H2[(L1+1):(L1+L2*K),1:L1] <- t(H2_12)
  }
  for(k in 1:K) {
    #X2_Q <- sweep(X2,1, d1_lam[,k], '*')
    X2_Q <- sweep3(tX2, d1_lam, k)
    H2[(L1 + 1 + (k-1)*L2):(L1 + k*L2), (L1 + 1 + (k-1)*L2):(L1 + k*L2)] <- t(X2_Q) %*% X2
  }
  H <- (H1 - H2)
}

#grad_1 <- d1_d2$d1_vec %*% sweep(X1, 1,lambda,'*')
#d1_lam <- sweep(lambda_k, 1, d1_d2$d1_vec, '*')
#grad_2 <- t(X2) %*% d1_lam

construct_hess_fast <- function(X1,X2,lambda,lambda_k, K, d1_d2) {
  tX1 <- t(X1); tX2 <- t(X2)
  g_1 <- sweep1t(tX1, lambda)
  tlk <- t(lambda_k)
  tdl <- Rfast::eachrow(tlk, d1_d2$d1_vec, '*')
  #d1_lam <- t(tdl)
  g_2 <- mysweept(tX2,tlk, K)
  grad <- rbind(g_1, g_2)
  grad_Q <- Rfast::eachrow(grad, d1_d2$d2_vec, '*')
  H1 <- -grad_Q %*% t(grad)
  X1_Q <- Rfast::eachrow(tX1, lambda*d1_d2$d1_vec, '*')
  #d1_lam <- sweep(lambda_k, 1, d1_d2$d1_vec, '*')
  grad_1 <- matrix(rowSums(X1_Q), 1, dim(X1)[2])
  H2_11 <- X1_Q %*% X1
  L1 <- dim(X1)[2]; L2 <- dim(X2)[2]
  H2_12 <- matrix(0,nrow = L1, ncol = L2*K)
  for(k in 1:K) {
    X1_Q <- sweep2t(tX1, tdl, k)
    H2_12[,(1 + L2*(k-1)):(L2*k)] <- X1_Q %*% X2
  }
  H2 <- matrix(0, nrow = L1 + L2*K, ncol = L1 + L2*K)
  if(L1 > 0) {
    H2[(1:L1),1:L1] <- H2_11
    H2[1:L1,(L1+1):(L1+L2*K)] <- H2_12
    H2[(L1+1):(L1+L2*K),1:L1] <- t(H2_12)
  }
  if(QUORP) {
    X2_Q <- sweep3t_all(tX2, tdl, K)
    grad_2 <- matrix(rowSums(X2_Q),dim(X2)[2],K)
    H2m <- X2_Q %*% X2
    for(k in 1:K) {
      #X2_Q <- sweep3t(tX2, d1_lam, k)
      H2[(L1 + 1 + (k-1)*L2):(L1 + k*L2), (L1 + 1 + (k-1)*L2):(L1 + k*L2)] <-
        H2m[(1 + (k-1)*L2):(k * L2),]# X2_Q %*% X2
    }
  } else {
    for(k in 1:K) {
      X2_Q <- sweep3t(tX2, tdl, k)
      H2[(L1 + 1 + (k-1)*L2):(L1 + k*L2), (L1 + 1 + (k-1)*L2):(L1 + k*L2)] <-
        X2_Q %*% X2
    }
  }
  H <- (H1 - H2)
  return(list(H = H, grad_1 = grad_1, grad_2 = grad_2))
}

solveIRWLS.effects_trust <-function(Y, X1, X2, my_beta, test_mode, verbose = FALSE,
                                    n.iter = 200, MIN_CHANGE = .01, PRECISION.THRESHOLD = .05,
                                    alpha1_init = NULL, alpha2_init = NULL){
  lam_threshold = 1e-8; MIN_ITERATIONS <- 15
  beta_succ <- 1.1; beta_fail <- 0.5;
  gamma <- 0.25 # prev gamma = 0.1
  epsilon_2 <- 1e-6 * length(Y);
  #epsilon_2 <- 1e-10;
  delta <- 0.1 #trust region radius
  init_val <- min(-5, log(10/median(rowSums(my_beta))))
  L1 <- dim(X1)[2]; L2 <- dim(X2)[2]
  n_cell_types = dim(my_beta)[2]
  if(is.null(alpha1_init))
    alpha1 <- numeric(dim(X1)[2])
  else
    alpha1 <- alpha1_init
  if(is.null(alpha2_init)) {
    alpha2 <- matrix(0,nrow = dim(X2)[2], ncol = n_cell_types) # initialize it to be the previous cell type means
    alpha2[1,] <- init_val
    if(test_mode == 'categorical') {
      alpha2[,] <- init_val # multi mode
    }
  } else {
    alpha2 <- alpha2_init
  }
  pred_decrease_vals <- numeric(n.iter)
  Y[Y > K_val] <- K_val
  K <- dim(my_beta)[2]
  lambda_k <- exp(sweep(X2 %*% (alpha2), 1, X1 %*% (alpha1),'+'))*my_beta #J by K
  lambda <- rowSums(lambda_k)
  prev_ll <- calc_log_l_vec(lambda,Y)
  error_vec <- (1:dim(my_beta)[2]) == 0
  for(itera in 1:n.iter) {
    lambda[lambda < lam_threshold] <- lam_threshold
    d1_d2 <- get_d1_d2(Y, lambda)
    #grad_1 <- d1_d2$d1_vec %*% sweep(X1, 1,lambda,'*')
    #d1_lam <- sweep(lambda_k, 1, d1_d2$d1_vec, '*')
    #grad_2 <- t(X2) %*% d1_lam
    H_list <- construct_hess_fast(X1,X2,lambda,lambda_k, K, d1_d2)
    H <- H_list$H
    grad_1 <- H_list$grad_1; grad_2 <- H_list$grad_2
    d_vec_o <- c(grad_1, grad_2)
    D_mat_o <- psd(H)
    norm_factor <- norm(D_mat_o,"2")
    D_mat <- D_mat_o / norm_factor
    d_vec <- d_vec_o / norm_factor
    epsilon <- 1e-7; D_mat <- D_mat + epsilon * diag(length(d_vec))
    A <- cbind(diag(dim(D_mat)[2]), -diag(dim(D_mat)[2]))
    bzero <- rep(-delta,2*dim(D_mat)[2])
    D_mat <- D_mat + diag(dim(D_mat)[1])*1e-10 # avoid numerical errors
    solution <-  quadprog::solve.QP(D_mat,d_vec,A,bzero,meq=0)$solution #* 0.01 # CHANGE THIS
    ### SCRATCH
      if(FALSE) {
      ind <- 26
      solution_new <- solution*0
      solution_new[ind] <- solution[ind]
      solution <- solution_new

      log_l_v <- calc_log_l_vec(lambda_new,Y, return_vec = T)
      plog_l_v <- calc_log_l_vec(lambda,Y, return_vec = T)
      diff_v <- log_l_v - plog_l_v
      print(sum(log_l_v) - sum(plog_l_v))
      print(sum(plog_l_v) - prev_ll)
      print(sum(log_l_v) - log_l)
      print(prev_ll - log_l)
      pd_v <- -X2[,2] * d1_lam[,13] * solution[ind]
      d1_lam[wind,13]
      error <- pd_v - diff_v
      pd_v[which(error != 0)[200]]
      diff_v[which(error != 0)[200]]
      pd_v[wind]
      diff_v[wind]
      errorn <- error[error != 0]
      tail(errorn[order(abs(errorn))])
      wind <- which.max(abs(error))
      print(abs(true_decrease - predicted_decrease))
      print(abs(true_decrease - predicted_decrease) / predicted_decrease)
    }
    ### END SCRATCH
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
    if(abs(true_decrease - predicted_decrease) <= (1-gamma) * predicted_decrease) {
      delta = beta_succ*delta
      alpha1 <- alpha1_new
      alpha2 <- alpha2_new
      lambda_k <- lambda_k_new
      lambda <- lambda_new
      prev_ll <- log_l
    } else {
      delta = min(1,beta_fail*delta)
    }
    if(F && itera >= MIN_ITERATIONS) {
      print(itera)
      print(delta)
      print(max(pred_decrease_vals[(itera - MIN_ITERATIONS+1):itera]))
    }
    if(delta < MIN_CHANGE || (itera >= MIN_ITERATIONS &&
       max(pred_decrease_vals[(itera - MIN_ITERATIONS+1):itera]) < min(epsilon_2))) {
      break
    }
    if(F) { # SCRATCH TEMP
      precision <- abs(solve(D_mat) %*% d_vec)
      if(itera >= MIN_ITERATIONS)
        print(paste(precision[28],max(pred_decrease_vals[(itera - MIN_ITERATIONS+1):itera])))
    }
  }
  lambda[lambda < lam_threshold] <- lam_threshold
  d1_d2 <- get_d1_d2(Y, lambda)
  #d1_lam <- sweep(lambda_k, 1, d1_d2$d1_vec, '*')
  H <- construct_hess_fast(X1,X2,lambda,lambda_k, K, d1_d2)$H
  eps = 1e-6
  if(min(eigen(H)$values) < eps) {
    I <- solve(psd(H, epsilon = eps))
  } else {
    I <- solve(H)
  }
  precision <- abs(solve(D_mat) %*% d_vec)
  precision[is.na(diag(I))] <- PRECISION.THRESHOLD + 100
  converged_vec <- check_converged_vec(X1,X2,my_beta, itera, n.iter,
                                       error_vec, precision, PRECISION.THRESHOLD)
  names(error_vec) <- colnames(my_beta)
  return(list(alpha1 = alpha1, alpha2 = alpha2, converged = any(converged_vec), I = I, d = d_vec_o,
              n.iter = itera, log_l = prev_ll, precision = precision, prediction = lambda,
              converged_vec = converged_vec, error_vec = error_vec))
}

estimate_gene_wrapper <- function(Y,X1,X2,my_beta, nUMI, sigma_init, test_mode, verbose = F, n.iter = 200, MIN_CHANGE = 1e-3, sigma_gene = T, PRECISION.THRESHOLD = 0.05) {
  if(sigma_gene)
    return(choose_sigma_gene(sigma_init, Y, X1, X2, my_beta, nUMI,test_mode, verbose = verbose, n.iter = n.iter, MIN_CHANGE = MIN_CHANGE, PRECISION.THRESHOLD = PRECISION.THRESHOLD))
  else {
    res <- estimate_effects_trust(Y,X1,X2,my_beta, nUMI,test_mode, verbose = verbose, n.iter = n.iter, MIN_CHANGE = MIN_CHANGE, PRECISION.THRESHOLD = PRECISION.THRESHOLD)
    return(list(sigma_s_best = -1, res = res))
  }
}

estimate_effects_trust <- function(Y, X1, X2, my_beta, nUMI, test_mode, verbose = F,
                                   n.iter = 200, MIN_CHANGE = 1e-3, PRECISION.THRESHOLD = 0.05,
                                   alpha1_init = NULL, alpha2_init = NULL) {
  my_beta<- sweep(my_beta,1, nUMI, '*')
  solveIRWLS.effects_trust(Y,X1,X2, my_beta, test_mode, verbose = verbose,
                           n.iter = n.iter, MIN_CHANGE = MIN_CHANGE,
                           PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                           alpha1_init = alpha1_init, alpha2_init = alpha2_init)
}
