estimate_effects <- function(Y, X1, X2, beta, nUMI, verbose = F, n.iter = 100, MIN_CHANGE = 0.001) {
  beta <- sweep(beta,1, nUMI, '*')
  solveIRWLS.effects(Y, X1, X2, beta, verbose = verbose, n.iter = n.iter, MIN_CHANGE = MIN_CHANGE)
}


solveIRWLS.effects <-function(Y, X1, X2, beta, verbose = FALSE, n.iter = 50, MIN_CHANGE = .001){
  n_cell_types <- dim(beta)[2]
  Y[Y > K_val] <- K_val
  alpha1 <- numeric(dim(X1)[2])
  alpha2 <- matrix(0,nrow = dim(X2)[2], ncol = n_cell_types) # initialize it to be the previous cell type means
  alpha2[1,] <- -5
  #alpha2[1,] <- c(-6.593454, -7.47053, -10.06267, -10.4595, -8.44947)

  iterations<-0 #now use dampened WLS, iterate weights until convergence
  change<-1;
  d <- rep(0, length(alpha1) + length(alpha2))
  while(change > MIN_CHANGE && iterations<n.iter){
    results<-solveWLS.effects(Y,X1,X2,alpha1,alpha2,beta,d)
    change<-max(abs(results$d))
    if(verbose) {
      print(paste("alpha1:",results$alpha1))
      print(results$alpha2)
    }
    alpha1 <- results$alpha1
    alpha2 <- results$alpha2
    d <- results$d
    iterations<-iterations+1
  }
  #calc H
  K <- dim(beta)[2]
  lambda_k <- exp(sweep(X2 %*% alpha2, 1, X1 %*% (alpha1),'+'))*beta #J by K
  lambda <- rowSums(lambda_k)
  threshold = 1e-8
  lambda[lambda < threshold] <- threshold
  d1_d2 <- get_d1_d2(Y, lambda)
  grad_1 <- d1_d2$d1_vec %*% sweep(X1, 1,lambda,'*')
  d1_lam <- sweep(lambda_k, 1, d1_d2$d1_vec, '*')
  grad_2 <- t(X2) %*% d1_lam
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
  return(list(alpha1 = alpha1, alpha2 = alpha2, converged = (change <= MIN_CHANGE), I = solve(H), d = d))
}


solveWLS.effects<-function(Y,X1,X2,alpha1,alpha2,beta,d) {
  K <- dim(beta)[2]
  if(min(alpha2) < -11) {
    ca = 1
  }
  lambda_k <- exp(sweep(X2 %*% alpha2, 1, X1 %*% (alpha1),'+'))*beta #J by K
  lambda <- rowSums(lambda_k)
  #print(calc_log_l_vec(lambda,Y))
  threshold = 1e-8
  lambda[lambda < threshold] <- threshold
  d1_d2 <- get_d1_d2(Y, lambda)
  grad_1 <- d1_d2$d1_vec %*% sweep(X1, 1,lambda,'*')
  d1_lam <- sweep(lambda_k, 1, d1_d2$d1_vec, '*')
  grad_2 <- t(X2) %*% d1_lam
  if(FALSE) {
    #all_grad <- cbind(sweep(X1, 1,d1_d2$d1_vec * lambda,'*'), sweep(X2,1,d1_lam,'*'))

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
    H2[(1:L1),1:L1] <- H2_11
    H2[1:L1,(L1+1):(L1+L2*K)] <- H2_12
    H2[(L1+1):(L1+L2*K),1:L1] <- t(H2_12)
    for(k in 1:K) {
      X2_Q <- sweep(X2,1, d1_lam[,k], '*')
      H2[(L1 + 1 + (k-1)*L2):(L1 + k*L2), (L1 + 1 + (k-1)*L2):(L1 + k*L2)] <- t(X2_Q) %*% X2
    }
    H <- (H1 - H2)
    #d <- solve(H) %*% c(grad_1, grad_2)
  }
  d <- 0.001*c(grad_1, grad_2)
  L1 <- dim(X1)[2]; L2 <- dim(X2)[2]
  alpha1 <- alpha1 + d[1:L1]
  alpha2 <- alpha2 + matrix(d[(L1 + 1):length(d)], nrow = L2, ncol = K)
  #print(calc_log_l_vec(lambda,Y))
  return(list(alpha1 = alpha1, alpha2 = alpha2, d = d, H = H))
}
