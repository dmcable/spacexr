
solveOLS<-function(S,B, constrain = T){
  D<-t(S)%*%S
  d<-t(S)%*%B
  norm_factor <- norm(D,"2")
  D <- D / norm_factor
  d <- d / norm_factor
  epsilon <- 1e-7; D <- D + epsilon * diag(length(d))
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  if(constrain) {
    A_const = t(rbind(1,A))
    b_const <-c(1 - sum(solution),bzero)
    solution <- quadprog::solve.QP(D,d,A_const,b_const,meq=1)$solution
  } else {
    solution <- quadprog::solve.QP(D,d,A,bzero,meq=0)$solution
  }
  names(solution)<-colnames(S)
  return(solution)
}

#solve using WLS with weights dampened by a certain dampening constant
#if constrain, constrain the weights to sum up to 1
#eta is alpha in the sparsity paper
solveIRWLS.weights <-function(S,B,nUMI, OLS=FALSE, constrain = TRUE, verbose = FALSE, n.iter = 50, MIN_CHANGE = .001){
  B[B > K_val] <- K_val
  if(OLS) {
    solution<-solveOLS(S,B, constrain = constrain) #first solve OLS, use this solution to find a starting point for the weights
    return(list(weights = solution, converged = T))
  }
  solution <- numeric(dim(S)[2])
  solution[] <- 1/length(solution) #actually, just ignore the OLS solution
  names(solution) <- colnames(S)

  S_mat <<- matrix(0,nrow = dim(S)[1],ncol = dim(S)[2]*(dim(S)[2] + 1)/2)
  counter = 1
  for(i in 1:dim(S)[2])
    for(j in i:dim(S)[2]) {
      S_mat[,counter] <<- S[,i] * S[,j] # depends on n^2
      counter <- counter + 1
    }

  iterations<-0 #now use dampened WLS, iterate weights until convergence
  changes<-c()
  change<-1;
  while(change > MIN_CHANGE && iterations<n.iter){
    new_solution<-solveWLS(S,B,solution, nUMI,TRUE, constrain=constrain)
    change<-norm(as.matrix(new_solution-solution))
    if(verbose) {
      print(paste("Change:",change))
      print(solution)
    }
    solution <- new_solution
    iterations<-iterations+1
  }
  return(list(weights = solution, converged = (change <= MIN_CHANGE)))
}

#solve WLS given a dampening constant
#for ..., think of alpha, lambda, constrain = TRUE
#either bead_mode is true and nUMI is scalar
#or bead_mode is false and nUMI is vector
solveWLS<-function(S,B,initialSol, nUMI,bead_mode,...){
  my_args = list(...)
  solution<-pmax(initialSol,0)
  prediction = abs(S%*%solution)
  threshold = max(1e-4, nUMI * 1e-7)
  prediction[prediction < threshold] <- threshold
  gene_list = rownames(S)
  derivatives <- get_der_fast(S, B, gene_list, prediction)
  d_vec <- -derivatives$grad
  D_mat <- psd(derivatives$hess)
  norm_factor <- norm(D_mat,"2")
  D_mat <- D_mat / norm_factor
  d_vec <- d_vec / norm_factor
  epsilon <- 1e-7; D_mat <- D_mat + epsilon * diag(length(d_vec))
  A<-cbind(diag(dim(S)[2]))
  bzero<- (-solution)
  alpha = 0.3
  if('constrain' %in% names(my_args) && my_args$constrain) {
    A_const = t(rbind(1,A))
    b_const <-c(1 - sum(solution),bzero)
    solution <- solution + alpha*quadprog::solve.QP(D_mat,d_vec,A_const,b_const,meq=1)$solution
  } else {
    solution <- solution + alpha*quadprog::solve.QP(D_mat,d_vec,A,bzero,meq=0)$solution
  }
  names(solution)<-colnames(S)
  return(solution)
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
