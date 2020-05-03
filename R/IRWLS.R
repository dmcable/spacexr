
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
  solution<-solveOLS(S,B, constrain = constrain) #first solve OLS, use this solution to find a starting point for the weights
  if(OLS)
    return(list(weights = solution, converged = T))
  solution<- (solution*0) + 1/length(solution) #actually, just ignore the OLS solution
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
  prediction = pmax(prediction, threshold)
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
