
#return cell number, not proportion
#do not print output
solveOLS<-function(S,B){
  D<-t(S)%*%S
  d<-t(S)%*%B
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  Ap = t(rbind(1,A))
  bp <-c(1,rep(0,dim(S)[2]))
  norm_factor <- norm(D,"2")
  D <- D / norm_factor
  d <- d / norm_factor
  solution<-quadprog::solve.QP(D,d,Ap,bp,meq=1)$solution
  names(solution)<-colnames(S)
  return(solution)
}

#solve using WLS with weights dampened by a certain dampening constant
#if constrain, constrain the weights to sum up to 1
#eta is alpha in the sparsity paper
solveIRWLS.weights <-function(S,B,nUMI, OLS=FALSE, constrain = TRUE, verbose = FALSE, n.iter = 50){
  solution<-solveOLS(S,B) #first solve OLS, use this solution to find a starting point for the weights
  if(OLS)
    return(solution)
  solution<- (solution*0) + 1/length(solution) #actually, just ignore the OLS solution
  iterations<-0 #now use dampened WLS, iterate weights until convergence
  changes<-c()
  j<-19; change<-1;
  MIN_CHANGE <- .001
  while(change > MIN_CHANGE && iterations<n.iter){
    new_solution<-solveWLS(S,B,solution,j, nUMI,TRUE, constrain=constrain)
    change<-norm(as.matrix(new_solution-solution))
    if(verbose)
      print(paste("Change:",change))
    solution <- new_solution
    iterations<-iterations+1
  }
  converged <- T
  if(change > MIN_CHANGE)
    converged <- F
  return(list(solution = new_solution, converged = converged))
}

#solve WLS given a dampening constant
#for ..., think of alpha, lambda, constrain = TRUE
#either bead_mode is true and nUMI is scalar
#or bead_mode is false and nUMI is vector
solveWLS<-function(S,B,initialSol,j, nUMI,bead_mode,...){
  my_args = list(...)
  solution<-pmax(initialSol,0)
  prediction = abs(S%*%solution)
  V <- get_V(prediction)
  sigma_bar <- pmax(sqrt(V),1)
  scaled_residual = -(S%*%solution - B) / sigma_bar
  phi_d <- phi_der(scaled_residual)
  weight <- as.vector(phi_d/V)
  weight_b = as.vector(sigma_bar/V*(phi(scaled_residual) - scaled_residual*phi_d))
  W<-diag(weight)
  D_mat<-t(S)%*%W%*%S
  d_vec<- t(S)%*%W%*%B + t(S)%*%weight_b
  norm_factor <- norm(D_mat,"2")
  D_mat <- D_mat / norm_factor
  d_vec <- d_vec / norm_factor
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  if('constrain' %in% names(my_args) && my_args$constrain) {
    A_const = t(rbind(1,A))
    b_const <-c(1,rep(0,dim(S)[2]))
    solution <- quadprog::solve.QP(D_mat,d_vec,A_const,b_const,meq=1)$solution
  } else {
    solution <- quadprog::solve.QP(D_mat,d_vec,A,bzero,meq=0)$solution
  }
  names(solution)<-colnames(S)
  return(solution)
}

#derivative of phi
phi_der <- function(r, c = 1) {
  r <- abs(r)
  x <- numeric(length(r))
  x[r <= c] = 1
  x[r > c] = c*(2*r[r>c] - c)/(r[r > c]^2)
  return(x)
}

phi <- function(r, c = 1) {
  neg = r < 0
  r <- abs(r)
  x <- numeric(length(r))
  x[r <= c] = r[r <= c]
  x[r > c] = c*(2*log(r[r>c]/c) + c/r[r>c])
  x[neg] = -x[neg]
  return(x)
}

get_V <- function(x, epsilon = 1e-4) {
  V = x
  V[V < epsilon] <- epsilon
  V[V > 1] <- V[V > 1]^2
  return(V)
}

