
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
solveIRWLS.weights <-function(S,B,nUMI, OLS=FALSE, constrain = TRUE, delta = NULL, q = 0.5, s = 1, eta = 0.9, verbose = FALSE){
  solution<-solveOLS(S,B) #first solve OLS, use this solution to find a starting point for the weights
  if(OLS)
    return(solution)
  solution<- (solution*0) + 1/length(solution) #actually, just ignore the OLS solution
  iterations<-0 #now use dampened WLS, iterate weights until convergence
  changes<-c()
  j<-19; change<-1; epsilon = 1; n.iter = 50
  while(change>.001 && iterations<n.iter){
    new_solution<-solveWLS(S,B,solution,j, nUMI,TRUE, constrain=constrain, delta = delta, q=q, epsilon = epsilon)
    change<-norm(as.matrix(new_solution-solution))
    if(verbose)
      print(paste("Change:",change)) #print(paste("Granule:",solution["Granule"]))
    solution <- new_solution
    #print(names(which.max(solution)))
    #print(max(solution))
    if(! is.null(delta)) {
      epsilon <- min(epsilon, eta*(sort(solution)[length(solution) - s]))
      if(verbose)
        print(paste("Epsilon:",epsilon))
      if(epsilon < 1e-1)
        break
    }
    iterations<-iterations+1
  }
  return(new_solution)
}

#solve WLS given a dampening constant
#for ..., think of alpha, lambda, constrain = TRUE
#either bead_mode is true and nUMI is scalar
#or bead_mode is false and nUMI is vector
solveWLS<-function(S,B,initialSol,j, nUMI,bead_mode,...){
  my_args = list(...)
  solution<-pmax(initialSol,0)
  if(bead_mode) {
    prediction = S%*%solution
    threshold<-1e-4
    prediction[which(prediction < threshold)] = threshold
  } else {
    throw("Not currently supporting non bead_mode")
    #v[which(v < threshold)] = threshold[which(v < threshold)]
  }
  #robust regression
  prev_prediction = prediction
  scaled_residual = abs(S%*%solution - B)
  scaled_residual[prediction > 1] = scaled_residual[prediction > 1] / prediction[prediction > 1]
  exceed_sr = scaled_residual[scaled_residual > 1]
  prediction[scaled_residual > 1] = prediction[scaled_residual > 1] * (exceed_sr^2/(2*exceed_sr - 1))
  prediction[prev_prediction > 1] = prediction[prev_prediction > 1] * prev_prediction[prev_prediction > 1]
  #end robust region of code
  W<-diag(as.vector(1/prediction))
  D_mat<-t(S)%*%W%*%S
  d_vec<- t(S)%*%W%*%B
  norm_factor = sum(abs(B))
  if(! bead_mode) {
    norm_factor = norm_factor / sum(abs(my_args$alpha))
  }
  D_mat <- D_mat / norm_factor
  d_vec <- d_vec / norm_factor
  if('lambda' %in% names(my_args) && my_args$lambda > 0) {
    if(! 'alpha' %in% names(my_args) && max(my_args$alpha) > 0)
      throw("Used regularization option lambda without specifying alpha")
    sum_solution <- sum(solution)
    D_mat <- D_mat + diag(dim(S)[2]) * my_args$lambda  / (sum_solution ^ 2)
    d_vec <- d_vec + (my_args$lambda / sum_solution) * my_args$alpha
  }
  if(FALSE) { #'delta' %in% names(my_args) && my_args$delta > 0
    if(! ('q' %in% names(my_args) || 'epsilon' %in% names(my_args)))
      throw("Used sparsity option delta without specifying both epsilon and q")
    D_mat <- D_mat + diag(my_args$q*my_args$delta/((sol**2 + my_args$epsilon**2)**(1-my_args$q/2)))
  }
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
