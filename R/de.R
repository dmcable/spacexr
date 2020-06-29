estimate_effects <- function(Y, X1, X2, beta, nUMI, verbose = F, n.iter = 100, MIN_CHANGE = 0.001) {
  beta <- sweep(beta,1, nUMI, '*')
  solveIRWLS.effects(Y, X1, X2, beta, verbose = verbose, n.iter = n.iter, MIN_CHANGE = MIN_CHANGE)
}
