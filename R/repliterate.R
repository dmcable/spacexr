repliterate <- function(init_mean, likelihood_fn, epsilon = 1e-8, MAX.ITER = 20) {
  mean_est_prev <- -100
  posterior_mean <- init_mean
  for(i in 1:MAX.ITER) {
    ld <- likelihood_fn(posterior_mean)
    samp_dist <- likelihood_to_dist(posterior_mean, ld$d1, ld$d2)
    posterior_mean <- get_posterior_mean(samp_dist$mu, samp_dist$sigma)
    rep_results <- replintegrate(samp_dist$mu, samp_dist$sigma)
    mean_est <- rep_results$mean_est
    if(abs(mean_est - mean_est_prev) < epsilon)
      break
    mean_est_prev <- mean_est
  }
  return(rep_results)
}

get_posterior_mean <- function(means, sigmas, m = 1000) {
  sigma_p_est <- estimate_tau(means, sigmas)
  G <- length(means)
  A <- diag(c(rep(sigma_p_est^2,G))) + m^2
  B <- diag(c(rep(sigma_p_est^2,G) + sigmas^2)) + m^2
  C <- A %*% solve(B)
  posterior_mean <- C %*% means
}

likelihood_to_dist <- function(x, d1, d2) {
  sigma <- sqrt(-2 / d2)
  mu <- x - d1/d2
  return(list(mu = mu, sigma = sigma))
}
