replintegrate <- function(means, sds) {
  sig_p <- estimate_tau(means, sds)
  var_t <- sds^2 + sig_p^2
  var_est <- 1/sum(1 / var_t)
  mean_est <- sum(means / var_t)*var_est
  sd_est <- sqrt(var_est)
  return(list(mean_est = mean_est, sd_est = sd_est, sig_p = sig_p))
}

replintegrate_batch <- function(mean_mat, sd_mat) {
  D <- nrow(mean_mat)
  results <- matrix(0, D, 3)
  colnames(results) <- c('mean', 'sd', 'sig_p')
  for(i in 1:D) {
    results[i,1:3] <- unlist(replintegrate(mean_mat[i,], sd_mat[i,]))
    results
  }
  return(data.frame(results))
}

### mean estimates and sds to apply shrinkage
replintegrate_shrink_estimates <- function(estimates, sds) {
  # qqnorm(estimates) # if looks pretty gaussian -> use a gaussian prior
  var_prior <- var(estimates) #max(var(estimates) - mean(sds^2), 0)
  print(paste0('Var_prior: ', var_prior))
  shrink_ratio <- var_prior / (var_prior + sds^2)
  mean_adjusted <- estimates * shrink_ratio
  var_adjusted <- shrink_ratio * (sds^2)
  Z_adj <- mean_adjusted / sqrt(var_adjusted)
  p_adj <- 2*pnorm(-abs(Z_adj)) # posterior p value
  results_df_adj <- data.frame('mean' = mean_adjusted, 'sd' = sqrt(var_adjusted), 'Z' = Z_adj, 'p' = p_adj)
  return(results_df_adj)
}

replintegrate_two_groups <- function(mean_list, sd_list, shrink_eb = T) {
  pop_est_1 <- replintegrate_batch(mean_list[[1]], sd_list[[1]])
  pop_est_2 <- replintegrate_batch(mean_list[[2]], sd_list[[2]])
  #plot(pop_est_1$mean, pop_est_2$mean)
  est_diff <- pop_est_1$mean - pop_est_2$mean
  sd_est <- sqrt(pop_est_1$sd^2 + pop_est_2$sd^2)
  if(!shrink_eb) {
    Z_est <- est_diff / sd_est
    p_val <- 2*pnorm(-abs(Z_est))
    q_val <- p.adjust(p_val, method='BH')
    results_df <- data.frame(est_diff, sd_est, Z_est, p_val, q_val,
                             pop_est_1$mean, pop_est_2$mean, pop_est_1$sd, pop_est_2$sd,
                             pop_est_1$sig_p, pop_est_2$sig_p)
    colnames(results_df) <- c('mean', 'sd', 'Z', 'p', 'q',  'm1' ,'m2', 's1', 's2', 'sig_p1', 'sig_p2')
  } else {
    results_df <- replintegrate_shrink_estimates(est_diff, sd_est)
    results_df$m1 <- pop_est_1$mean
    results_df$m2 <- pop_est_2$mean
    results_df$s1 <- pop_est_1$sd
    results_df$s2 <- pop_est_2$sd
  }
  return(results_df)
}


