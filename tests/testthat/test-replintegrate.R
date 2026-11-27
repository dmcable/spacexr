test_that("estimate_tau is zero when variation is explained by the standard errors", {
  expect_equal(estimate_tau(c(1, 1, 1), c(1, 1, 1)), 0)
  expect_equal(estimate_tau(c(0.1, -0.1, 0.05, 0), rep(1, 4)), 0)
})

test_that("estimate_tau recovers between-replicate variance", {
  withr::local_seed(30)
  x <- rnorm(200, 0, 2) + rnorm(200, 0, 0.5)
  tau <- estimate_tau(x, rep(0.5, 200))
  tau_small <- estimate_tau(rnorm(200, 0, 0.5) + rnorm(200, 0, 0.5), rep(0.5, 200))
  expect_gt(tau, 1)
  expect_lt(tau, 2.5)
  expect_lt(tau_small, tau)
})

test_that("estimate_tau warns when it does not converge", {
  withr::local_seed(31)
  expect_warning(estimate_tau(rnorm(10, 0, 5), rep(0.1, 10), n.iter = 1), "not converged")
})

test_that("estimate_tau_group averages per-group estimates", {
  x <- c(1, 1, 1, 5, -5, 0)
  s <- rep(1, 6)
  expect_equal(estimate_tau_group(x, s), estimate_tau(x, s))
  g <- c(1, 1, 1, 2, 2, 2)
  expect_equal(estimate_tau_group(x, s, group_ids = g),
               mean(c(estimate_tau(x[1:3], s[1:3]), estimate_tau(x[4:6], s[4:6]))))
})

test_that("get_p_qf returns a valid p-value that is small for heterogeneous estimates", {
  p_null <- get_p_qf(c(0.1, -0.1, 0.05, 0), rep(1, 4))
  p_alt <- get_p_qf(c(5, -5, 4, -4), rep(0.5, 4))
  expect_true(p_null >= 0 && p_null <= 1)
  expect_gt(p_null, 0.5)
  expect_lt(p_alt, 1e-3)
  expect_gt(get_p_qf(c(5, -5, 4, -4), rep(0.5, 4), delta = 5), p_alt)
})

test_that("replintegrate computes an inverse-variance weighted mean", {
  res <- replintegrate(c(1, 1, 1), c(1, 1, 1))
  expect_named(res, c("mean_est", "sd_est", "sig_p"))
  expect_equal(res$mean_est, 1)
  expect_equal(res$sd_est, 1 / sqrt(3))
  expect_equal(res$sig_p, 0)
  res2 <- replintegrate(c(0, 1), c(1, 0.5))
  expect_equal(res2$mean_est, (0 / 1 + 1 / 0.25) / (1 + 4))
})

test_that("replintegrate inflates the standard error under heterogeneity", {
  het <- replintegrate(c(-3, 0, 3), c(0.5, 0.5, 0.5))
  expect_gt(het$sig_p, 0)
  expect_gt(het$sd_est, 0.5 / sqrt(3))
  expect_equal(het$mean_est, 0)
})

test_that("replintegrate_batch applies replintegrate row-wise", {
  means <- rbind(c(1, 1, 1), c(-3, 0, 3))
  sds <- rbind(c(1, 1, 1), c(0.5, 0.5, 0.5))
  res <- replintegrate_batch(means, sds)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("mean", "sd", "sig_p"))
  expect_equal(unlist(res[2, ]), unlist(replintegrate(means[2, ], sds[2, ])), ignore_attr = TRUE)
})

test_that("replintegrate_shrink_estimates shrinks toward zero", {
  est <- c(2, -1, 0.5, 3)
  sds <- c(1, 1, 2, 0.5)
  res <- expect_output(replintegrate_shrink_estimates(est, sds), "Var_prior")
  ratio <- var(est) / (var(est) + sds^2)
  expect_equal(res$mean, est * ratio)
  expect_equal(res$sd, sqrt(ratio * sds^2))
  expect_equal(res$p, 2 * pnorm(-abs(res$Z)))
  expect_true(all(abs(res$mean) <= abs(est)))
})

test_that("replintegrate_two_groups compares two groups of replicates", {
  m1 <- rbind(c(1, 1.1, 0.9), c(0, 0.1, -0.1))
  m2 <- rbind(c(0, 0.1, -0.1), c(0, 0.1, -0.1))
  s <- matrix(0.1, 2, 3)
  res <- replintegrate_two_groups(list(m1, m2), list(s, s), shrink_eb = FALSE)
  expect_named(res, c("mean", "sd", "Z", "p", "q", "m1", "m2", "s1", "s2", "sig_p1", "sig_p2"))
  expect_equal(res$mean, c(1, 0), tolerance = 1e-8)
  expect_lt(res$p[1], 1e-6)
  expect_gt(res$p[2], 0.5)
  expect_equal(res$q, p.adjust(res$p, "BH"))
  res_eb <- suppressWarnings(capture.output(r <- replintegrate_two_groups(list(m1, m2), list(s, s))))
  expect_true(all(c("mean", "sd", "Z", "p", "m1", "m2", "s1", "s2") %in% names(r)))
})

test_that("likelihood_to_dist converts derivatives to a normal approximation", {
  # derivatives of a quadratic with optimum mu and curvature -2 / sigma^2, evaluated at x
  x <- 1; mu <- 3; sigma <- 2
  d2 <- -2 / sigma^2
  d1 <- (x - mu) * d2
  res <- likelihood_to_dist(x, d1, d2)
  expect_equal(res$mu, mu)
  expect_equal(res$sigma, sigma)
  res2 <- likelihood_to_dist(c(0, 0), c(2, 4), c(-2, -8))
  expect_equal(res2$mu, c(1, 0.5))
  expect_equal(res2$sigma, c(1, 0.5))
})

test_that("get_posterior_mean shrinks replicate means toward each other", {
  means <- c(-1, 0, 1)
  sigmas <- c(1, 1, 1)
  pm <- get_posterior_mean(means, sigmas)
  expect_equal(dim(pm), c(3L, 1L))
  expect_equal(mean(pm), 0, tolerance = 1e-6)
  expect_true(all(abs(pm) <= abs(means) + 1e-8))
})

test_that("repliterate converges for a quadratic log-likelihood", {
  targets <- c(1, 1.5, 2)
  lik <- function(x) list(d1 = -2 * (x - targets), d2 = rep(-2, length(x)))
  res <- repliterate(c(0, 0, 0), lik)
  expect_named(res, c("mean_est", "sd_est", "sig_p"))
  expect_equal(res$mean_est, 1.5, tolerance = 1e-6)
})
