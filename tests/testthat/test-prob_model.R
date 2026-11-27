test_that("ht_pdf_norm is a symmetric, continuous probability density", {
  expect_equal(integrate(ht_pdf_norm, -Inf, Inf)$value, 1, tolerance = 1e-4)
  x <- c(0.1, 1, 2.5, 3, 5, 50)
  expect_equal(ht_pdf_norm(x), ht_pdf_norm(-x))
  expect_equal(ht_pdf_norm(3 - 1e-9), ht_pdf_norm(3 + 1e-9), tolerance = 1e-6)
  expect_true(all(ht_pdf_norm(x) > 0))
  expect_true(all(diff(ht_pdf_norm(x)) < 0))
  expect_length(ht_pdf_norm(numeric(0)), 0)
})

test_that("ht_pdf_norm has heavier tails than the normal density", {
  expect_gt(ht_pdf_norm(6), dnorm(6))
  expect_gt(ht_pdf_norm(0), 0.9 * dnorm(0))
})

test_that("ht_pdf rescales ht_pdf_norm by sigma", {
  z <- c(-2, 0, 0.5, 4)
  expect_equal(ht_pdf(z, 2), ht_pdf_norm(z / 2) / 2)
  expect_equal(ht_pdf(z, 1), ht_pdf_norm(z))
  expect_equal(integrate(function(z) ht_pdf(z, 0.5), -Inf, Inf)$value, 1, tolerance = 1e-4)
})

test_that("psd returns a symmetric positive definite matrix", {
  H <- matrix(c(1, 2, 2, 1), 2)
  res <- psd(H)
  expect_equal(res, t(res))
  expect_true(all(eigen(res)$values >= 1e-3 - 1e-12))

  P <- matrix(c(2, 0.5, 0.5, 1), 2)
  expect_equal(psd(P), P)
  expect_equal(psd(matrix(-1)), matrix(1e-3))
  expect_equal(psd(matrix(5)), matrix(5))
})

test_that("set_likelihood_vars sets global likelihood state", {
  use_test_likelihood("100")
  Q <- get_test_Q_mat("100")
  expect_equal(K_val, nrow(Q) - 3)
  expect_equal(N_X, ncol(Q))
  expect_equal(X_vals, get_test_X_vals())
  expect_equal(dim(SQ_mat), dim(Q))
})

test_that("solve_sq returns natural spline coefficients with zero boundaries", {
  X <- c(0, 1, 2, 3, 4)
  Q <- rbind(X^2, X)
  SQ <- solve_sq(Q, X)
  expect_equal(dim(SQ), dim(Q))
  expect_equal(unname(SQ[, 1]), c(0, 0))
  expect_equal(unname(SQ[, 5]), c(0, 0))
  expect_equal(SQ[2, ], rep(0, 5))
  ref <- splinefun(X, X^2, method = "natural")
  expect_equal(SQ[1, 2:4], ref(X[2:4], deriv = 2))
})

test_that("calc_Q_all interpolates the likelihood table at the knots", {
  use_test_likelihood("100")
  Q <- get_test_Q_mat("100")
  knots <- c(100, 200, 300, 400)
  Y <- c(0, 3, 10, 25)
  res <- calc_Q_all(Y, X_vals[knots])
  expect_named(res, c("d0_vec", "d1_vec", "d2_vec"))
  expect_equal(res$d0_vec, Q[cbind(Y + 1, knots)], tolerance = 1e-8)
})

test_that("calc_Q_all derivatives match finite differences", {
  use_test_likelihood("100")
  Y <- c(0, 1, 5, 20, 50)
  lambda <- c(0.3, 2.1, 4.7, 18.2, 47.5)
  h <- 1e-5
  res <- calc_Q_all(Y, lambda)
  up <- calc_Q_all(Y, lambda + h)
  down <- calc_Q_all(Y, lambda - h)
  expect_equal(res$d1_vec, (up$d0_vec - down$d0_vec) / (2 * h), tolerance = 1e-4)
  expect_equal(res$d2_vec, (up$d1_vec - down$d1_vec) / (2 * h), tolerance = 1e-3)
})

test_that("calc_Q_all caps counts at K_val and clamps lambda to the table range", {
  use_test_likelihood("100")
  expect_equal(calc_Q_all(K_val + 50, 30), calc_Q_all(K_val, 30))
  expect_equal(calc_Q_all(2, -5), calc_Q_all(2, 1e-4))
  expect_equal(calc_Q_all(2, max(X_vals) + 10), calc_Q_all(2, max(X_vals) - 1e-4))
})

test_that("the likelihood peaks near the observed count", {
  use_test_likelihood("100")
  grid <- seq(1, 60, by = 0.25)
  for (y in c(5, 20, 40)) {
    nll <- sapply(grid, function(l) calc_log_l_vec(l, y))
    best <- grid[which.min(nll)]
    expect_lt(abs(best - y) / y, 0.3)
  }
})

test_that("calc_log_l_vec sums the negative log likelihood", {
  use_test_likelihood("100")
  lambda <- c(1.5, 3.2, 10)
  Y <- c(1, 4, 9)
  v <- calc_log_l_vec(lambda, Y, return_vec = TRUE)
  expect_length(v, 3)
  expect_equal(v, -calc_Q_all(Y, lambda)$d0_vec)
  expect_equal(calc_log_l_vec(lambda, Y), sum(v))
  expect_true(all(v > 0))
})

test_that("calc_log_l_vec_fast approximates calc_log_l_vec", {
  use_test_likelihood("100")
  lambda <- c(0.7, 1.5, 3.2, 10, 33.3)
  Y <- c(0, 1, 4, 9, 30)
  expect_equal(calc_log_l_vec_fast(lambda, Y), calc_log_l_vec(lambda, Y), tolerance = 1e-3)
  expect_equal(calc_log_l_vec_fast(X_vals[200], 3), calc_log_l_vec(X_vals[200], 3), tolerance = 1e-8)
})

test_that("get_d1_d2 returns the derivative components of calc_Q_all", {
  use_test_likelihood("100")
  res <- get_d1_d2(c(1, 2), c(1.5, 2.5))
  full <- calc_Q_all(c(1, 2), c(1.5, 2.5))
  expect_equal(res, full[c("d1_vec", "d2_vec")])
})

build_S_mat <- function(S) {
  K <- ncol(S)
  out <- matrix(0, nrow(S), K * (K + 1) / 2)
  counter <- 1
  for (i in 1:K) for (j in i:K) { out[, counter] <- S[, i] * S[, j]; counter <- counter + 1 }
  out
}

test_that("get_der_fast gradient and Hessian match the Poisson-lognormal likelihood", {
  use_test_likelihood("100")
  withr::local_seed(4)
  S <- matrix(runif(30, 0.5, 5), 10, 3, dimnames = list(paste0("g", 1:10), c("a", "b", "c")))
  w <- c(0.5, 0.3, 0.2)
  B <- rpois(10, S %*% w)
  assign("S_mat", build_S_mat(S), envir = globalenv())
  f <- function(w) calc_log_l_vec(S %*% w, B)
  der <- get_der_fast(S, B, rownames(S), S %*% w)
  h <- 1e-5
  num_grad <- sapply(1:3, function(k) { e <- replace(numeric(3), k, h); (f(w + e) - f(w - e)) / (2 * h) })
  expect_equal(as.vector(der$grad), num_grad, tolerance = 1e-4)
  expect_equal(der$hess, t(der$hess))
  num_hess <- sapply(1:3, function(k) {
    e <- replace(numeric(3), k, h)
    (get_der_fast(S, B, rownames(S), S %*% (w + e))$grad -
       get_der_fast(S, B, rownames(S), S %*% (w - e))$grad) / (2 * h)
  })
  expect_equal(der$hess, unname(num_hess), tolerance = 1e-3)
})

test_that("get_der_fast in bulk mode differentiates the squared log error", {
  withr::local_seed(5)
  S <- matrix(runif(24, 1, 10), 8, 3)
  w <- c(0.2, 0.5, 0.3)
  B <- as.vector(S %*% c(0.3, 0.4, 0.3))
  assign("S_mat", build_S_mat(S), envir = globalenv())
  f <- function(w) sum((log(S %*% w) - log(B))^2)
  der <- get_der_fast(S, B, NULL, S %*% w, bulk_mode = TRUE)
  h <- 1e-6
  num_grad <- sapply(1:3, function(k) { e <- replace(numeric(3), k, h); (f(w + e) - f(w - e)) / (2 * h) })
  expect_equal(as.vector(der$grad), num_grad, tolerance = 1e-5)
  expect_equal(dim(der$hess), c(3L, 3L))
})
