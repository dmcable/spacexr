irwls_problem <- function(w_true, n_genes = 90, umi = 2000, seed = 10) {
  P <- make_profiles()
  S <- P * umi
  B <- withr::with_seed(seed, rpois(nrow(S), S %*% w_true))
  names(B) <- rownames(S)
  list(S = S, B = B, umi = umi)
}

test_that("solveOLS recovers exact non-negative weights", {
  withr::local_seed(1)
  S <- matrix(runif(60, 0, 10), 20, 3, dimnames = list(NULL, c("a", "b", "c")))
  w <- c(0.6, 0.4, 0)
  B <- S %*% w
  res <- solveOLS(S, B, numeric(3), constrain = FALSE)
  expect_named(res, c("a", "b", "c"))
  expect_equal(unname(res), w, tolerance = 1e-4)
})

test_that("solveOLS enforces non-negativity", {
  withr::local_seed(2)
  S <- matrix(runif(60, 0, 10), 20, 3)
  B <- S %*% c(1, -0.5, 0.5)
  res <- solveOLS(S, B, numeric(3), constrain = FALSE)
  expect_true(all(res >= -1e-10))
})

test_that("solveOLS with constrain makes weights sum to one", {
  withr::local_seed(3)
  S <- matrix(runif(60, 0, 10), 20, 3)
  B <- S %*% c(2, 1, 1)
  res <- solveOLS(S, B, numeric(3), constrain = TRUE)
  expect_equal(sum(res), 1, tolerance = 1e-8)
  expect_true(all(res >= -1e-10))
})

test_that("solveIRWLS.weights in OLS mode returns converged OLS solution", {
  use_test_likelihood("100")
  p <- irwls_problem(c(0.5, 0.5, 0))
  res <- solveIRWLS.weights(p$S, p$B, p$umi, OLS = TRUE, constrain = FALSE)
  expect_true(res$converged)
  expect_equal(res$weights, solveOLS(p$S, pmin(p$B, K_val), numeric(3) + 1 / 3, constrain = FALSE))
})

test_that("solveIRWLS.weights recovers mixture proportions", {
  use_test_likelihood("100")
  w_true <- c(0.7, 0.3, 0)
  p <- irwls_problem(w_true)
  res <- solveIRWLS.weights(p$S, p$B, p$umi, constrain = FALSE)
  expect_named(res, c("weights", "converged"))
  expect_true(res$converged)
  expect_named(res$weights, colnames(p$S))
  expect_true(all(res$weights >= 0))
  expect_lt(max(abs(res$weights - w_true)), 0.15)
})

test_that("solveIRWLS.weights with constrain sums weights to one", {
  use_test_likelihood("100")
  p <- irwls_problem(c(0, 1, 0))
  res <- solveIRWLS.weights(p$S, p$B, p$umi, constrain = TRUE)
  expect_equal(sum(res$weights), 1, tolerance = 1e-6)
  expect_equal(names(which.max(res$weights)), "TypeB")
})

test_that("solveIRWLS.weights in bulk mode fits a bulk profile", {
  P <- make_profiles()
  w_true <- c(0.2, 0.5, 0.3)
  S <- P * 1e5
  B <- as.vector(S %*% w_true)
  res <- solveIRWLS.weights(S, B, 1e5, constrain = FALSE, bulk_mode = TRUE, n.iter = 100, MIN_CHANGE = 1e-6)
  expect_equal(unname(res$weights), w_true, tolerance = 1e-3)
})

test_that("solveIRWLS.weights reports non-convergence when iterations run out", {
  use_test_likelihood("100")
  p <- irwls_problem(c(0.7, 0.3, 0))
  res <- solveIRWLS.weights(p$S, p$B, p$umi, constrain = FALSE, n.iter = 1, MIN_CHANGE = 1e-12)
  expect_false(res$converged)
})

test_that("solveWLS takes a step that keeps weights non-negative", {
  use_test_likelihood("100")
  p <- irwls_problem(c(1, 0, 0))
  init <- c(0.01, 0.5, 0.49)
  names(init) <- colnames(p$S)
  invisible(solveIRWLS.weights(p$S, p$B, p$umi, n.iter = 0))  # sets S_mat
  step <- solveWLS(p$S, p$B, init, p$umi)
  expect_named(step, colnames(p$S))
  expect_true(all(step >= 0))
  expect_gt(step[["TypeA"]], init[["TypeA"]])
})
