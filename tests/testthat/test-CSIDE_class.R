grid_puck <- function(n = 30) {
  withr::with_seed(7, {
    counts <- matrix(rpois(5 * n, 20), 5, n, dimnames = list(paste0("g", 1:5), paste0("bc", 1:n)))
    coords <- data.frame(x = runif(n, 0, 100), y = runif(n, 0, 100), row.names = colnames(counts))
  })
  SpatialRNA(coords, counts)
}

test_that("check_designmatrix standardizes non-intercept columns", {
  X <- cbind(1, c(2, 4, 6, 10))
  rownames(X) <- paste0("bc", 1:4)
  res <- check_designmatrix(X, "f")
  expect_equal(res[, 1], c(bc1 = 1, bc2 = 1, bc3 = 1, bc4 = 1))
  expect_equal(unname(res[, 2]), c(0, 0.25, 0.5, 1))

  X2 <- cbind(c(1, 2, 3), c(0, 1, 0)); rownames(X2) <- c("a", "b", "c")
  res2 <- check_designmatrix(X2, "f")
  expect_equal(unname(res2[, 1]), c(0, 0.5, 1))
})

test_that("check_designmatrix accepts data frames and empty designs", {
  df <- data.frame(a = 1, b = c(0, 1, 2), row.names = c("x", "y", "z"))
  expect_true(is.matrix(check_designmatrix(df, "f")))
  empty <- matrix(numeric(0), nrow = 3, ncol = 0, dimnames = list(c("a", "b", "c"), NULL))
  expect_equal(dim(check_designmatrix(empty, "f")), c(3L, 0L))
  expect_error(check_designmatrix(empty, "f", require_2d = TRUE), "no covariates are present")
})

test_that("check_designmatrix rejects invalid designs", {
  expect_error(check_designmatrix(matrix(1, 1, 1, dimnames = list("a", NULL)), "f"), "only one pixel present")
  expect_error(check_designmatrix(matrix(c("a", "b"), 2, 1, dimnames = list(c("x", "y"), NULL)), "f"),
               "elements of X are not numeric")
  expect_error(check_designmatrix(matrix(1:4, 2), "f"), "rownames\\(X\\) is null")
  expect_error(check_designmatrix(matrix(1:4, 2, dimnames = list(c("a", "a"), NULL)), "f"), "duplicated")
  collinear <- cbind(1, 1:3, 2 * (1:3)); rownames(collinear) <- letters[1:3]
  expect_error(check_designmatrix(collinear, "f"), "not full rank")
  const <- cbind(c(0, 1, 2), 5); rownames(const) <- letters[1:3]
  expect_error(check_designmatrix(const, "f"), "column, 2 of X is a constant")
  expect_error(check_designmatrix(sum, "f"), "could not convert X to matrix")
})

test_that("check_cell_type_specific validates its input", {
  expect_silent(check_cell_type_specific(c(TRUE, FALSE), 2, "f"))
  expect_error(check_cell_type_specific(list(TRUE), 1, "f"), "not an atomic vector")
  expect_error(check_cell_type_specific(c(1, 0), 2, "f"), "is not numeric")
  expect_error(check_cell_type_specific(TRUE, 2, "f"), "not currently equal")
})

test_that("build.designmatrix.single creates an intercept and standardized covariate", {
  rctd <- make_rctd_shell(grid_puck())
  ev <- setNames(seq(10, 300, length.out = 30), paste0("bc", 1:30))
  ev <- c(ev, other = 5)
  X <- build.designmatrix.single(rctd, ev)
  expect_equal(dim(X), c(30L, 2L))
  expect_equal(rownames(X), paste0("bc", 1:30))
  expect_true(all(X[, 1] == 1))
  expect_equal(range(X[, 2]), c(0, 1))
  expect_equal(unname(X[, 2]), (ev[1:30] - 10) / 290, ignore_attr = TRUE)
})

test_that("build.designmatrix.single validates the explanatory variable", {
  rctd <- make_rctd_shell(grid_puck())
  expect_error(build.designmatrix.single(rctd, c(1, 2)), "names\\(explanatory.variable\\) is null")
  expect_error(build.designmatrix.single(rctd, c(bc1 = 1, x = 2)), "1 common barcode names found")
  expect_error(build.designmatrix.single(rctd, c(bc1 = 1, bc2 = 1)), "values of explanatory.variable are constant")
  expect_error(build.designmatrix.single(rctd, c(bc1 = "a", bc2 = "b")), "is not numeric")
})

test_that("build.designmatrix.intercept returns a column of ones", {
  rctd <- make_rctd_shell(grid_puck())
  X <- build.designmatrix.intercept(rctd)
  expect_equal(dim(X), c(30L, 1L))
  expect_true(all(X == 1))
  X2 <- build.designmatrix.intercept(rctd, barcodes = c("bc3", "bc1", "zzz"))
  expect_equal(rownames(X2), c("bc3", "bc1"))
  expect_error(build.designmatrix.intercept(rctd, barcodes = "bc1"), "1 common barcode names")
})

test_that("build.designmatrix.nonparam returns a standardized spline basis", {
  rctd <- make_rctd_shell(grid_puck())
  X <- build.designmatrix.nonparam(rctd, df = 6)
  expect_equal(dim(X), c(30L, 6L))
  expect_true(all(abs(X[, 1] - 1) < 1e-8))
  expect_equal(unname(colMeans(X[, 2:6])), rep(0, 5), tolerance = 1e-8)
  expect_equal(unname(apply(X[, 2:6], 2, sd)), rep(1, 5), tolerance = 1e-8)
  sub <- build.designmatrix.nonparam(rctd, barcodes = c("bc2", "bc5"), df = 6)
  expect_equal(rownames(sub), c("bc2", "bc5"))
  expect_error(build.designmatrix.nonparam(rctd, barcodes = "none"), "0 common barcode names")
})

test_that("build.designmatrix.regions builds region indicator columns", {
  rctd <- make_rctd_shell(grid_puck())
  regions <- list(paste0("bc", 1:5), paste0("bc", 6:10), paste0("bc", 11:12))
  X <- build.designmatrix.regions(rctd, regions)
  expect_equal(dim(X), c(12L, 3L))
  expect_true(all(rowSums(X) == 1))
  expect_equal(unname(colSums(X)), c(5, 5, 2))
  expect_equal(unname(X["bc7", ]), c(0, 1, 0))
})

test_that("build.designmatrix.regions validates region_list", {
  rctd <- make_rctd_shell(grid_puck())
  r <- list(paste0("bc", 1:3), paste0("bc", 4:6), paste0("bc", 7:9))
  expect_error(build.designmatrix.regions(rctd, r[[1]]), "region_list must be a list")
  expect_error(build.designmatrix.regions(rctd, r[1:2]), "Must be at least 3")
  expect_error(build.designmatrix.regions(rctd, list(1:3, r[[2]], r[[3]])), "atomic character vectors")
  expect_warning(expect_error(build.designmatrix.regions(rctd, list(c("bc1", "nope"), r[[2]], r[[3]])),
                              "length at least 2"), "not found in RCTD object")
  expect_warning(X <- build.designmatrix.regions(rctd, list(c(r[[1]], "nope"), r[[2]], r[[3]])),
                 "not found in RCTD object")
  expect_equal(nrow(X), 9)
})
