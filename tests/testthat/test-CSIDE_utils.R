doublet_results_fixture <- function() {
  barcodes <- paste0("bc", 1:5)
  types <- c("A", "B", "C")
  results_df <- data.frame(
    spot_class = factor(c("singlet", "doublet_certain", "doublet_uncertain", "reject", "doublet_certain"),
                        levels = c("reject", "singlet", "doublet_certain", "doublet_uncertain")),
    first_type = factor(c("A", "B", "C", "A", "A"), levels = types),
    second_type = factor(c("B", "C", "A", "B", "C"), levels = types),
    row.names = barcodes)
  weights_doublet <- matrix(c(0.9, 0.6, 0.7, 0.5, 0.3,
                              0.1, 0.4, 0.3, 0.5, 0.7), ncol = 2,
                            dimnames = list(barcodes, c("first_type", "second_type")))
  list(barcodes = barcodes, types = types, results_df = results_df, weights_doublet = weights_doublet)
}

test_that("get_beta_doublet converts doublet results into a weight matrix", {
  f <- doublet_results_fixture()
  beta <- get_beta_doublet(f$barcodes, f$types, f$results_df, f$weights_doublet)
  expect_equal(dim(beta), c(5L, 3L))
  expect_equal(beta["bc1", ], c(A = 1, B = 0, C = 0))
  expect_equal(beta["bc2", ], c(A = 0, B = 0.6, C = 0.4))
  expect_equal(beta["bc3", ], c(A = 0, B = 0, C = 0))
  expect_equal(beta["bc4", ], c(A = 0, B = 0, C = 0))
  expect_equal(beta["bc5", ], c(A = 0.3, B = 0, C = 0.7))
})

test_that("get_doublet_weights uses the RCTD results", {
  f <- doublet_results_fixture()
  puck <- SpatialRNA(small_coords(5), small_counts(n_cols = 5))
  rctd <- make_rctd_shell(puck, results = list(results_df = f$results_df, weights_doublet = f$weights_doublet),
                          cell_type_info = list(info = list(NULL, f$types, 3)))
  expect_equal(get_doublet_weights(rctd), get_beta_doublet(f$barcodes, f$types, f$results_df, f$weights_doublet))
})

multi_results_fixture <- function() {
  list(list(cell_type_list = c("A", "C"), sub_weights = c(0.6, 0.4)),
       list(cell_type_list = "B", sub_weights = 1),
       list(cell_type_list = c("C", "B", "A"), sub_weights = c(0.5, 0.3, 0.2)))
}

test_that("get_beta_multi converts multi-mode results into a weight matrix", {
  coords <- small_coords(3)
  beta <- get_beta_multi(c("bc3", "bc1"), c("A", "B", "C"), multi_results_fixture(), coords)
  expect_equal(rownames(beta), c("bc3", "bc1"))
  expect_equal(beta["bc1", ], c(A = 0.6, B = 0, C = 0.4))
  expect_equal(beta["bc3", ], c(A = 0.2, B = 0.3, C = 0.5))
  expect_error(get_beta_multi("bc1", c("A", "B", "C"), multi_results_fixture(), small_coords(2)),
               "must be the same length")
})

test_that("set_cell_types_assigned flags the RCTD object", {
  rctd <- make_rctd_shell(SpatialRNA(small_coords(), small_counts()))
  expect_null(rctd@internal_vars$cell_types_assigned)
  expect_true(set_cell_types_assigned(rctd)@internal_vars$cell_types_assigned)
})

test_that("filter_barcodes_cell_types keeps pixels dominated by selected types", {
  beta <- matrix(c(0.5, 0.5, 0,
                   0.2, 0.3, 0.5,
                   1, 0, 0), ncol = 3, byrow = TRUE,
                 dimnames = list(c("p1", "p2", "p3"), c("A", "B", "C")))
  res <- filter_barcodes_cell_types(c("p1", "p2", "p3"), c("A", "B"), beta)
  expect_equal(res$barcodes, c("p1", "p3"))
  expect_equal(dim(res$my_beta), c(2L, 2L))
  res2 <- filter_barcodes_cell_types(c("p1", "p2", "p3"), "A", beta, thresh = 0.4)
  expect_equal(res2$barcodes, c("p1", "p3"))
  expect_equal(colnames(res2$my_beta), "A")
})

test_that("filter_genes keeps genes above an expression threshold and drops mt genes", {
  counts <- matrix(c(90, 80, 70,
                     1, 0, 0,
                     9, 20, 30), nrow = 3, byrow = TRUE,
                   dimnames = list(c("g1", "mt-g2", "g3"), paste0("bc", 1:3)))
  counts2 <- rbind(counts, g4 = c(0, 0, 0))
  puck <- SpatialRNA(small_coords(), counts2)
  genes <- suppressMessages(filter_genes(puck, threshold = 1e-3))
  expect_equal(genes, c("g1", "g3"))
  expect_equal(suppressMessages(filter_genes(puck, threshold = 0.5)), "g1")
  expect_equal(suppressMessages(filter_genes(puck, threshold = 1e-3, batch_size = 2)), c("g1", "g3"))
})

test_that("get_param_names and get_cell_type_ind label CSIDE parameters", {
  X1 <- matrix(1, 4, 1)
  X2 <- matrix(1, 4, 2)
  expect_equal(get_param_names(X1, X2, c("A", "B")), c("1_1", "2_1_A", "2_2_A", "2_1_B", "2_2_B"))
  expect_equal(get_cell_type_ind(X1, X2, 2), c(0, 1, 1, 2, 2))
  X1_empty <- matrix(numeric(0), 4, 0)
  expect_equal(get_param_names(X1_empty, X2, "A"), c("2_1_A", "2_2_A"))
  expect_equal(get_cell_type_ind(X1_empty, X2, 1), c(1, 1))
})

test_that("check_converged_vec flags imprecise or errored cell types", {
  X1 <- matrix(numeric(0), 4, 0)
  X2 <- matrix(1, 4, 2)
  beta <- matrix(0, 4, 3, dimnames = list(NULL, c("A", "B", "C")))
  precision <- c(0.01, 0.01, 0.5, 0.01, 0.01, 0.01)
  res <- check_converged_vec(X1, X2, beta, 1, 10, c(FALSE, FALSE, TRUE), precision, 0.05)
  expect_equal(res, c(A = TRUE, B = FALSE, C = FALSE))
})

test_that("fdr_sig_genes applies Benjamini-Hochberg", {
  genes <- paste0("g", 1:6)
  p <- c(0.001, 0.2, 0.004, 0.03, 0.9, 0.012)
  res <- fdr_sig_genes(genes, p, 0.05)
  expect_setequal(res, genes[p.adjust(p, "BH") < 0.05])
  expect_equal(res, c("g1", "g3", "g6", "g4"))
  expect_null(fdr_sig_genes(genes, rep(0.5, 6), 0.05))
})

test_that("fdr_sig_genes supports local fdr", {
  skip_if_not_installed("locfdr")
  # locfdr is listed in Imports but not imported in NAMESPACE, so it must be attached.
  withr::local_package("locfdr")
  withr::local_seed(11)
  genes <- paste0("g", 1:500)
  Z <- c(rnorm(480), rnorm(20, 6))
  res <- suppressWarnings(fdr_sig_genes(genes, NULL, 0.1, Z = Z, method = "locfdr"))
  expect_true(all(res %in% genes))
  expect_gt(mean(res %in% genes[481:500]), 0.8)
})

test_that("get_l_chi and get_p_chi compute chi-squared tail probabilities", {
  S_inv <- diag(2)
  points <- cbind(c(1, 0), c(0, 1), c(1, 1))
  expect_equal(get_l_chi(3, 2, S_inv, points), 0)
  expect_equal(as.numeric(get_l_chi(3, 2, S_inv, points, delta = 0.5)), 2 * 0.25)
  expect_equal(get_p_chi(3, 2, S_inv, points), 1 - pchisq(3, 2))
  expect_equal(get_p_chi(3, 2, S_inv, points, delta = 0.5), 1 - pchisq(3, 2, ncp = 0.5))
})

test_that("get_con_regions extracts per-region convergence for a cell type", {
  gene_fits <- list(con_all = matrix(c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE), 1, dimnames = list("g1", NULL)))
  expect_equal(get_con_regions(gene_fits, "g1", 3, 1, 2), c(TRUE, FALSE, TRUE))
  expect_equal(get_con_regions(gene_fits, "g1", 3, 2, 2), c(TRUE, TRUE, FALSE))
})

test_that("normalize_ev rescales to [0, 1] using the 85th percentile", {
  ev <- c(a = 10, b = 20, c = 30, d = 40, e = 100)
  res <- normalize_ev(ev)
  expect_equal(min(res), 0)
  expect_equal(max(res), 1)
  shifted <- ev - 10
  expect_equal(res, pmin(shifted / quantile(shifted, 0.85), 1))
  expect_named(res, names(ev))
})

test_that("exvar.point.density is higher near the points", {
  coords <- data.frame(x = c(0, 10, 50, 100), y = 0, row.names = paste0("bc", 1:4))
  puck <- SpatialRNA(coords, small_counts(n_cols = 4))
  rctd <- make_rctd_shell(puck)
  ev <- exvar.point.density(rctd, rownames(coords), data.frame(x = 0, y = 0), radius = 20)
  expect_named(ev, rownames(coords))
  expect_true(all(diff(ev) <= 0))
  expect_equal(unname(ev[1]), 1)
  expect_equal(unname(ev[4]), 0)
})

test_that("exvar.celltocell.interactions scores proximity to a cell type", {
  coords <- data.frame(x = c(0, 5, 50, 100, 200), y = 0, row.names = paste0("bc", 1:5))
  puck <- SpatialRNA(coords, small_counts(n_cols = 5))
  f <- doublet_results_fixture()
  rctd <- make_rctd_shell(puck, results = list(results_df = f$results_df, weights_doublet = f$weights_doublet))
  ev <- exvar.celltocell.interactions(rctd, f$barcodes, "B", radius = 10)
  expect_named(ev, f$barcodes)
  expect_true(all(ev >= 0 & ev <= 1))
  expect_equal(unname(ev[5]), 0)
  expect_gt(ev[["bc1"]], ev[["bc3"]])
})

full_mode_rctd <- function() {
  barcodes <- paste0("bc", 1:6)
  weights <- matrix(c(0.9, 0.1, 0,
                      0.85, 0.15, 0,
                      0.1, 0.9, 0,
                      0.2, 0.2, 0.6,
                      0.05, 0.9, 0.05,
                      1, 0, 0), ncol = 3, byrow = TRUE,
                    dimnames = list(barcodes, c("A", "B", "C")))
  puck <- SpatialRNA(small_coords(6), small_counts(n_cols = 6))
  make_rctd_shell(puck, results = list(weights = weights), config = list(RCTDmode = "full"),
                  cell_type_info = list(info = list(NULL, c("A", "B", "C"), 3)))
}

test_that("aggregate_cell_types counts doublet-mode assignments", {
  f <- doublet_results_fixture()
  puck <- SpatialRNA(small_coords(5), small_counts(n_cols = 5))
  rctd <- make_rctd_shell(puck, results = list(results_df = f$results_df), config = list(RCTDmode = "doublet"))
  counts <- aggregate_cell_types(rctd, f$barcodes)
  expect_equal(as.vector(counts[c("A", "B", "C")]), c(2, 1, 2))
  expect_error(aggregate_cell_types(full_mode_rctd(), paste0("bc", 1:6)), "RCTD was not run in doublet mode")
})

test_that("aggregate_cell_types sums weights in full and multi mode", {
  rctd <- full_mode_rctd()
  expect_equal(aggregate_cell_types(rctd, paste0("bc", 1:3), doublet_mode = FALSE),
               c(A = 1.85, B = 1.15, C = 0))
  puck <- SpatialRNA(small_coords(3), small_counts(n_cols = 3))
  multi <- make_rctd_shell(puck, results = multi_results_fixture(), config = list(RCTDmode = "multi"),
                           cell_type_info = list(info = list(NULL, c("A", "B", "C"), 3)))
  expect_equal(aggregate_cell_types(multi, c("bc1", "bc2"), doublet_mode = FALSE), c(A = 0.6, B = 1, C = 0.4))
})

test_that("choose_cell_types selects cell types above the threshold", {
  rctd <- full_mode_rctd()
  barcodes <- paste0("bc", 1:6)
  beta <- rctd@results$weights
  expect_message(res <- choose_cell_types(rctd, barcodes, FALSE, 2, NULL, beta, 0.8, NULL), "detected the following")
  expect_equal(res, c("A", "B"))
  expect_equal(choose_cell_types(rctd, barcodes, FALSE, 2, "A", beta, 0.8, NULL), "A")
  expect_error(choose_cell_types(rctd, barcodes, FALSE, 2, "C", beta, 0.8, NULL), "less than the minimum")
  expect_error(suppressMessages(choose_cell_types(rctd, barcodes, FALSE, 100, NULL, beta, 0.8, NULL)),
               "no cell types occured greater")
  filt <- c(A = TRUE, B = FALSE)
  expect_message(res2 <- choose_cell_types(rctd, barcodes, FALSE, 2, "A", beta, 0.8, filt), NA)
  expect_equal(suppressMessages(choose_cell_types(rctd, barcodes, FALSE, 2, NULL, beta, 0.8, filt)), "A")
})

test_that("count_cell_types sums weights over pixels passing the threshold", {
  rctd <- full_mode_rctd()
  barcodes <- paste0("bc", 1:6)
  res <- suppressMessages(count_cell_types(rctd, barcodes, c("A", "B"), cell_type_threshold = 2, doublet_mode = FALSE))
  expect_equal(res, c(A = 2.9, B = 2.05), tolerance = 1e-8)
  res2 <- count_cell_types(rctd, barcodes, c("A", "B"), cell_type_threshold = 2, doublet_mode = FALSE, weight_threshold = 0.99)
  expect_equal(res2, c(A = 2.85, B = 1.15), tolerance = 1e-8)
})

test_that("get_spline_matrix builds a df-column standardized basis", {
  withr::local_seed(12)
  n <- 40
  counts <- matrix(rpois(3 * n, 10), 3, n, dimnames = list(paste0("g", 1:3), paste0("bc", 1:n)))
  puck <- SpatialRNA(data.frame(x = runif(n), y = runif(n), row.names = colnames(counts)), counts)
  X <- get_spline_matrix(puck, df = 5)
  expect_equal(dim(X), c(n, 5L))
  expect_equal(rownames(X), colnames(counts))
  expect_true(all(abs(X[, 1] - 1) < 1e-8))
  expect_equal(Matrix::rankMatrix(X)[1], 5)
})

test_that("get_standard_errors reshapes the standard error matrix", {
  s_mat <- matrix(1:8, 2, 4, dimnames = list(c("g1", "g2"), NULL))
  rctd <- make_rctd_shell(SpatialRNA(small_coords(), small_counts()))
  rctd@de_results <- list(gene_fits = list(s_mat = s_mat))
  rctd@internal_vars_de <- list(cell_types = c("A", "B"))
  res <- get_standard_errors(rctd)
  expect_equal(dim(res), c(2L, 2L, 2L))
  expect_equal(dimnames(res)[[3]], c("A", "B"))
  expect_equal(res["g2", 2, "B"], 8)
})
