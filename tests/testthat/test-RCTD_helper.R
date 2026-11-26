helper_profiles <- function(umi) data.matrix(make_profiles() * umi)

sim_bead <- function(w, umi = 1500, seed = 20) {
  P <- make_profiles()
  b <- withr::with_seed(seed, rpois(nrow(P), umi * (P %*% w)))
  names(b) <- rownames(P)
  b
}

test_that("decompose_full returns weights for every cell type", {
  use_test_likelihood("100")
  bead <- sim_bead(c(0, 0, 1))
  res <- decompose_full(helper_profiles(1500), 1500, bead, constrain = FALSE)
  expect_named(res, c("weights", "converged"))
  expect_named(res$weights, test_cell_types)
  expect_equal(names(which.max(res$weights)), "TypeC")
})

test_that("decompose_sparse returns normalized weights for a pair of types", {
  use_test_likelihood("100")
  bead <- sim_bead(c(0.5, 0.5, 0))
  res <- decompose_sparse(helper_profiles(1500), 1500, bead, "TypeA", "TypeB")
  expect_named(res$weights, c("TypeA", "TypeB"))
  expect_equal(sum(res$weights), 1)
  expect_equal(unname(res$weights), c(0.5, 0.5), tolerance = 0.2)
})

test_that("decompose_sparse score_mode scores the correct pair best", {
  use_test_likelihood("100")
  bead <- sim_bead(c(0.5, 0.5, 0))
  prof <- helper_profiles(1500)
  good <- decompose_sparse(prof, 1500, bead, "TypeA", "TypeB", score_mode = TRUE)
  bad <- decompose_sparse(prof, 1500, bead, "TypeA", "TypeC", score_mode = TRUE)
  expect_type(good, "double")
  expect_length(good, 1)
  expect_lt(good, bad)
})

test_that("decompose_sparse accepts a custom list of cell types", {
  use_test_likelihood("100")
  bead <- sim_bead(c(0.4, 0.3, 0.3))
  res <- decompose_sparse(helper_profiles(1500), 1500, bead, custom_list = test_cell_types)
  expect_named(res$weights, test_cell_types)
  expect_equal(sum(res$weights), 1)
})

test_that("get_prediction_sparse mixes two profiles", {
  prof <- matrix(c(1, 2, 3, 4), 2, dimnames = list(NULL, c("a", "b")))
  expect_equal(as.vector(get_prediction_sparse(prof, 1, 0.25, "a", "b")), c(0.25 * 1 + 0.75 * 3, 0.25 * 2 + 0.75 * 4))
  expect_equal(as.vector(get_prediction_sparse(prof, 1, 1, "b", "a")), c(3, 4))
})

test_that("get_singlet_score equals the likelihood of a pure profile", {
  use_test_likelihood("100")
  bead <- sim_bead(c(0, 1, 0))
  prof <- helper_profiles(1500)
  score <- get_singlet_score(prof, bead, 1500, "TypeB", constrain = TRUE)
  expect_equal(score, calc_log_l_vec(prof[, "TypeB"], bead))
  expect_lt(score, get_singlet_score(prof, bead, 1500, "TypeA", constrain = TRUE))
  vec <- get_singlet_score(prof, bead, 1500, "TypeB", constrain = TRUE, return_vec = TRUE)
  expect_length(vec, length(bead))
  expect_equal(sum(vec), score)
  expect_type(get_singlet_score(prof, bead, 1500, "TypeB", constrain = FALSE), "double")
})

test_that("decompose_doublet_fast splits counts between two cell types", {
  genes <- c("g1", "g2", "g3", "g4")
  means <- data.frame(A = c(0.4, 0.1, 0.25, 0.25), B = c(0.1, 0.4, 0.25, 0), row.names = genes)
  cti <- list(means, c("A", "B"), 2)
  bead <- c(g1 = 10, g2 = 8, g3 = 0, g4 = 6)
  res <- decompose_doublet_fast(bead, c(0.5, 0.5), genes, cti, "A", "B")
  expect_equal(res$expect_1 + res$expect_2, bead)
  expect_equal(res$expect_1[["g1"]], 10 * 0.4 / 0.5, tolerance = 1e-8)
  expect_equal(res$expect_1[["g3"]], 0)
  expect_equal(res$expect_1[["g4"]], 6, tolerance = 1e-8)
  expect_true(all(res$variance >= 0))
  expect_equal(res$variance[["g1"]], 10 * 0.8 * 0.2, tolerance = 1e-8)
})

test_that("check_pairs_type identifies whether all good pairs include a type", {
  score_mat <- matrix(c(0, 10, 50,
                        10, 0, 60,
                        50, 60, 0), 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  class_df <- data.frame(class = c("A", "B", "C"), row.names = c("A", "B", "C"))
  singlet <- c(A = 30, B = 40, C = 80)
  res <- check_pairs_type(NULL, NULL, NULL, score_mat, 10, "A", class_df, 5, TRUE, singlet)
  expect_true(res$all_pairs)
  expect_true(res$all_pairs_class)
  expect_equal(res$singlet_score, c(A = 30))
  res_c <- check_pairs_type(NULL, NULL, NULL, score_mat, 10, "C", class_df, 5, TRUE, singlet)
  expect_false(res_c$all_pairs)
  expect_false(res_c$all_pairs_class)
})

test_that("check_pairs_type uses class membership to lower the singlet score", {
  score_mat <- matrix(c(0, 10, 12,
                        10, 0, 60,
                        12, 60, 0), 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  class_df <- data.frame(class = c("X", "Y", "Y"), row.names = c("A", "B", "C"))
  singlet <- c(A = 30, B = 40, C = 20)
  res <- check_pairs_type(NULL, NULL, NULL, score_mat, 10, "B", class_df, 5, TRUE, singlet)
  expect_false(res$all_pairs)
  expect_true(res$all_pairs_class)
  expect_equal(unname(res$singlet_score), 20)
})

ct_info_from_profiles <- function() {
  P <- make_profiles()
  list(as.data.frame(P), colnames(P), ncol(P))
}

test_that("process_bead_doublet classifies a singlet pixel", {
  use_test_likelihood("100")
  cti <- ct_info_from_profiles()
  bead <- sim_bead(c(0, 1, 0), umi = 1500)
  res <- process_bead_doublet(cti, rownames(cti[[1]]), 1500, bead, class_df = get_class_df(cti[[2]]))
  expect_s3_class(res$spot_class, "factor")
  expect_equal(as.character(res$spot_class), "singlet")
  expect_equal(res$first_type, "TypeB")
  expect_named(res$all_weights, test_cell_types)
  expect_equal(sum(res$doublet_weights), 1)
  expect_true(res$conv_all)
  expect_lt(res$singlet_score - res$min_score, 20)
})

test_that("process_bead_doublet classifies a doublet pixel", {
  use_test_likelihood("100")
  cti <- ct_info_from_profiles()
  bead <- sim_bead(c(0.5, 0, 0.5), umi = 2000, seed = 21)
  res <- process_bead_doublet(cti, rownames(cti[[1]]), 2000, bead, class_df = get_class_df(cti[[2]]))
  expect_equal(as.character(res$spot_class), "doublet_certain")
  expect_setequal(c(res$first_type, res$second_type), c("TypeA", "TypeC"))
  expect_equal(unname(res$doublet_weights), c(0.5, 0.5), tolerance = 0.25)
  expect_true(all(c("TypeA", "TypeC") %in% rownames(res$score_mat)))
  expect_gt(res$singlet_score - res$min_score, 20)
})

test_that("process_bead_multi selects the cell types present", {
  use_test_likelihood("100")
  cti <- ct_info_from_profiles()
  bead <- sim_bead(c(0.4, 0.3, 0.3), umi = 3000, seed = 22)
  res <- process_bead_multi(cti, rownames(cti[[1]]), 3000, bead, class_df = get_class_df(cti[[2]]))
  expect_setequal(res$cell_type_list, test_cell_types)
  expect_named(res$conf_list, res$cell_type_list)
  expect_equal(sum(res$sub_weights), 1)
  single <- process_bead_multi(cti, rownames(cti[[1]]), 1500, sim_bead(c(1, 0, 0)), class_df = get_class_df(cti[[2]]))
  expect_equal(single$cell_type_list, "TypeA")
  limited <- process_bead_multi(cti, rownames(cti[[1]]), 3000, bead, MAX.TYPES = 2)
  expect_length(limited$cell_type_list, 2)
})
