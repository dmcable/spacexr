.pipeline_cache <- new.env()

pipeline_inputs <- function() {
  if (is.null(.pipeline_cache$inputs)) {
    puck_data <- make_puck_data()
    .pipeline_cache$inputs <- list(
      reference = make_reference(),
      puck = SpatialRNA(puck_data$coords, puck_data$counts),
      truth = puck_data$truth
    )
  }
  .pipeline_cache$inputs
}

pipeline_rctd <- function() {
  if (is.null(.pipeline_cache$rctd)) {
    inp <- pipeline_inputs()
    .pipeline_cache$rctd <- suppressMessages(capture_output(
      r <- withr::with_seed(1, create.RCTD(inp$puck, inp$reference, max_cores = 1))))
    .pipeline_cache$rctd <- r
  }
  .pipeline_cache$rctd
}

run_pipeline <- function(mode) {
  key <- paste0("run_", mode)
  if (is.null(.pipeline_cache[[key]])) {
    capture_output(r <- suppressMessages(withr::with_seed(1, run.RCTD(pipeline_rctd(), doublet_mode = mode))))
    .pipeline_cache[[key]] <- r
  }
  .pipeline_cache[[key]]
}

test_that("check_vector validates named numeric vectors", {
  expect_silent(check_vector(c(a = 1, b = 2), "v", "f"))
  expect_error(check_vector(list(a = 1, b = 2), "v", "f"), "v is not an atomic vector")
  expect_error(check_vector(c(a = "x", b = "y"), "v", "f"), "v is not numeric")
  expect_error(check_vector(c(1, 2), "v", "f"), "names\\(v\\) is null")
  expect_error(check_vector(c(a = 1), "v", "f"), "the length of v is 1")
  expect_error(check_vector(c(a = 1.5, b = 2), "v", "f", require_int = TRUE), "does not contain integers")
})

test_that("process_cell_type_info requires CELL_MIN cells per type", {
  ref <- pipeline_inputs()$reference
  capture_output(info <- suppressMessages(process_cell_type_info(ref, levels(ref@cell_types))))
  expect_equal(info[[2]], test_cell_types)
  capture_output(expect_error(suppressMessages(process_cell_type_info(ref, levels(ref@cell_types), CELL_MIN = 31)),
                              "need a minimum of 31 cells"))
})

test_that("normalize_weights makes rows sum to one", {
  w <- matrix(c(1, 3, 2, 2, 0, 5), 3, dimnames = list(paste0("p", 1:3), c("A", "B")))
  res <- normalize_weights(w)
  expect_equal(unname(rowSums(res)), c(1, 1, 1))
  expect_equal(res["p1", ], c(A = 1 / 3, B = 2 / 3))
  expect_equal(dimnames(res), dimnames(w))
})

test_that("create.RCTD builds an RCTD object with DE gene lists", {
  inp <- pipeline_inputs()
  rctd <- pipeline_rctd()
  expect_s4_class(rctd, "RCTD")
  expect_equal(rctd@cell_type_info$info[[2]], test_cell_types)
  expect_null(rctd@cell_type_info$renorm)
  expect_true(all(rctd@internal_vars$gene_list_reg %in% paste0("Gene", 1:90)))
  expect_gte(length(rctd@internal_vars$gene_list_reg), 10)
  expect_gte(length(rctd@internal_vars$gene_list_bulk), length(rctd@internal_vars$gene_list_reg))
  expect_false(rctd@internal_vars$cell_types_assigned)
  expect_equal(rctd@internal_vars$class_df$class, test_cell_types)
  expect_equal(rctd@config$max_cores, 1)
  expect_setequal(rownames(rctd@spatialRNA@counts), rctd@internal_vars$gene_list_bulk)
  expect_equal(ncol(rctd@originalSpatialRNA@counts), ncol(inp$puck@counts))
  expect_true(all(table(rctd@reference@cell_types) == 5))
})

test_that("create.RCTD honors test_mode, keep_reference and cell_type_profiles", {
  inp <- pipeline_inputs()
  capture_output(rctd <- suppressMessages(create.RCTD(inp$puck, inp$reference, test_mode = TRUE, keep_reference = TRUE,
                                                     UMI_min = 100)))
  expect_equal(rctd@config$N_fit, 50)
  expect_equal(rctd@config$UMI_min, 1000)
  expect_equal(ncol(rctd@reference@counts), ncol(inp$reference@counts))
  profiles <- as.data.frame(make_profiles())
  rctd2 <- suppressMessages(create.RCTD(inp$puck, inp$reference, cell_type_profiles = profiles))
  expect_equal(rctd2@cell_type_info$info[[1]], profiles)
  expect_equal(rctd2@cell_type_info$info[[3]], 3)
})

test_that("create.RCTD errors when too few DE genes are found", {
  inp <- pipeline_inputs()
  capture_output(expect_error(suppressMessages(create.RCTD(inp$puck, inp$reference, fc_cutoff_reg = 10)),
                              "fewer than 10 regression differentially expressed genes"))
})

test_that("run.RCTD and fitPixels reject invalid modes", {
  expect_error(run.RCTD(pipeline_rctd(), doublet_mode = "bad"), "doublet_mode=bad is not a valid choice")
  rctd <- pipeline_rctd()
  rctd@internal_vars$Q_mat <- get_test_Q_mat("100")
  rctd@internal_vars$X_vals <- get_test_X_vals()
  expect_error(fitPixels(rctd, doublet_mode = "bad"), "doublet_mode=bad is not a valid choice")
})

test_that("fitBulk estimates proportions and renormalized profiles", {
  rctd <- suppressMessages(fitBulk(pipeline_rctd()))
  props <- rctd@internal_vars$proportions
  expect_named(props, test_cell_types)
  expect_equal(unname(props / sum(props)), c(140 / 3, 140 / 3, 50 / 3) / 110 * c(1, 1, 1) * 1, tolerance = 0.25)
  expect_equal(dim(rctd@cell_type_info$renorm[[1]]),
               c(length(rctd@internal_vars$gene_list_bulk), 3L))
})

test_that("run.RCTD in doublet mode classifies singlets and doublets", {
  rctd <- run_pipeline("doublet")
  truth <- pipeline_inputs()$truth
  res <- rctd@results$results_df
  expect_true(rctd@internal_vars$cell_types_assigned)
  expect_equal(rctd@config$RCTDmode, "doublet")
  expect_true(rctd@internal_vars$sigma > 0)
  expect_equal(rownames(res), colnames(rctd@spatialRNA@counts))
  expect_equal(levels(res$spot_class), c("reject", "singlet", "doublet_certain", "doublet_uncertain"))
  singlets <- names(truth)[truth != "doublet"]
  expect_gt(mean(as.character(res[singlets, "first_type"]) == truth[singlets]), 0.95)
  expect_gt(mean(res[singlets, "spot_class"] == "singlet"), 0.8)
  doublets <- names(truth)[truth == "doublet"]
  doublet_calls <- res[doublets, ]
  certain <- doublet_calls$spot_class == "doublet_certain"
  expect_gt(mean(certain), 0.6)
  pairs <- paste(pmin(as.character(doublet_calls$first_type[certain]), as.character(doublet_calls$second_type[certain])),
                 pmax(as.character(doublet_calls$first_type[certain]), as.character(doublet_calls$second_type[certain])))
  expect_gt(mean(pairs == "TypeA TypeB"), 0.9)
  expect_equal(dim(rctd@results$weights), c(nrow(res), 3L))
  expect_equal(colnames(rctd@results$weights_doublet), c("first_type", "second_type"))
  expect_equal(unname(Matrix::rowSums(rctd@results$weights_doublet)), rep(1, nrow(res)), tolerance = 1e-8)
  expect_length(rctd@results$score_mat, nrow(res))
})

test_that("get_doublet_weights and aggregate_cell_types work on doublet results", {
  rctd <- run_pipeline("doublet")
  beta <- get_doublet_weights(rctd)
  expect_equal(dim(beta), c(ncol(rctd@spatialRNA@counts), 3L))
  expect_true(all(rowSums(beta) <= 1 + 1e-8))
  counts <- aggregate_cell_types(rctd, rownames(beta))
  expect_true(all(counts[test_cell_types] > 30))
})

test_that("run.RCTD in full mode returns per-pixel weights", {
  rctd <- run_pipeline("full")
  truth <- pipeline_inputs()$truth
  w <- as.matrix(normalize_weights(rctd@results$weights))
  expect_equal(names(rctd@results), "weights")
  expect_equal(colnames(w), test_cell_types)
  singlets <- names(truth)[truth != "doublet"]
  pred <- colnames(w)[apply(w[singlets, ], 1, which.max)]
  expect_gt(mean(pred == truth[singlets]), 0.95)
  doublets <- names(truth)[truth == "doublet"]
  expect_equal(unname(colMeans(w[doublets, ])), c(0.5, 0.5, 0), tolerance = 0.15)
})

test_that("run.RCTD in multi mode returns a list of results", {
  rctd <- run_pipeline("multi")
  truth <- pipeline_inputs()$truth
  expect_length(rctd@results, ncol(rctd@spatialRNA@counts))
  first <- vapply(rctd@results, function(r) r$cell_type_list[1], character(1))
  singlets <- which(truth[colnames(rctd@spatialRNA@counts)] != "doublet")
  expect_gt(mean(first[singlets] == truth[colnames(rctd@spatialRNA@counts)][singlets]), 0.95)
  beta <- get_beta_multi(colnames(rctd@spatialRNA@counts), test_cell_types, rctd@results, rctd@spatialRNA@coords)
  expect_equal(unname(rowSums(beta)), rep(1, nrow(beta)), tolerance = 1e-8)
})

test_that("gather_results assembles per-pixel doublet results", {
  rctd <- pipeline_rctd()
  rctd@cell_type_info$renorm <- rctd@cell_type_info$info
  rctd <- restrict_puck_rctd <- rctd
  rctd@spatialRNA <- restrict_puck(rctd@spatialRNA, colnames(rctd@spatialRNA@counts)[1:2])
  make_res <- function(ft, st, cls) list(
    all_weights = c(TypeA = 0.6, TypeB = 0.4, TypeC = 0), spot_class = factor(cls, c("reject", "singlet", "doublet_certain", "doublet_uncertain")),
    first_type = ft, second_type = st, first_class = FALSE, second_class = FALSE, min_score = 10,
    singlet_score = 40, conv_all = TRUE, conv_doublet = TRUE, doublet_weights = c(0.6, 0.4),
    score_mat = matrix(1), singlet_scores = c(TypeA = 1))
  res <- gather_results(rctd, list(make_res("TypeA", "TypeB", "doublet_certain"), make_res("TypeC", "TypeA", "singlet")))
  df <- res@results$results_df
  expect_equal(rownames(df), colnames(rctd@spatialRNA@counts))
  expect_equal(as.character(df$spot_class), c("doublet_certain", "singlet"))
  expect_equal(as.character(df$first_type), c("TypeA", "TypeC"))
  expect_equal(levels(df$first_type), test_cell_types)
  expect_equal(df$singlet_score, c(40, 40))
  expect_equal(as.matrix(res@results$weights)[1, ], c(TypeA = 0.6, TypeB = 0.4, TypeC = 0))
  expect_equal(unname(as.matrix(res@results$weights_doublet)[2, ]), c(0.6, 0.4))
})

test_that("import_weights validates and stores external weights", {
  rctd <- pipeline_rctd()
  barcodes <- colnames(rctd@spatialRNA@counts)
  w <- matrix(1, length(barcodes), 3)
  expect_error(import_weights(rctd, w), "must contain rownames and colnames")
  dimnames(w) <- list(barcodes, test_cell_types)
  expect_error(import_weights(rctd, w), "rows sum to 1")
  res <- suppressMessages(import_weights(rctd, normalize_weights(w)))
  expect_equal(res@results$weights, normalize_weights(w))
  expect_equal(unname(res@internal_vars$proportions), rep(1 / 3, 3))
  expect_true(res@internal_vars$cell_types_assigned)
  expect_equal(res@internal_vars$sigma, 1)
  expect_equal(res@config$RCTDmode, "full")
  expect_equal(res@internal_vars$Q_mat, get_test_Q_mat("100"))
})

test_that("convert.old.RCTD preserves an up-to-date RCTD object", {
  rctd <- run_pipeline("full")
  res <- convert.old.RCTD(rctd)
  expect_s4_class(res, "RCTD")
  expect_equal(res@results, rctd@results)
  expect_equal(res@spatialRNA, rctd@spatialRNA)
  expect_equal(res@config, rctd@config)
})

test_that("run.CSIDE.general validates its inputs", {
  rctd <- run_pipeline("doublet")
  barcodes <- colnames(rctd@spatialRNA@counts)
  X1 <- build.designmatrix.intercept(rctd)
  X2 <- cbind(1, setNames(rctd@spatialRNA@coords$x, barcodes))
  rownames(X2) <- barcodes
  run <- function(...) suppressMessages(run.CSIDE.general(rctd, X1, X2, barcodes, cell_type_threshold = 20,
                                                          weight_threshold = 0.8, test_error = TRUE, ...))
  expect_s4_class(run(), "RCTD")
  expect_error(run(test_mode = "bad"), "not valid test_mode")
  expect_error(run(params_to_test = 5), "params_to_test must be a vector of integers")
  expect_error(run(test_mode = "categorical"), "must have values 0 or 1")
  expect_error(run(normalize_expr = TRUE, test_mode = "categorical"), "normalize_expr = TRUE is only valid")
  expect_error(run(gene_threshold = 1), "no genes past threshold")
  expect_error(suppressMessages(run.CSIDE.general(rctd, X1, X2, barcodes, cell_type_threshold = 20,
                                                  weight_threshold = 0.8, doublet_mode = FALSE,
                                                  test_error = TRUE)), NA)
  full <- run_pipeline("full")
  expect_error(suppressMessages(run.CSIDE.general(full, X1, X2, barcodes, cell_type_threshold = 20,
                                                  weight_threshold = 0.8, test_error = TRUE)),
               "RCTD was not run in doublet mode")
  unassigned <- rctd; unassigned@internal_vars$cell_types_assigned <- FALSE
  expect_error(suppressMessages(run.CSIDE.general(unassigned, X1, X2, barcodes, cell_type_threshold = 20,
                                                  weight_threshold = 0.8, test_error = TRUE)),
               "cell types have been assigned")
  expect_warning(suppressMessages(run.CSIDE.general(rctd, X1, X2, c(barcodes, "missing"), cell_type_threshold = 20,
                                                    weight_threshold = 0.8, test_error = TRUE)),
                 "Downsampling barcodes")
  expect_warning(suppressMessages(run.CSIDE.general(rctd, X1, X2, barcodes, cell_type_threshold = 5,
                                                    weight_threshold = 0.8, test_error = TRUE)),
                 "CSIDE vignette values")
  expect_error(suppressMessages(run.CSIDE.general(rctd, X1[1:10, , drop = FALSE], X2, barcodes, cell_type_threshold = 20,
                                                  weight_threshold = 0.8, test_error = TRUE)),
               "do not appear in the rownames of X1 or X2")
})

test_that("merge_RCTD_objects validates and combines replicates", {
  rctd <- pipeline_rctd()
  merged <- merge_RCTD_objects(list(rctd, rctd), c("r1", "r2"))
  expect_s4_class(merged, "RCTD.replicates")
  expect_length(merged@RCTD.reps, 2)
  expect_equal(merged@group_ids, c(r1 = 1, r2 = 1))
  expect_error(merge_RCTD_objects(list(rctd, "x"), c("r1", "r2")), "must be a list of RCTD objects")
  expect_error(merge_RCTD_objects(list(rctd), "r1"), "at least two RCTD objects")
  expect_error(merge_RCTD_objects(list(rctd, rctd), c("r1", "r2", "r3")), "must both be the same length")
  expect_error(merge_RCTD_objects(list(rctd, rctd), c("r1", "r2"), group_ids = c(1, 2)),
               "each group in group_ids must contain at least two replicates")
  expect_error(merge_RCTD_objects(list(rctd, rctd), c("r1", "r2"), group_ids = c(1.5, 1.5)),
               "does not contain integers")
})

test_that("create.RCTD.replicates builds one RCTD per replicate", {
  inp <- pipeline_inputs()
  capture_output(reps <- suppressMessages(create.RCTD.replicates(list(inp$puck, inp$puck), inp$reference,
                                                                 c("a", "b"), max_cores = 1)))
  expect_s4_class(reps, "RCTD.replicates")
  expect_length(reps@RCTD.reps, 2)
  expect_true(all(vapply(reps@RCTD.reps, is, logical(1), "RCTD")))
  expect_equal(reps@group_ids, c(a = 1, b = 1))
  capture_output({
    expect_error(suppressMessages(create.RCTD.replicates(list(inp$puck), inp$reference, "a")), "at least two SpatialRNA")
    expect_error(suppressMessages(create.RCTD.replicates(list(inp$puck, 1), inp$reference, c("a", "b"))),
                 "must be a list of SpatialRNA objects")
    expect_error(suppressMessages(create.RCTD.replicates(list(inp$puck, inp$puck), inp$reference, "a")),
                 "must both be the same length")
  })
})
