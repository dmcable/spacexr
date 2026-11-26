test_that("UMI_cutoff decreases with nUMI and is floored at 0.25", {
  expect_equal(UMI_cutoff(32), 1)
  expect_equal(UMI_cutoff(1024), 0.25)
  expect_equal(UMI_cutoff(1e6), 0.25)
  expect_equal(UMI_cutoff(c(1, 32)), c(2, 1))
  expect_true(all(diff(UMI_cutoff(2^(0:10))) <= 0))
})

test_that("get_class_df maps each cell type to its own class by default", {
  df <- get_class_df(c("A", "B"))
  expect_equal(colnames(df), "class")
  expect_equal(rownames(df), c("A", "B"))
  expect_equal(df$class, c("A", "B"))
})

test_that("get_class_df optionally groups cerebellum cell types", {
  types <- c("Bergmann", "Astrocytes", "Fibroblast", "Endothelial", "MLI1", "MLI2",
             "Macrophages", "Microglia", "Polydendrocytes", "Oligodendrocytes")
  df <- get_class_df(types, use_classes = TRUE)
  expect_equal(df["Bergmann", "class"], "Astrocytes")
  expect_equal(df["Fibroblast", "class"], "Endothelial")
  expect_equal(df["MLI2", "class"], "MLI1")
  expect_equal(df["Macrophages", "class"], "Microglia")
  expect_equal(df["Polydendrocytes", "class"], "Oligodendrocytes")
  expect_equal(df["Astrocytes", "class"], "Astrocytes")
})

test_that("subset_rctd subsets spatial data and weights", {
  puck <- SpatialRNA(small_coords(4), small_counts(n_cols = 4))
  weights <- normalize_weights(matrix(1:8, 4, dimnames = list(paste0("bc", 1:4), c("A", "B"))))
  rctd <- make_rctd_shell(puck, results = list(weights = weights))
  expect_message(subset_rctd(rctd, c("bc1", "bc2")), "before subset: 4")
  res <- suppressMessages(subset_rctd(rctd, c("bc4", "bc2", "missing")))
  expect_equal(rownames(res@results$weights), c("bc4", "bc2"))
  expect_equal(colnames(res@spatialRNA@counts), c("bc4", "bc2"))
  expect_equal(names(res@spatialRNA@nUMI), c("bc4", "bc2"))
  expect_equal(rownames(res@spatialRNA@coords), c("bc4", "bc2"))
  expect_equal(colnames(res@originalSpatialRNA@counts), c("bc4", "bc2"))
  expect_error(suppressMessages(subset_rctd(rctd, "missing")), "No spot to keep")
})
