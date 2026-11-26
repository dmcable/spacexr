test_that("SpatialRNA builds a valid object from a dense matrix", {
  counts <- small_counts()
  puck <- SpatialRNA(small_coords(), counts)
  expect_s4_class(puck, "SpatialRNA")
  expect_s4_class(puck@counts, "dgCMatrix")
  expect_equal(as.matrix(puck@counts), counts)
  expect_equal(puck@nUMI, colSums(counts))
  expect_named(puck@nUMI, colnames(counts))
  expect_named(puck@coords, c("x", "y"))
  expect_equal(rownames(puck@coords), colnames(counts))
})

test_that("SpatialRNA accepts data.frame counts and renames coordinate columns", {
  coords <- small_coords()
  colnames(coords) <- c("xcoord", "ycoord")
  puck <- SpatialRNA(coords, as.data.frame(small_counts()))
  expect_s4_class(puck@counts, "dgCMatrix")
  expect_named(puck@coords, c("x", "y"))
})

test_that("SpatialRNA keeps a user-provided nUMI and warns on mismatch", {
  counts <- small_counts()
  nUMI <- colSums(counts) * 2
  expect_warning(puck <- SpatialRNA(small_coords(), counts, nUMI = nUMI),
                 "nUMI does not match colSums")
  expect_equal(puck@nUMI, nUMI)
})

test_that("SpatialRNA keeps only mutually shared barcodes", {
  counts <- small_counts(n_cols = 4)
  coords <- small_coords(3)
  expect_warning(puck <- SpatialRNA(coords, counts), "not mutually shared")
  expect_equal(colnames(puck@counts), paste0("bc", 1:3))
  expect_equal(names(puck@nUMI), paste0("bc", 1:3))
  expect_equal(rownames(puck@coords), paste0("bc", 1:3))
})

test_that("SpatialRNA errors when no barcodes are shared", {
  expect_error(SpatialRNA(small_coords(prefix = "other"), small_counts()),
               "do not share any barcode names")
})

test_that("SpatialRNA can generate fake coordinates", {
  puck <- SpatialRNA(NULL, small_counts(), use_fake_coords = TRUE)
  expect_equal(dim(puck@coords), c(3L, 2L))
  expect_true(all(puck@coords == 0))
  expect_equal(rownames(puck@coords), colnames(puck@counts))
})

test_that("fake_coords returns zero x/y coordinates named by barcode", {
  coords <- fake_coords(small_counts(n_cols = 5))
  expect_s3_class(coords, "data.frame")
  expect_named(coords, c("x", "y"))
  expect_equal(rownames(coords), paste0("bc", 1:5))
  expect_true(all(coords == 0))
})

test_that("SpatialRNA validates count matrices", {
  coords <- small_coords()
  one_gene <- small_counts(n_genes = 1)
  expect_error(SpatialRNA(coords, one_gene), "only one gene present")

  frac <- small_counts() + 0.5
  expect_error(SpatialRNA(coords, frac), "counts does not contain integers")
  expect_s4_class(SpatialRNA(coords, frac, require_int = FALSE), "SpatialRNA")

  no_genes <- small_counts(); rownames(no_genes) <- NULL
  expect_error(SpatialRNA(coords, no_genes), "rownames\\(counts\\) is null")

  no_bc <- small_counts(); colnames(no_bc) <- NULL
  expect_error(SpatialRNA(coords, no_bc), "colnames\\(counts\\) is null")

  dup_genes <- small_counts(); rownames(dup_genes)[2] <- "g1"
  expect_error(SpatialRNA(coords, dup_genes), "rownames\\(counts\\) contain duplicated")

  dup_bc <- small_counts(); colnames(dup_bc)[2] <- "bc1"
  expect_error(SpatialRNA(coords, dup_bc), "colnames\\(counts\\) contain duplicated")

  na_genes <- small_counts(); rownames(na_genes)[2] <- NA
  expect_error(SpatialRNA(coords, na_genes), "rownames\\(counts\\) contain NA")
})

test_that("SpatialRNA supports a single pixel", {
  expect_warning(puck <- SpatialRNA(small_coords(1), small_counts(n_cols = 1)), "only one cell/pixel present")
  expect_s4_class(puck@counts, "dgCMatrix")
  expect_equal(dim(puck@counts), c(4L, 1L))
  expect_equal(rownames(puck@coords), "bc1")
  expect_equal(puck@nUMI, c(bc1 = 10))
})

test_that("SpatialRNA accepts matrix coordinates", {
  puck <- SpatialRNA(as.matrix(small_coords()), small_counts())
  expect_s3_class(puck@coords, "data.frame")
  expect_named(puck@coords, c("x", "y"))
  expect_equal(puck@coords$y, c(2, 4, 6))
})

test_that("check_counts warns for a single pixel", {
  expect_warning(check_counts(small_counts(n_cols = 1), "SpatialRNA"), "only one cell/pixel present")
})

test_that("SpatialRNA validates coordinates", {
  counts <- small_counts()
  three_col <- cbind(small_coords(), z = 1:3)
  expect_error(SpatialRNA(three_col, counts), "second dimension of coords is not 2")

  char_coords <- data.frame(x = letters[1:3], y = 1:3, row.names = paste0("bc", 1:3))
  expect_error(SpatialRNA(char_coords, counts), "coords is not numeric")

  expect_error(SpatialRNA(sum, counts), "could not convert coords")
})

test_that("check_UMI validates nUMI vectors", {
  expect_error(check_UMI(list(a = 1, b = 2), "f"), "not an atomic vector")
  expect_error(check_UMI(c(a = "1", b = "2"), "f"), "nUMI is not numeric")
  expect_error(check_UMI(c(a = 1.5, b = 2), "f"), "does not contain integers")
  expect_silent(check_UMI(c(a = 1.5, b = 2), "f", require_int = FALSE))
  expect_error(check_UMI(c(1, 2), "f"), "names\\(nUMI\\) is null")
  expect_error(check_UMI(c(a = 1, a = 2), "f"), "duplicated elements")
  expect_warning(check_UMI(c(a = 1), "f"), "length of nUMI is 1")
  expect_error(check_UMI(c(a = 1), "f", require_2d = TRUE), "length of nUMI is 1")
  expect_error(check_UMI(c(a = 1, b = 2), "f", min_UMI = 10), "all less than min_UMI")
  expect_warning(check_UMI(c(a = 1, b = 20), "f", min_UMI = 10), "will be removed")
  expect_silent(check_UMI(c(a = 10, b = 20), "f", min_UMI = 10))
})

test_that("check_counts errors on a single cell when 2d is required", {
  expect_error(check_counts(small_counts(n_cols = 1), "f", require_2d = TRUE),
               "only one cell present")
})

test_that("restrict_counts filters genes and pixels by UMI thresholds", {
  counts <- matrix(c(10, 0, 0,
                     0, 5, 0,
                     1, 1, 100,
                     0, 0, 3), nrow = 4, byrow = TRUE,
                   dimnames = list(paste0("g", 1:4), paste0("bc", 1:3)))
  puck <- SpatialRNA(small_coords(), counts)
  # nUMI = c(11, 6, 103)
  res <- restrict_counts(puck, c("g1", "g2"), UMI_thresh = 5, UMI_max = 50)
  expect_equal(rownames(res@counts), c("g1", "g2"))
  expect_equal(colnames(res@counts), c("bc1", "bc2"))
  expect_equal(unname(res@nUMI), c(11, 6))

  res2 <- restrict_counts(puck, c("g1", "g3"), counts_thresh = 6, UMI_max = 1000)
  expect_equal(colnames(res2@counts), c("bc1", "bc3"))
})

test_that("restrict_counts keeps matrix shape for one gene or one pixel", {
  puck <- SpatialRNA(small_coords(), small_counts())
  one_gene <- restrict_counts(puck, "g2")
  expect_s4_class(one_gene@counts, "dgCMatrix")
  expect_equal(dim(one_gene@counts), c(1L, 3L))
  one_pixel <- restrict_counts(puck, c("g1", "g2"), UMI_max = 20)
  expect_equal(dim(one_pixel@counts), c(2L, 1L))
  expect_equal(colnames(one_pixel@counts), "bc1")
})

test_that("restrict_puck subsets all slots consistently and ignores unknown barcodes", {
  puck <- SpatialRNA(small_coords(), small_counts())
  res <- restrict_puck(puck, c("bc3", "bc1", "missing"))
  expect_equal(colnames(res@counts), c("bc1", "bc3"))
  expect_equal(names(res@nUMI), c("bc1", "bc3"))
  expect_equal(rownames(res@coords), c("bc1", "bc3"))
  expect_equal(res@nUMI, puck@nUMI[c("bc1", "bc3")])
  one <- restrict_puck(puck, "bc2")
  expect_equal(dim(one@counts), c(4L, 1L))
  expect_equal(rownames(one@coords), "bc2")
})

test_that("coerce_old rebuilds a SpatialRNA with the same data", {
  puck <- SpatialRNA(small_coords(), small_counts())
  res <- coerce_old(puck)
  expect_s4_class(res, "SpatialRNA")
  expect_equal(res@counts, puck@counts)
  expect_equal(res@nUMI, puck@nUMI)
  expect_equal(res@coords, puck@coords)
})

test_that("save.SpatialRNA and read.SpatialRNA round-trip", {
  puck <- SpatialRNA(small_coords(), small_counts())
  dir <- file.path(withr::local_tempdir(), "puck")
  save.SpatialRNA(puck, dir)
  expect_true(all(file.exists(file.path(dir, c("coords.csv", "nUMI.csv", "counts.csv")))))

  coords <- read.csv(file.path(dir, "coords.csv"))
  colnames(coords)[1] <- "barcodes"
  write.csv(coords, file.path(dir, "coords.csv"), row.names = FALSE)
  res <- suppressMessages(read.SpatialRNA(dir))
  expect_s4_class(res, "SpatialRNA")
  expect_equal(as.matrix(res@counts), as.matrix(puck@counts))
  expect_equal(res@nUMI, puck@nUMI)
  expect_equal(res@coords$x, puck@coords$x)
})
