ref_inputs <- function() {
  counts <- matrix(c(100, 200, 50, 300,
                     20, 10, 150, 40,
                     5, 15, 25, 35), nrow = 3, byrow = TRUE,
                   dimnames = list(paste0("g", 1:3), paste0("c", 1:4)))
  cell_types <- factor(c("A", "A", "B", "B"))
  names(cell_types) <- colnames(counts)
  list(counts = counts, cell_types = cell_types)
}

test_that("Reference builds a valid object", {
  d <- ref_inputs()
  ref <- Reference(d$counts, d$cell_types)
  expect_s4_class(ref, "Reference")
  expect_s4_class(ref@counts, "dgCMatrix")
  expect_equal(ref@nUMI, colSums(d$counts))
  expect_equal(ref@cell_types, d$cell_types)
  expect_equal(levels(ref@cell_types), c("A", "B"))
})

test_that("Reference removes cells below min_UMI", {
  d <- ref_inputs()
  # nUMI = c(125, 225, 225, 375)
  ref <- Reference(d$counts, d$cell_types, min_UMI = 200)
  expect_equal(colnames(ref@counts), c("c2", "c3", "c4"))
  expect_equal(names(ref@cell_types), c("c2", "c3", "c4"))
  expect_error(Reference(d$counts, d$cell_types, min_UMI = 1000),
               "no barcodes were included")
})

test_that("Reference warns about cell types with no remaining cells", {
  d <- ref_inputs()
  ct <- factor(c("A", "A", "B", "B"), levels = c("A", "B", "C"))
  names(ct) <- colnames(d$counts)
  expect_warning(ref <- Reference(d$counts, ct), "missing cell types with no occurences:  C")
  expect_equal(levels(ref@cell_types), c("A", "B", "C"))
})

test_that("Reference keeps only mutually shared barcodes", {
  d <- ref_inputs()
  ct <- d$cell_types[1:3]
  expect_warning(ref <- Reference(d$counts, ct), "not mutually shared")
  expect_equal(colnames(ref@counts), c("c1", "c2", "c3"))
  ct2 <- d$cell_types; names(ct2) <- paste0("x", 1:4)
  expect_error(Reference(d$counts, ct2), "do not share any barcode names")
})

test_that("Reference respects a provided nUMI", {
  d <- ref_inputs()
  nUMI <- colSums(d$counts) + 1
  expect_warning(ref <- Reference(d$counts, d$cell_types, nUMI = nUMI), "nUMI does not match")
  expect_equal(ref@nUMI, nUMI)
})

test_that("Reference downsamples to n_max_cells per cell type", {
  d <- make_reference_data(n_cells_per_type = 30)
  expect_warning(ref <- withr::with_seed(1, Reference(d$counts, d$cell_types, n_max_cells = 10)),
                 "Downsampling number of cells to: 10")
  expect_true(all(table(ref@cell_types) == 10))
  expect_equal(ncol(ref@counts), 30)
  expect_equal(names(ref@nUMI), colnames(ref@counts))
  expect_equal(names(ref@cell_types), colnames(ref@counts))
})

test_that("Reference requires at least two cells", {
  d <- ref_inputs()
  expect_error(Reference(d$counts[, 1, drop = FALSE], d$cell_types[1]), "only one cell present")
})

test_that("check_cell_types validates cell type factors", {
  expect_error(check_cell_types(c(a = "A", b = "B")), "not a factor")
  expect_error(check_cell_types(factor(c(a = "A"))), "length\\(cell_types\\) < 2")
  expect_error(check_cell_types(factor(c(a = "A", b = "A"))), "multiple levels")
  expect_error(check_cell_types(factor(c("A", "B"))), "names\\(cell_types\\) is null")
  expect_error(check_cell_types(factor(c(a = "A", b = ""))), "empty name")
  expect_error(check_cell_types(factor(c(a = "A", b = "B/C"))), "prohibited character /")
  expect_silent(check_cell_types(factor(c(a = "A", b = "B"))))
})

test_that("restrict_reference subsets all slots", {
  d <- ref_inputs()
  ref <- Reference(d$counts, d$cell_types)
  res <- restrict_reference(ref, c("c4", "c1"))
  expect_equal(colnames(res@counts), c("c4", "c1"))
  expect_equal(names(res@nUMI), c("c4", "c1"))
  expect_equal(as.character(res@cell_types), c("B", "A"))
})

test_that("restrict_reference_cell_types keeps selected types and drops levels", {
  ref <- make_reference()
  res <- restrict_reference_cell_types(ref, c("TypeA", "TypeC"))
  expect_equal(levels(res@cell_types), c("TypeA", "TypeC"))
  expect_true(all(res@cell_types %in% c("TypeA", "TypeC")))
  expect_equal(ncol(res@counts), 60)
})

test_that("create_downsampled_data samples cells per type", {
  ref <- make_reference()
  res <- withr::with_seed(3, create_downsampled_data(ref, cell_types_keep = c("TypeA", "TypeB"), n_samples = 5))
  expect_equal(as.vector(table(res@cell_types)), c(5, 5))
  expect_equal(levels(res@cell_types), c("TypeA", "TypeB"))
  expect_equal(names(res@cell_types), colnames(res@counts))
  big <- create_downsampled_data(ref, n_samples = 100)
  expect_equal(ncol(big@counts), ncol(ref@counts))
})

test_that("coerce_deglam_reference rebuilds a Reference with all cells", {
  d <- ref_inputs()
  ref <- Reference(d$counts, d$cell_types)
  res <- coerce_deglam_reference(ref)
  expect_s4_class(res, "Reference")
  expect_equal(res@counts, ref@counts)
})

test_that("save.Reference writes counts, nUMI and cell types", {
  d <- ref_inputs()
  ref <- Reference(d$counts, d$cell_types)
  dir <- file.path(withr::local_tempdir(), "ref")
  save.Reference(ref, dir)
  saved_counts <- read.csv(file.path(dir, "counts.csv"), row.names = 1)
  expect_equal(as.matrix(saved_counts), d$counts)
  saved_types <- read.csv(file.path(dir, "cell_types.csv"), row.names = 1)
  expect_equal(saved_types[, 1], as.character(d$cell_types))
})

test_that("remap_celltypes maps cluster identifiers to names", {
  path <- withr::local_tempfile(fileext = ".csv")
  write.csv(data.frame(Cluster = c(1, 2, 3), Name = c("A", "B", "A")), path, row.names = FALSE)
  res <- remap_celltypes(path, c(3, 1, 2))
  expect_equal(as.character(res), c("A", "A", "B"))
})
