test_that("get_cell_type_info averages UMI-normalized profiles per type", {
  counts <- matrix(c(1, 3, 2, 2,
                     3, 1, 8, 0,
                     0, 0, 0, 8), nrow = 3, byrow = TRUE,
                   dimnames = list(c("g1", "g2", "g3"), paste0("c", 1:4)))
  cell_types <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))
  nUMI <- colSums(counts)
  info <- get_cell_type_info(counts, cell_types, nUMI)
  expect_length(info, 3)
  expect_equal(info[[2]], c("A", "B"))
  expect_equal(info[[3]], 2)
  expect_equal(info[[1]]$A, c(0.5, 0.5, 0))
  expect_equal(info[[1]]$B, c(0.2, 0.4, 0.4))
  expect_equal(rownames(info[[1]]), c("g1", "g2", "g3"))
  expect_equal(unname(colSums(info[[1]])), c(1, 1))
})

test_that("get_cell_type_info respects cell_type_names ordering", {
  d <- make_reference_data()
  info <- get_cell_type_info(d$counts, d$cell_types, colSums(d$counts), cell_type_names = c("TypeC", "TypeA"))
  expect_equal(colnames(info[[1]]), c("TypeC", "TypeA"))
  expect_equal(info[[3]], 2)
})

test_that("get_cell_type_info recovers simulated profiles", {
  d <- make_reference_data()
  info <- get_cell_type_info(d$counts, d$cell_types, colSums(d$counts))
  expect_lt(max(abs(as.matrix(info[[1]]) - d$profiles)), 0.005)
})

test_that("get_norm_ref matches the bulk expression of the target data", {
  puck <- make_puck()
  P <- as.data.frame(make_profiles())
  genes <- rownames(P)[1:50]
  props <- c(TypeA = 1, TypeB = 1, TypeC = 1)
  res <- get_norm_ref(puck, P, genes, props)
  expect_equal(dim(res), c(50L, 3L))
  bulk <- rowSums(puck@counts)[genes] / sum(puck@nUMI)
  weighted <- as.matrix(res) %*% (props / sum(props))
  expect_equal(as.vector(weighted), unname(bulk))
})

test_that("get_de_genes finds marker genes", {
  puck <- make_puck()
  d <- make_reference_data()
  info <- get_cell_type_info(d$counts, d$cell_types, colSums(d$counts))
  genes <- suppressMessages(get_de_genes(info, puck, fc_thresh = 0.75, expr_thresh = 0.0002))
  expect_setequal(genes, paste0("Gene", 1:90))
})

test_that("get_de_genes removes mitochondrial genes", {
  puck <- make_puck()
  d <- make_reference_data()
  rownames(d$counts)[1] <- "mt-Gene1"
  info <- get_cell_type_info(d$counts, d$cell_types, colSums(d$counts))
  counts <- puck@counts; rownames(counts)[1] <- "mt-Gene1"
  puck2 <- SpatialRNA(puck@coords, counts)
  genes <- suppressMessages(get_de_genes(info, puck2, fc_thresh = 0.75, expr_thresh = 0.0002))
  expect_false("mt-Gene1" %in% genes)
  expect_length(genes, 89)
})

test_that("get_de_genes handles two cell types", {
  puck <- make_puck()
  d <- make_reference_data()
  keep <- d$cell_types != "TypeC"
  ct <- droplevels(d$cell_types[keep])
  info <- get_cell_type_info(d$counts[, keep], ct, colSums(d$counts[, keep]))
  genes <- suppressMessages(get_de_genes(info, puck, fc_thresh = 0.75, expr_thresh = 0.0002))
  expect_setequal(genes, paste0("Gene", 1:60))
})

test_that("get_de_genes errors when genes do not overlap", {
  puck <- make_puck()
  d <- make_reference_data()
  rownames(d$counts) <- paste0("Other", seq_len(nrow(d$counts)))
  info <- get_cell_type_info(d$counts, d$cell_types, colSums(d$counts))
  expect_error(get_de_genes(info, puck), "0 common genes")
  rownames(d$counts)[1:5] <- paste0("Gene", 1:5)
  info <- get_cell_type_info(d$counts, d$cell_types, colSums(d$counts))
  expect_error(get_de_genes(info, puck), "At least 90% of genes do not match")
})

test_that("prepareBulkData scales profiles to total UMI and filters low genes", {
  puck <- make_puck()
  P <- as.data.frame(make_profiles())
  counts <- puck@counts
  low <- names(which(rowSums(counts) < 1000))[1:2]
  res <- prepareBulkData(P, puck, rownames(P), MIN_OBS = 1000)
  expect_named(res, c("X", "b"))
  expect_false(any(low %in% rownames(res$X)))
  expect_equal(res$b, rowSums(counts)[rownames(res$X)])
  expect_equal(res$X, as.matrix(P[rownames(res$X), ] * sum(puck@nUMI)))
})
