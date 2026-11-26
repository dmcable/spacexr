# Synthetic data generators and shared fixtures for the spacexr test suite.

test_cell_types <- c("TypeA", "TypeB", "TypeC")

# Genes x cell types matrix of normalized expression profiles (columns sum to 1).
# Each cell type has `n_markers` marker genes that are 20x enriched.
make_profiles <- function(n_markers = 30, n_background = 60) {
  n_types <- length(test_cell_types)
  n_genes <- n_types * n_markers + n_background
  genes <- paste0("Gene", seq_len(n_genes))
  P <- matrix(1, nrow = n_genes, ncol = n_types, dimnames = list(genes, test_cell_types))
  for (k in seq_len(n_types))
    P[((k - 1) * n_markers + 1):(k * n_markers), k] <- 20
  sweep(P, 2, colSums(P), "/")
}

make_reference_data <- function(n_cells_per_type = 30, umi_range = c(800, 1500), seed = 1) {
  withr::with_seed(seed, {
    P <- make_profiles()
    n_cells <- n_cells_per_type * ncol(P)
    types <- rep(colnames(P), each = n_cells_per_type)
    umi <- sample(umi_range[1]:umi_range[2], n_cells, replace = TRUE)
    counts <- vapply(seq_len(n_cells), function(i) rpois(nrow(P), umi[i] * P[, types[i]]),
                     numeric(nrow(P)))
    barcodes <- paste0("cell", seq_len(n_cells))
    dimnames(counts) <- list(rownames(P), barcodes)
    cell_types <- factor(types, levels = colnames(P))
    names(cell_types) <- barcodes
    list(counts = counts, cell_types = cell_types, profiles = P)
  })
}

make_reference <- function(...) {
  d <- make_reference_data(...)
  suppressWarnings(Reference(d$counts, d$cell_types))
}

# Simulated spatial data: singlet pixels of each type followed by A/B doublets.
make_puck_data <- function(n_singlets_per_type = 40, n_doublets = 30, umi_range = c(600, 1200), seed = 2) {
  withr::with_seed(seed, {
    P <- make_profiles()
    types <- colnames(P)
    W <- rbind(
      do.call(rbind, lapply(seq_along(types), function(k) {
        w <- matrix(0, n_singlets_per_type, length(types)); w[, k] <- 1; w
      })),
      matrix(rep(c(0.5, 0.5, 0), each = n_doublets), ncol = length(types))
    )
    colnames(W) <- types
    n_pix <- nrow(W)
    umi <- sample(umi_range[1]:umi_range[2], n_pix, replace = TRUE)
    counts <- vapply(seq_len(n_pix), function(i) rpois(nrow(P), umi[i] * (P %*% W[i, ])[, 1]),
                     numeric(nrow(P)))
    barcodes <- paste0("pix", seq_len(n_pix))
    dimnames(counts) <- list(rownames(P), barcodes)
    rownames(W) <- barcodes
    coords <- data.frame(x = runif(n_pix, 0, 100), y = runif(n_pix, 0, 100), row.names = barcodes)
    truth <- c(rep(types, each = n_singlets_per_type), rep("doublet", n_doublets))
    names(truth) <- barcodes
    list(counts = counts, coords = coords, weights = W, truth = truth, profiles = P)
  })
}

make_puck <- function(...) {
  d <- make_puck_data(...)
  SpatialRNA(d$coords, d$counts)
}

# Small, fully-specified matrices for constructor/validation tests.
small_counts <- function(n_genes = 4, n_cols = 3, prefix = "bc") {
  m <- matrix(seq_len(n_genes * n_cols), nrow = n_genes, ncol = n_cols)
  dimnames(m) <- list(paste0("g", seq_len(n_genes)), paste0(prefix, seq_len(n_cols)))
  m
}

small_coords <- function(n = 3, prefix = "bc") {
  data.frame(x = seq_len(n), y = seq_len(n) * 2, row.names = paste0(prefix, seq_len(n)))
}

# Loads (and caches) the precomputed likelihood tables shipped with the package.
.likelihood_cache <- new.env()

get_test_X_vals <- function() {
  if (is.null(.likelihood_cache$X_vals))
    .likelihood_cache$X_vals <- readRDS(system.file("extdata", "Qmat/X_vals.rds", package = "spacexr"))
  .likelihood_cache$X_vals
}

get_test_Q_mat <- function(sigma = "100") {
  key <- paste0("Q_", sigma)
  if (is.null(.likelihood_cache[[key]])) {
    for (i in 1:5) {
      Q_list <- readRDS(system.file("extdata", paste0("Qmat/Q_mat_", i, ".rds"), package = "spacexr"))
      if (sigma %in% names(Q_list)) {
        .likelihood_cache[[key]] <- Q_list[[sigma]]
        break
      }
    }
  }
  .likelihood_cache[[key]]
}

use_test_likelihood <- function(sigma = "100") {
  set_likelihood_vars(get_test_Q_mat(sigma), get_test_X_vals())
}

# A minimal RCTD object holding the given SpatialRNA (and optional results).
make_rctd_shell <- function(puck, results = list(), config = list(), internal_vars = list(),
                            cell_type_info = list(), reference = NULL) {
  if (is.null(reference)) {
    ct <- factor(c(c1 = "A", c2 = "B"))
    reference <- Reference(matrix(100, 2, 2, dimnames = list(c("g1", "g2"), c("c1", "c2"))), ct)
  }
  new("RCTD", spatialRNA = puck, originalSpatialRNA = puck, reference = reference, config = config,
      results = results, internal_vars = internal_vars, cell_type_info = cell_type_info)
}
