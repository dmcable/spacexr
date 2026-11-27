plot_rctd <- function(X2_fun) {
  puck <- make_puck()
  barcodes <- colnames(puck@counts)
  my_beta <- matrix(1, length(barcodes), 1, dimnames = list(barcodes, "TypeA"))
  X2 <- X2_fun(seq(0, 1, length.out = length(barcodes)))
  rownames(X2) <- barcodes
  rctd <- make_rctd_shell(puck)
  rctd@internal_vars_de <- list(barcodes = barcodes, my_beta = my_beta, X2 = X2)
  rctd@de_results <- list(gene_fits = list(mean_val = matrix(0.5, 1, 1, dimnames = list("Gene1", "TypeA"))),
                          sig_gene_list = list(TypeA = NULL))
  rctd
}

test_that("CSIDE gene plots work without ggplot2 attached", {
  skip_if("package:ggplot2" %in% search())
  two <- plot_rctd(function(x) cbind(1, x))
  p <- plot_gene_two_regions(two, "Gene1", "TypeA")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))

  regions <- plot_rctd(function(x) cbind(as.numeric(x <= 0.5), as.numeric(x > 0.5)))
  p2 <- plot_gene_regions(regions, "TypeA", "Gene1")
  expect_s3_class(p2, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p2))
})
