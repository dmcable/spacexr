#' Runs CSIDE on a \code{\linkS4class{RCTD}} object with a single explanatory variable
#'
#' Identifies cell type specific differential expression (DE) as a function of the explanatory variable.
#' The design matrix contains an intercept column and a column of the explanatory variable. Uses maximum
#' likelihood estimation to estimate DE and standard errors for each gene and each cell type. Selects
#' genes with significant nonzero DE.
#'
#' @param myRCTD an \code{\linkS4class{RCTD}} object with annotated cell types e.g. from the \code{\link{run.RCTD}} function.
#' @param explanatory.variable a named numeric vector representing the explanatory variable used for explaining differential expression in CSIDE. Names of the variable
#' are the \code{\linkS4class{SpatialRNA}} pixel names, and values should be standardized between 0 and 1.
#' @param cell_types the cell types used for CSIDE. If null, cell types will be chosen with aggregate occurrences of
#' at least `cell_type_threshold`, as aggregated by \code{\link{choose_cell_types}}
#' @param cell_type_threshold (default 125) min occurrence of number of cells for each cell type to be used, as aggregated by \code{\link{choose_cell_types}}
#' @param gene_threshold (default 5e-5) minimum average normalized expression required for selecting genes
#' @param doublet_mode (default TRUE) if TRUE, uses RCTD doublet mode weights. Otherwise, uses RCTD full mode weights
#' @param sigma_gene (default TRUE) if TRUE, fits gene specific overdispersion parameter. If FALSE, overdispersion parameter is same across all genes.
#' @param weight_threshold (default NULL) the threshold of total normalized weights across all cell types
#' in \code{cell_types} per pixel to be included in the model. Default 0.99 for doublet_mode or 0.8 for full_mode.
#' @param PRECISION.THRESHOLD (default 0.01) for checking for convergence, the maximum parameter change per algorithm step
#' @param cell_types_present cell types (a superset of `cell_types`) to be considered as occurring often enough
#' to consider for gene expression contamination during the step filtering out marker genes of other cell types.
#' @param fdr (default 0.01) false discovery rate for hypothesis testing
#' @param test_genes_sig (default TRUE) logical controlling whether genes will be tested for significance
#' @param normalize_expr (default FALSE) if TRUE, constrains total gene expression to sum to 1 in each condition.
#' @param logs (default FALSE) if TRUE, writes progress to logs/de_logs.txt
#' @return an \code{\linkS4class{RCTD}} object containing the results of the CSIDE algorithm. Contains objects \code{de_results},
#' which contain the results of the CSIDE algorithm including `gene_fits`, which contains the results of fits on individual genes,
#' in addition `sig_gene_list`, a list, for each cell type, of significant genes detected by CSIDE.
#' Additionally, the object contains `internal_vars_de` a list of variables that are used internally by CSIDE
#' @export
run.CSIDE.single <- function(myRCTD, explanatory.variable,  cell_types = NULL, cell_type_threshold = 125,
                          gene_threshold = 5e-5, doublet_mode = T, weight_threshold = NULL,
                          sigma_gene = T, PRECISION.THRESHOLD = 0.01, cell_types_present = NULL, fdr = .01,
                          test_genes_sig = T, normalize_expr = F, logs=F) {
  X2 <- build.designmatrix.single(myRCTD, explanatory.variable)
  barcodes <- rownames(X2)
  return(run.CSIDE(myRCTD, X2, barcodes, cell_types, gene_threshold = gene_threshold, cell_type_threshold = cell_type_threshold,
                       doublet_mode = doublet_mode, test_mode = 'individual', params_to_test = 2,
                       weight_threshold = weight_threshold, sigma_gene = sigma_gene, test_genes_sig = test_genes_sig,
                       PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                       cell_types_present = cell_types_present, fdr = fdr, normalize_expr = normalize_expr,
                   logs=logs))
}

#' Runs CSIDE on a \code{\linkS4class{RCTD}} object to detect nonparametric smooth gene expression patterns
#'
#' Identifies cell type specific smooth gene expression patterns. The design matrix contains thin plate spline
#' basis functions spanning the space of smooth functions. Uses maximum likelihood estimation to estimate
#' DE and standard errors for each gene and each cell type. Selects genes with significant nonzero DE.
#'
#' @param myRCTD an \code{\linkS4class{RCTD}} object with annotated cell types e.g. from the \code{\link{run.RCTD}} function.
#' @param df (default 15) the degrees of freedom, or number of basis functions to be used in the model.
#' @param barcodes the barcodes, or pixel names, of the \code{\linkS4class{SpatialRNA}} object to be used when fitting the model.
#' @param cell_types the cell types used for CSIDE. If null, cell types will be chosen with aggregate occurences of
#' at least `cell_type_threshold`, as aggregated by \code{\link{choose_cell_types}}
#' @param cell_type_threshold (default 125) min occurence of number of cells for each cell type to be used, as aggregated by \code{\link{choose_cell_types}}
#' @param gene_threshold (default 5e-5) minimum average normalized expression required for selecting genes
#' @param doublet_mode (default TRUE) if TRUE, uses RCTD doublet mode weights. Otherwise, uses RCTD full mode weights
#' @param sigma_gene (default TRUE) if TRUE, fits gene specific overdispersion parameter. If FALSE, overdispersion parameter is same across all genes.
#' @param weight_threshold (default NULL) the threshold of total normalized weights across all cell types
#' in \code{cell_types} per pixel to be included in the model. Default 0.99 for doublet_mode or 0.8 for full_mode.
#' @param PRECISION.THRESHOLD (default 0.01) for checking for convergence, the maximum parameter change per algorithm step
#' @param cell_types_present cell types (a superset of `cell_types`) to be considered as occuring often enough
#' to consider for gene expression contamination during the step filtering out marker genes of other cell types.
#' @param fdr (default 0.01) false discovery rate for hypothesis testing
#' @param test_genes_sig (default TRUE) logical controlling whether genes will be tested for significance
#' @return an \code{\linkS4class{RCTD}} object containing the results of the CSIDE algorithm. Contains objects \code{de_results},
#' which contain the results of the CSIDE algorithm including `gene_fits`, which contains the results of fits on individual genes,
#' in addition `sig_gene_list`, a list, for each cell type, of significant genes detected by CSIDE.
#' Additionally, the object contains `internal_vars_de` a list of variables that are used internally by CSIDE
#' @param logs (default FALSE) if TRUE, writes progress to logs/de_logs.txt
#' @export
run.CSIDE.nonparam <- function(myRCTD, df = 15, barcodes = NULL, cell_types = NULL,
                            cell_type_threshold = 125, gene_threshold = 5e-5, doublet_mode = T,
                            weight_threshold = NULL, sigma_gene = T,
                            PRECISION.THRESHOLD = 0.01, cell_types_present = NULL, fdr = .01, test_genes_sig = T,
                            logs=F) {
  X2 <- build.designmatrix.nonparam(myRCTD, barcodes = barcodes, df = df)
  barcodes <- rownames(X2)
  return(run.CSIDE(myRCTD, X2, barcodes, cell_types, gene_threshold = gene_threshold,
                       doublet_mode = doublet_mode, test_mode = 'individual', cell_type_threshold = cell_type_threshold,
                       weight_threshold = weight_threshold, sigma_gene = sigma_gene,test_genes_sig = test_genes_sig,
                       PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                       cell_types_present = cell_types_present, params_to_test = 2:df, fdr = fdr, normalize_expr = F,
                   logs=logs))
}

#' Runs CSIDE on a \code{\linkS4class{RCTD}} object for DE across multiple discrete regions
#'
#' Identifies cell type specific differential expression (DE) across multiple discrete regions
#' The design matrix contains for each region a column of 0s and 1s representing membership in that region. Uses maximum
#' likelihood estimation to estimate DE and standard errors for each gene and each cell type. Selects
#' genes with significant nonzero DE. Tests for differences in gene expression across regions.
#'
#' @param myRCTD an \code{\linkS4class{RCTD}} object with annotated cell types e.g. from the \code{\link{run.RCTD}} function.
#' @param region_list a list of \code{character} vectors, where each vector contains pixel names, or barcodes, for a single region. These pixel names
#' should be a subset of the pixels in the \code{\linkS4class{SpatialRNA}} object
#' @param cell_types the cell types used for CSIDE. If null, cell types will be chosen with aggregate occurences of
#' at least `cell_type_threshold`, as aggregated by \code{\link{choose_cell_types}}
#' @param cell_type_threshold (default 125) min occurence of number of cells for each cell type to be used, as aggregated by \code{\link{choose_cell_types}}
#' @param gene_threshold (default 5e-5) minimum average normalized expression required for selecting genes
#' @param doublet_mode (default TRUE) if TRUE, uses RCTD doublet mode weights. Otherwise, uses RCTD full mode weights
#' @param sigma_gene (default TRUE) if TRUE, fits gene specific overdispersion parameter. If FALSE, overdispersion parameter is same across all genes.
#' @param weight_threshold (default NULL) the threshold of total normalized weights across all cell types
#' in \code{cell_types} per pixel to be included in the model. Default 0.99 for doublet_mode or 0.8 for full_mode.
#' @param PRECISION.THRESHOLD (default 0.01) for checking for convergence, the maximum parameter change per algorithm step
#' @param cell_types_present cell types (a superset of `cell_types`) to be considered as occuring often enough
#' to consider for gene expression contamination during the step filtering out marker genes of other cell types.
#' @param fdr (default 0.01) false discovery rate for hypothesis testing
#' @param test_genes_sig (default TRUE) logical controlling whether genes will be tested for significance
#' @param logs (default FALSE) if TRUE, writes progress to logs/de_logs.txt
#' @return an \code{\linkS4class{RCTD}} object containing the results of the CSIDE algorithm. Contains objects \code{de_results},
#' which contain the results of the CSIDE algorithm including `gene_fits`, which contains the results of fits on individual genes,
#' in addition `sig_gene_list`, a list, for each cell type, of significant genes detected by CSIDE.
#' Additionally, the object contains `internal_vars_de` a list of variables that are used internally by CSIDE
#' @export
run.CSIDE.regions <- function(myRCTD, region_list, cell_types = NULL,
                           cell_type_threshold = 125, gene_threshold = 5e-5, doublet_mode = T,
                          weight_threshold = NULL, sigma_gene = T,
                           PRECISION.THRESHOLD = 0.01, cell_types_present = NULL, fdr = 0.01, test_genes_sig = T,
                          logs=F) {
  X2 <- build.designmatrix.regions(myRCTD, region_list)
  barcodes <- rownames(X2)
  return(run.CSIDE(myRCTD, X2, barcodes, cell_types, cell_type_threshold = cell_type_threshold, gene_threshold = gene_threshold,
                       doublet_mode = doublet_mode, test_mode = 'categorical',
                       weight_threshold = weight_threshold, sigma_gene = sigma_gene, params_to_test = 1:dim(X2)[2],
                       PRECISION.THRESHOLD = PRECISION.THRESHOLD,test_genes_sig = test_genes_sig,
                       cell_types_present = cell_types_present, fdr = fdr, normalize_expr = F,
                   logs=logs))
}

#' Runs cell type specific CSIDE on a \code{\linkS4class{RCTD}} object with a general design matrix
#'
#' Identifies cell type specific differential expression (DE) across a general design matrix of covariates. Uses maximum
#' likelihood estimation to estimate DE and standard errors for each gene and each cell type. Selects
#' genes with significant nonzero DE. The type of test is determined by \code{test_mode}, and the parameters tested
#' is determined by \code{params_to_test}.
#'
#' @param myRCTD an \code{\linkS4class{RCTD}} object with annotated cell types e.g. from the \code{\link{run.RCTD}} function.
#' @param X a matrix containing the covariates for running CSIDE. The rownames represent pixel names and
#' should be a subset of the pixels in the \code{\linkS4class{SpatialRNA}} object. The columns each represent a covariate for
#' explaining differential expression and need to be linearly independent.
#' @param barcodes the barcodes, or pixel names, of the \code{\linkS4class{SpatialRNA}} object to be used when fitting the model.
#' @param cell_types the cell types used for CSIDE. If null, cell types will be chosen with aggregate occurences of
#' at least `cell_type_threshold`, as aggregated by \code{\link{choose_cell_types}}
#' @param cell_type_specific: (default TRUE for all covariates). A logical vector of length the number of covariates
#' indicating whether each covariate's DE parameters should be cell type-specific or shared across all cell types.
#' @param params_to_test: (default 2 for test_mode = 'individual', all parameters for test_mode = 'categorical'). An integer vector of parameter
#' indices to test. For example c(1,4,5) would test only parameters corresponding to columns 1, 4, and 5 of the design matrix.
#' @param cell_type_threshold (default 125) min occurence of number of cells for each cell type to be used, as aggregated by \code{\link{choose_cell_types}}
#' @param gene_threshold (default 5e-5) minimum average normalized expression required for selecting genes
#' @param doublet_mode (default TRUE) if TRUE, uses RCTD doublet mode weights. Otherwise, uses RCTD full mode weights
#' @param sigma_gene (default TRUE) if TRUE, fits gene specific overdispersion parameter. If FALSE, overdispersion parameter is same across all genes.
#' @param weight_threshold (default NULL) the threshold of total normalized weights across all cell types
#' in \code{cell_types} per pixel to be included in the model. Default 0.99 for doublet_mode or 0.8 for full_mode.
#' @param test_mode (default 'individual') if 'individual', tests for DE individually for each parameter. If 'categorical', then tests for differences
#' across multiple categorical parameters
#' @param PRECISION.THRESHOLD (default 0.01) for checking for convergence, the maximum parameter change per algorithm step
#' @param cell_types_present cell types (a superset of `cell_types`) to be considered as occuring often enough
#' to consider for gene expression contamination during the step filtering out marker genes of other cell types.
#' @param fdr (default 0.01) false discovery rate for hypothesis testing
#' @param normalize_expr (default FALSE) if TRUE, constrains total gene expression to sum to 1 in each condition.
#' Setting normalize_expr = TRUE is only valid for testing single parameters with test_mode = 'individual'.
#' @param test_genes_sig (default TRUE) logical controlling whether genes will be tested for significance
#' @param logs (default FALSE) if TRUE, writes progress to logs/de_logs.txt
#' @return an \code{\linkS4class{RCTD}} object containing the results of the CSIDE algorithm. Contains objects \code{de_results},
#' which contain the results of the CSIDE algorithm including `gene_fits`, which contains the results of fits on individual genes,
#' in addition `sig_gene_list`, a list, for each cell type, of significant genes detected by CSIDE.
#' Additionally, the object contains `internal_vars_de` a list of variables that are used internally by CSIDE
#' @export
run.CSIDE <- function(myRCTD, X, barcodes, cell_types, gene_threshold = 5e-5, cell_type_threshold = 125,
                          doublet_mode = T, test_mode = 'individual', weight_threshold = NULL,
                          sigma_gene = T, PRECISION.THRESHOLD = 0.01, cell_types_present = NULL,
                          test_genes_sig = T, fdr = .01, cell_type_specific = NULL,
                      params_to_test = NULL, normalize_expr = F, logs=F) {
  X <- check_designmatrix(X, 'run.CSIDE', require_2d = TRUE)
  if(is.null(cell_type_specific))
    cell_type_specific <- !logical(dim(X)[2])
  check_cell_type_specific(cell_type_specific, dim(X)[2])
  X1 <- X[,!cell_type_specific];
  if(any(!cell_type_specific))
    X2 <- X[,cell_type_specific]
  else
    X2 <- X
  return(run.CSIDE.general(myRCTD, X1, X2, barcodes, cell_types, cell_type_threshold = cell_type_threshold,
                           gene_threshold = gene_threshold,
                       doublet_mode = doublet_mode, test_mode = test_mode, weight_threshold = weight_threshold,
                       sigma_gene = sigma_gene, PRECISION.THRESHOLD = PRECISION.THRESHOLD, params_to_test = params_to_test,
                       cell_types_present = cell_types_present, test_genes_sig = test_genes_sig,
                       fdr = fdr, normalize_expr = normalize_expr, logs=logs))
}

#' Runs CSIDE on a \code{\linkS4class{RCTD}} object with a general design matrix
#'
#' Identifies differential expression (DE) across a general design matrix of covariates. DE parameters can be
#' cell type-specific or shared across all cell types. Uses maximum
#' likelihood estimation to estimate DE and standard errors for each gene and each cell type. Selects
#' genes with significant nonzero DE. The type of test is determined by \code{test_mode}, and the parameters tested
#' is determined by \code{params_to_test}.
#'
#' @param myRCTD an \code{\linkS4class{RCTD}} object with annotated cell types e.g. from the \code{\link{run.RCTD}} function.
#' @param X1 a matrix containing the covariates shared across all cell types. The rownames represent pixel names and
#' should be a subset of the pixels in the \code{\linkS4class{SpatialRNA}} object. The columns each represent a covariate for
#' explaining differential expression and need to be linearly independent.
#' @param X2 a matrix containing the cell type-specific covariates. The rownames represent pixel names and
#' should be a subset of the pixels in the \code{\linkS4class{SpatialRNA}} object. The columns each represent a covariate for
#' explaining differential expression and need to be linearly independent.
#' @param barcodes the barcodes, or pixel names, of the \code{\linkS4class{SpatialRNA}} object to be used when fitting the model.
#' @param cell_types the cell types used for CSIDE. If null, cell types will be chosen with aggregate occurences of
#' at least `cell_type_threshold`, as aggregated by \code{\link{choose_cell_types}}
#' @param params_to_test: (default 2 for test_mode = 'individual', all parameters for test_mode = 'categorical'). An integer vector of parameter
#' indices to test. For example c(1,4,5) would test only parameters corresponding to columns 1, 4, and 5 of the design matrix X2.
#' @param cell_type_threshold (default 125) min occurence of number of cells for each cell type to be used, as aggregated by \code{\link{choose_cell_types}}
#' @param gene_threshold (default 5e-5) minimum average normalized expression required for selecting genes
#' @param doublet_mode (default TRUE) if TRUE, uses RCTD doublet mode weights. Otherwise, uses RCTD full mode weights
#' @param sigma_gene (default TRUE) if TRUE, fits gene specific overdispersion parameter. If FALSE, overdispersion parameter is same across all genes.
#' @param weight_threshold (default NULL) the threshold of total normalized weights across all cell types
#' in \code{cell_types} per pixel to be included in the model. Default 0.99 for doublet_mode or 0.8 for full_mode.
#' @param test_mode (default 'individual') if 'individual', tests for DE individually for each parameter. If 'categorical', then tests for differences
#' across multiple categorical parameters
#' @param PRECISION.THRESHOLD (default 0.01) for checking for convergence, the maximum parameter change per algorithm step
#' @param cell_types_present cell types (a superset of `cell_types`) to be considered as occuring often enough
#' to consider for gene expression contamination during the step filtering out marker genes of other cell types.
#' @param fdr (default 0.01) false discovery rate for hypothesis testing
#' @param test_genes_sig (default TRUE) logical controlling whether genes will be tested for significance
#' @param normalize_expr (default FALSE) if TRUE, constrains total gene expression to sum to 1 in each condition.
#' Setting normalize_expr = TRUE is only valid for testing single parameters with test_mode = 'individual'.
#' @param logs (default FALSE) if TRUE, writes progress to logs/de_logs.txt
#' @return an \code{\linkS4class{RCTD}} object containing the results of the CSIDE algorithm. Contains objects \code{de_results},
#' which contain the results of the CSIDE algorithm including `gene_fits`, which contains the results of fits on individual genes,
#' in addition `sig_gene_list`, a list, for each cell type, of significant genes detected by CSIDE, whereas
#' `all_gene_list` is the analogous list for all genes (including nonsignificant).
#' Additionally, the object contains `internal_vars_de` a list of variables that are used internally by CSIDE
#' @export
run.CSIDE.general <- function(myRCTD, X1, X2, barcodes, cell_types, gene_threshold = 5e-5, cell_type_threshold = 125,
                          doublet_mode = T, test_mode = 'individual', weight_threshold = NULL,
                          sigma_gene = T, PRECISION.THRESHOLD = 0.01, cell_types_present = NULL,
                          test_genes_sig = T, fdr = .01, params_to_test = NULL, normalize_expr = F,
                          logs=F) {
  if(doublet_mode && myRCTD@config$RCTDmode != 'doublet')
    stop('run.CSIDE.general: attempted to run CSIDE in doublet mode, but RCTD was not run in doublet mode. Please run CSIDE in full mode (doublet_mode = F) or run RCTD in doublet mode.')
  if(!any("cell_types_assigned" %in% names(myRCTD@internal_vars)) || !myRCTD@internal_vars$cell_types_assigned)
    stop('run.CSIDE.general: cannot run CSIDE unless cell types have been assigned. If cell types have been assigned, you may run "myRCTD <- set_cell_types_assigned(myRCTD)".')
  cell_types <- choose_cell_types(myRCTD, barcodes, doublet_mode, cell_type_threshold, cell_types)
  X1 <- check_designmatrix(X1, 'run.CSIDE.general')
  X2 <- check_designmatrix(X2, 'run.CSIDE.general', require_2d = TRUE)
  if(!(test_mode %in% c('individual', 'categorical')))
    stop(c('run.CSIDE.general: not valid test_mode = ',test_mode,'. Please set test_mode = "categorical" or "individual".'))
  if(is.null(params_to_test))
    if(test_mode == 'individual')
      params_to_test <- 2
    else
      params_to_test <- 1:dim(X2)[2]
  if(normalize_expr && (test_mode != 'individual' || length(params_to_test) > 1))
    stop('run.CSIDE.general: Setting normalize_expr = TRUE is only valid for testing single parameters with test_mode = individual')
  message(paste0("run.CSIDE.general: configure params_to_test = ",params_to_test))
  if(any(!(params_to_test %in% 1:dim(X2)[2])))
    stop(c('run.CSIDE.general: params_to_test must be a vector of integers from 1 to dim(X2)[2] = ', dim(X2)[2],
           'please make sure that tested parameters are in the required range.'))
  if(test_mode == 'categorical' && any(!(X2[,params_to_test] %in% c(0,1))))
    stop(c('run.CSIDE.general: for test_mode = categorical, colums params_to_test, ',params_to_test,', must have values 0 or 1.'))
  if(is.null(cell_types_present))
    cell_types_present <- cell_types
  if(any(!(barcodes %in% rownames(X1))) || any(!(barcodes %in% rownames(X2))))
    stop('run.CSIDE.general: some barcodes do not appear in the rownames of X1 or X2.')
  puck = myRCTD@originalSpatialRNA
  gene_list_tot <- filter_genes(puck, threshold = gene_threshold)
  nUMI <- puck@nUMI[barcodes]
  cell_type_info <- myRCTD@cell_type_info$info
  if(doublet_mode) {
    my_beta <- get_beta_doublet(barcodes, cell_type_info[[2]], myRCTD@results$results_df, myRCTD@results$weights_doublet)
    thresh = 0.999
  } else {
    my_beta <- as.matrix(sweep(myRCTD@results$weights, 1, rowSums(myRCTD@results$weights), '/'))
    thresh = 0.8
  }
  if(!is.null(weight_threshold))
    thresh = weight_threshold
  res <- filter_barcodes_cell_types(barcodes, cell_types, my_beta, thresh = thresh)
  barcodes <- res$barcodes; my_beta <- res$my_beta
  set_likelihood_vars(myRCTD@internal_vars$Q_mat, myRCTD@internal_vars$X_vals)
  if(sigma_gene)
    set_global_Q_all()
  sigma_init <- as.character(100*myRCTD@internal_vars$sigma)
  gene_fits <- get_de_gene_fits(X1[barcodes, , drop = FALSE],X2[barcodes, , drop = FALSE],my_beta, nUMI[barcodes], gene_list_tot,
                                cell_types, restrict_puck(puck, barcodes), barcodes, sigma_init,
                                test_mode, numCores = myRCTD@config$max_cores, sigma_gene = sigma_gene,
                                PRECISION.THRESHOLD = PRECISION.THRESHOLD, params_to_test = params_to_test,
                                logs=logs)
  if(normalize_expr)
    myRCTD <- normalize_de_estimates(myRCTD, normalize_expr = normalize_expr,
                                     param_position = params_to_test)
  if(test_genes_sig) {
    both_gene_list <- get_sig_genes(puck, myRCTD, gene_list_tot, cell_types, my_beta, barcodes, nUMI,
                            gene_fits, cell_types_present, X2, test_mode, fdr = fdr,
                            params_to_test = params_to_test, normalize_expr = normalize_expr)
    sig_gene_list <- both_gene_list$sig_gene_list; all_gene_list <- both_gene_list$all_gene_list
  } else {
    sig_gene_list <- NULL
    all_gene_list <- NULL
  }
  myRCTD@internal_vars_de <- list(barcodes = barcodes, cell_types = cell_types, doublet_mode = doublet_mode,
                                  cell_types_present = cell_types_present,
                                  my_beta = my_beta, X1 = X1, X2 = X2,
                                  test_mode = test_mode, params_to_test = params_to_test)
  myRCTD@de_results <- list(gene_fits = gene_fits, sig_gene_list = sig_gene_list, all_gene_list = all_gene_list)
  return(myRCTD)
}

get_sig_genes <- function(puck, myRCTD, gene_list_tot, cell_types, my_beta, barcodes, nUMI,
                          gene_fits, cell_types_present, X2, test_mode, params_to_test = 2,
                          fdr = .01, p_thresh = 1, log_fc_thresh = 0.4, normalize_expr = F) {
  cti_renorm <- get_norm_ref(puck, myRCTD@cell_type_info$info[[1]], intersect(gene_list_tot,rownames(myRCTD@cell_type_info$info[[1]])), myRCTD@internal_vars$proportions)
  sig_gene_list <- list(); all_gene_list <- list()
  for(cell_type in cell_types) {
    gene_list_type <- get_gene_list_type(my_beta, barcodes, cell_type, nUMI, gene_list_tot,
                                         cti_renorm, cell_types_present, gene_fits, test_mode = test_mode)
    if(test_mode == 'individual')
      both_genes <- find_sig_genes_individual(cell_type, cell_types, gene_fits, gene_list_type, X2,
                                         params_to_test = params_to_test, fdr = fdr, p_thresh = p_thresh,
                                         log_fc_thresh = log_fc_thresh, normalize_expr = normalize_expr)
    else if(test_mode == 'categorical') {
      sig_genes <- find_sig_genes_categorical(cell_type, cell_types, gene_fits, gene_list_type, X2,
                                        p_thresh = p_thresh, log_fc_thresh = log_fc_thresh,
                                        params_to_test = params_to_test)
    }
    sig_genes <- both_genes$sig_genes; all_genes <- both_genes$all_genes
    sig_gene_list[[cell_type]] <- sig_genes
    all_gene_list[[cell_type]] <- all_genes
  }
  return(list(sig_gene_list = sig_gene_list, all_gene_list = all_gene_list))
}

test_genes_sig_post <- function(myRCTD, params_to_test = NULL, fdr = .01, p_thresh = 1,
                                log_fc_thresh = 0.4, normalize_expr = F) {
  puck <- myRCTD@originalSpatialRNA
  gene_list_tot <- rownames(myRCTD@de_results$gene_fits$s_mat)
  cell_types <- myRCTD@internal_vars_de$cell_types
  my_beta <- myRCTD@internal_vars_de$my_beta
  X2 <- myRCTD@internal_vars_de$X2
  barcodes <- myRCTD@internal_vars_de$barcodes
  nUMI <- puck@nUMI[barcodes]
  test_mode <- myRCTD@internal_vars_de$test_mode
  cell_types_present <- myRCTD@internal_vars_de$cell_types_present
  gene_fits <- myRCTD@de_results$gene_fits
  if(is.null(params_to_test))
    params_to_test <- myRCTD@internal_vars_de$params_to_test
  both_gene_list <- get_sig_genes(puck, myRCTD, gene_list_tot, cell_types, my_beta, barcodes, nUMI,
                                  gene_fits, cell_types_present, X2, test_mode,
                                  params_to_test = params_to_test, fdr = fdr,
                                 p_thresh = p_thresh, log_fc_thresh = log_fc_thresh, normalize_expr = normalize_expr)
  myRCTD@de_results$sig_gene_list <- both_gene_list$sig_gene_list
  myRCTD@de_results$all_gene_list <- both_gene_list$all_gene_list
  return(myRCTD)
}

find_sig_genes_categorical <- function(cell_type, cell_types, gene_fits, gene_list_type, X2, fdr = 0.01,
                                p_thresh = 1, log_fc_thresh = 0.4, params_to_test = NULL) {
  if(length(gene_list_type) == 0)
    stop(paste0('find_sig_genes_categorical: cell type ', cell_type,
                ' has not converged on any genes. Consider removing this cell type from the model using the cell_types option.'))
  if(is.null(params_to_test))
    params_to_test <- 1:dim(X2)[2]
  n_regions <- length(params_to_test); n_cell_types <- length(cell_types)
  cell_ind = (which(cell_types == cell_type))
  s_mat_ind <- (1:dim(X2)[2]) + (n_regions*(cell_ind - 1))
  p_val_sig_pair <- numeric(length(gene_list_type)); names(p_val_sig_pair) <- gene_list_type
  log_fc_best_pair <- numeric(length(gene_list_type)); names(log_fc_best_pair) <- gene_list_type
  sd_vec <- numeric(length(gene_list_type)); names(sd_vec) <- gene_list_type
  sd_lfc_vec <- numeric(length(gene_list_type)); names(sd_lfc_vec) <- gene_list_type
  i1_vec <- numeric(length(gene_list_type)); names(i1_vec) <- gene_list_type
  i2_vec <- numeric(length(gene_list_type)); names(i2_vec) <- gene_list_type
  for(gene in gene_list_type) {
    con_regions <- get_con_regions(gene_fits, gene, dim(X2)[2], cell_ind, n_cell_types) &
      (params_to_test %in% 1:dim(X2)[2])
    n_regions_con <- sum(con_regions)
    x <- gene_fits$all_vals[gene, con_regions,cell_ind]
    s_mat_ind_cur <- s_mat_ind[con_regions]
    var_vals <- (gene_fits$s_mat[gene, s_mat_ind_cur])^2
    ovr_best_p_val <- 1
    best_log_fc <- 0; best_sd <- 0
    best_i1 <- 0; best_i2 <- 0
    for(i1 in 1:(n_regions_con-1))
      for(i2 in (i1+1):n_regions_con) {
        log_fc <- abs(x[i1] - x[i2])
        sd_cur <- sqrt(var_vals[i1] + var_vals[i2])
        z_score <- (log_fc) / sd_cur
        p_val <- 2*(pnorm(-z_score))
        if(p_val < ovr_best_p_val) {
          ovr_best_p_val <- p_val
          best_log_fc <- log_fc
          best_sd <- sd_cur
          best_i1 <- which(con_regions)[i1]
          best_i2 <- which(con_regions)[i2]
        }
      }
    p_val_sig_pair[gene] <- min(1, ovr_best_p_val * choose(n_regions_con, 2))
    log_fc_best_pair[gene] <- best_log_fc
    sd_vec[gene] <- best_sd
    sd_lfc_vec[gene] <- sd(x)
    i1_vec[gene] <- best_i1
    i2_vec[gene] <- best_i2
  }
  gene_list_sig <- fdr_sig_genes(gene_list_type, p_val_sig_pair, fdr)
  sig_genes <- data.frame(sd_lfc_vec[gene_list_type], i1_vec[gene_list_type], i2_vec[gene_list_type],
                          sd_vec[gene_list_type], p_val_sig_pair[gene_list_type],
                          log_fc_best_pair[gene_list_type])
  rownames(sig_genes) <- gene_list_type
  custom_names <- c('sd_lfc','sd_best','p_val_best','log_fc_best', 'paramindex1_best', 'paramindex2_best')
  colnames(sig_genes) <- custom_names
  all_genes <- sig_genes
  sig_genes <- sig_genes[abs(sig_genes$p_val < p_thresh) & abs(sig_genes$log_fc) >= log_fc_thresh, ]
  if(length(gene_list_sig) > 1)
    sig_genes <- data.frame(sig_genes, gene_fits$all_vals[rownames(sig_genes),params_to_test,cell_ind]) # add on the means
  else {
    if(length(gene_list_sig) == 1) {
      sig_genes <- data.frame(t(unlist((c(sig_genes, gene_fits$all_vals[rownames(sig_genes),params_to_test,cell_ind],
                                          gene_fits$s_mat[rownames(sig_genes),s_mat_ind[params_to_test]])))))
      rownames(sig_genes) <- gene_list_sig
      if(length(sig_genes) > 0)
        colnames(sig_genes)[(length(custom_names)+1):length(sig_genes)] <-
        c(lapply(params_to_test,function(x) paste('mean_',x)), lapply(params_to_test,function(x) paste('sd_',x)))
    } else
      sig_genes <- list()
  }
  return(list(sig_genes = sig_genes, all_genes = all_genes))
}

find_sig_genes_individual <- function(cell_type, cell_types, gene_fits, gene_list_type, X2, params_to_test = 2, fdr = 0.01, p_thresh = 1,
                                      log_fc_thresh = 0.4, normalize_expr = F) {
  if(length(gene_list_type) == 0)
    stop(paste0('find_sig_genes_individual: cell type ', cell_type,
                ' has not converged on any genes. Consider removing this cell type from the model using the cell_types option.'))
  ct_ind <- which(cell_types == cell_type)
  I_ind = dim(X2)[2]*(ct_ind - 1) + params_to_test
  if(normalize_expr) {
    log_fc <- gene_fits$mean_val_cor[[cell_type]][gene_list_type]
  } else {
    log_fc <- gene_fits$all_vals[gene_list_type,params_to_test, ct_ind]
  }
  s_vec <- gene_fits$s_mat[gene_list_type,I_ind]
  z_score <- abs(log_fc) / s_vec
  p_val <- 2*(pnorm(-z_score))
  if(length(params_to_test) > 1)
    p_val <- pmin(apply(p_val, 1, min)*length(params_to_test),1)
  names(p_val) <- gene_list_type
  gene_list_sig <- fdr_sig_genes(gene_list_type, p_val, fdr)
  if(length(gene_list_sig) > 0)
    p_thresh <- min(p_thresh, max(p_val[gene_list_sig]))
  if(length(params_to_test) > 1) {
    p_val <- 2*(pnorm(-z_score))
    best_mat <- function(gene) {
      index <- which(p_val[gene,]*length(params_to_test) < p_thresh)
      if(length(index) > 0) {
        best_ind <- which.max(abs(z_score[gene,index]))
        lfc <- log_fc[gene,index][best_ind]
        sd <- s_vec[gene,index][best_ind]
        z <- z_score[gene,index][best_ind]
        best_ind <- params_to_test[index[best_ind]]
        return(c(best_ind, lfc,sd,z))
      } else {
        return(c(0,0,0,0))
      }
    }
    best_ind <- function(gene) {
      best_mat(gene)[1]
    }
    best_log_fc <- function(gene) {
      best_mat(gene)[2]
    }
    best_sd <- function(gene) {
      best_mat(gene)[3]
    }
    best_Z <- function(gene) {
      best_mat(gene)[4]
    }
    best_indn <- unlist(lapply(gene_list_type, best_ind))
    names(best_indn) <- gene_list_type
    log_fcn <- unlist(lapply(gene_list_type, best_log_fc))
    names(log_fcn) <- gene_list_type
    z_scoren <- unlist(lapply(gene_list_type, best_Z))
    names(z_scoren) <- gene_list_type
    s_vec <- unlist(lapply(gene_list_type, best_sd))
    names(s_vec) <- gene_list_type
    best_ind <- best_indn; z_score <- z_scoren; log_fc <- log_fcn
    p_val <- pmin(apply(p_val, 1, min)*length(params_to_test),1)
  } else {
    best_ind <- rep(params_to_test, length(gene_list_type))
    names(best_ind) <- gene_list_type
  }
  sig_genes <- data.frame(z_score[gene_list_sig], log_fc[gene_list_sig], s_vec[gene_list_sig], best_ind[gene_list_sig])
  names(sig_genes) <- c('Z_score','log_fc', 'se', 'paramindex_best')
  sig_genes$conv <- gene_fits$con_mat[gene_list_sig, cell_type]
  sig_genes$p_val <- p_val[gene_list_sig]
  sig_genes <- sig_genes[abs(sig_genes$p_val < p_thresh) & abs(sig_genes$log_fc) >= log_fc_thresh, ]
  all_genes <- data.frame(z_score[gene_list_type], log_fc[gene_list_type], s_vec[gene_list_type], best_ind[gene_list_type])
  names(all_genes) <- c('Z_score','log_fc', 'se', 'paramindex_best')
  all_genes$conv <- gene_fits$con_mat[gene_list_type, cell_type]
  all_genes$p_val <- p_val[gene_list_type]
  return(list(sig_genes = sig_genes, all_genes = all_genes))
}

get_de_gene_fits <- function(X1,X2,my_beta, nUMI, gene_list, cell_types, puck, barcodes, sigma_init, test_mode,
                             numCores = 4, sigma_gene = T, PRECISION.THRESHOLD = 0.01, params_to_test = 2, logs=F) {
  results_list <- fit_de_genes(X1,X2,my_beta, nUMI, gene_list, puck, barcodes,
                               sigma_init, test_mode, numCores = numCores,
                               sigma_gene = sigma_gene,
                               PRECISION.THRESHOLD = PRECISION.THRESHOLD, logs = logs)
  N_genes <- length(results_list)
  intercept_val <- matrix(0,nrow = N_genes, ncol = length(cell_types))
  mean_val <- matrix(0,nrow = N_genes, ncol = length(cell_types))
  all_vals <- array(0, dim = c(N_genes, dim(X2)[2],length(cell_types)))
  dimnames(all_vals)[[1]] <- gene_list
  dimnames(all_vals)[[3]] <- cell_types
  con_val <- logical(N_genes)
  ll_val <- numeric(N_genes)
  n_val <- numeric(N_genes)
  sigma_g <- numeric(N_genes)
  names(sigma_g) <- gene_list
  I_val <- list()
  names(n_val) <- gene_list
  names(con_val) <- gene_list
  names(ll_val) <- gene_list
  rownames(mean_val) <- gene_list; colnames(mean_val) <- cell_types
  rownames(intercept_val) <- gene_list; colnames(intercept_val) <- cell_types
  d_vals <- matrix(0,nrow=N_genes,ncol=dim(X2)[2]*length(cell_types))
  s_mat <- matrix(0, nrow = N_genes, ncol = dim(X2)[2]*length(cell_types))
  precision_mat <- matrix(0, nrow = N_genes, ncol = dim(X2)[2]*length(cell_types))
  con_all <- matrix(FALSE, nrow = N_genes, ncol = dim(X2)[2]*length(cell_types))
  con_mat <- matrix(FALSE, nrow = N_genes, ncol = length(cell_types))
  error_mat <- matrix(FALSE, nrow = N_genes, ncol = length(cell_types))
  rownames(precision_mat) <- gene_list; rownames(con_all) <- gene_list
  rownames(s_mat) <- gene_list; rownames(con_mat) <- gene_list; rownames(error_mat) <- gene_list
  colnames(s_mat) <- get_param_names(X1,X2, cell_types)
  colnames(precision_mat) <- get_param_names(X1,X2, cell_types)
  colnames(con_all) <- get_param_names(X1,X2, cell_types)
  colnames(con_mat) <- cell_types
  colnames(error_mat) <- cell_types
  rownames(d_vals) <- gene_list
  for(i in 1:N_genes) {
    sigma_g[i] <- results_list[[i]]$sigma_s_best
    res <- results_list[[i]]$res
    d_vals[i,] <- res$d
    mean_val[i,] <- res$alpha2[params_to_test[1],]
    intercept_val[i,] <- res$alpha2[1,]
    all_vals[i, ,] <- res$alpha2
    con_val[i] <- res$converged
    precision_mat[i,] <- res$precision
    ll_val[i] <- res$log_l
    n_val[i] <- res$n.iter
    I_val[[i]] <- res$I
    s_mat[i,] <- sqrt(diag(I_val[[i]]))
    con_mat[i,] <- res$converged_vec
    con_all[i,] <- res$precision < PRECISION.THRESHOLD
    error_mat[i,] <- res$error_vec
  }
  return(list(mean_val = mean_val, con_val = con_val, ll_val = ll_val, I_val = I_val, s_mat = s_mat,
              n.iter = n_val,d_vals = d_vals, intercept_val = intercept_val, all_vals = all_vals,
              precision_mat = precision_mat, sigma_g = sigma_g, con_mat = con_mat, con_all = con_all, error_mat = error_mat))
}

fit_de_genes <- function(X1,X2,my_beta, nUMI, gene_list, puck, barcodes, sigma_init, test_mode, numCores = 4, sigma_gene = T, PRECISION.THRESHOLD = 0.01,
                         logs=F) {
  results_list <- list()
  if(numCores == 1) {
    for(i in 1:length(gene_list)) {
      message(i)
      gene <- gene_list[i]
      Y <- puck@counts[gene, barcodes]
      results_list[[i]] <- estimate_gene_wrapper(Y,X1,X2,my_beta, nUMI, sigma_init, test_mode, verbose = F, n.iter = 200, MIN_CHANGE = 1e-3, sigma_gene = sigma_gene, PRECISION.THRESHOLD = PRECISION.THRESHOLD)
    }
  } else {
    cl <- parallel::makeCluster(numCores,setup_strategy = "sequential",outfile="") #makeForkCluster
    doParallel::registerDoParallel(cl)
    environ = c('estimate_effects_trust', 'solveIRWLS.effects_trust', 'Q_mat', 'K_val','X_vals',
                'calc_log_l_vec', 'calc_Q_k','get_d1_d2', 'calc_Q_all','psd','construct_hess_fast',
                'choose_sigma_gene', 'estimate_gene_wrapper', 'check_converged_vec')
    if(sigma_gene)
      environ <- c(environ, 'Q_mat_all')
    if (logs) {
      out_file = "logs/de_log.txt"
      if(!dir.exists('logs'))
        dir.create('logs')
      if(file.exists(out_file))
        file.remove(out_file)
    }
    results_list <- foreach::foreach(i = 1:length(gene_list), .packages = c("quadprog", "spacexr"), .export = environ) %dopar% {
      if (logs) {
        if(i %% 10 == 0) {
          cat(paste0("Finished sample: ",i," gene ", gene_list[i],"\n"), file=out_file, append=TRUE)
        }
      }
      assign("Q_mat",Q_mat, envir = globalenv()); assign("X_vals",X_vals, envir = globalenv())
      assign("K_val",K_val, envir = globalenv());
      if(sigma_gene)
        assign("Q_mat_all",Q_mat_all, envir = globalenv());
      gene <- gene_list[i]
      Y <- puck@counts[gene, barcodes]
      res <- estimate_gene_wrapper(Y,X1,X2,my_beta, nUMI, sigma_init, test_mode, verbose = F, n.iter = 200, MIN_CHANGE = 1e-3, sigma_gene = sigma_gene)
    }
    parallel::stopCluster(cl)
  }
  return(results_list)
}
