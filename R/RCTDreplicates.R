#' Creates an \code{\linkS4class{RCTD.replicates}} object across multiple \code{\linkS4class{SpatialRNA}} replicates
#'
#' Applies the \code{\link{create.RCTD}} function for each \code{\linkS4class{SpatialRNA}} replicate inputted using a
#' scRNA-seq reference \code{Reference} object.
#'
#' @param spatialRNA.replicates a list of multiple \code{\linkS4class{SpatialRNA}} objects to run RCTD on
#' @param reference a \code{\linkS4class{Reference}} object scRNA-seq reference used for RCTD
#' @param replicate_names a \code{character} vector of names for each replicate provided in \code{spatialRNA.replicates}
#' @param group_ids (default constant across replicates) a named integer vector (length number of replicates) containing the group id for each replicate.
#' Names represent the replicate names, and replicates of the same group will be expected to be more similar than
#' replicates across groups
#' @param gene_cutoff minimum normalized gene expression for genes to be included in the platform effect normalization step.
#' @param fc_cutoff minimum log-fold-change (across cell types) for genes to be included in the platform effect normalization step.
#' @param gene_cutoff_reg minimum normalized gene expression for genes to be included in the RCTD step.
#' @param fc_cutoff_reg minimum log-fold-change (across cell types) for genes to be included in the RCTD step.
#' @param UMI_min minimum UMI per pixel included in the analysis
#' @param UMI_max maximum UMI per pixel included in the analysis
#' @param UMI_min_sigma minimum UMI per pixel for the \link{choose_sigma_c} function
#' @param max_cores for parallel processing, the number of cores used. If set to 1, parallel processing is not used. The system will additionally be checked for
#' number of available cores.
#' @param class_df (optional) if not NULL, then a dataframe mapping each cell type to a cell class, so that RCTD will report confidence on the class level.
#' @param CELL_MIN_INSTANCE minimum number of cells required per cell type. Default 25, can be lowered if desired.
#' @param cell_type_names A list of cell types to be included from the reference. If NULL, uses all cell types
#' @param MAX_MULTI_TYPES (multi-mode only) Default 4, max number of cell types per pixel
#' @param cell_type_info Default NULL, option to pass in \code{cell_type_info} directly
#' @param keep_reference (Default FALSE) if true, keeps the \code{reference} object stored within the \code{\linkS4class{RCTD}} object
#' @return an \code{\linkS4class{RCTD.replicates}} object, which is ready to run the \code{\link{run.RCTD.replicates}} function
#' @export
create.RCTD.replicates <- function(spatialRNA.replicates, reference, replicate_names, group_ids = NULL, max_cores = 4, test_mode = FALSE,
                                   gene_cutoff = 0.000125, fc_cutoff = 0.5, gene_cutoff_reg = 0.0002,
                                   fc_cutoff_reg = 0.75, UMI_min = 100, UMI_max = 20000000, UMI_min_sigma = 300,
                        class_df = NULL, CELL_MIN_INSTANCE = 25, cell_type_names = NULL, MAX_MULTI_TYPES = 4,
                        keep_reference = F) {
  if(is.null(cell_type_names))
    cell_type_names <- levels(reference@cell_types)
  cell_type_info <- process_cell_type_info(reference, cell_type_names = cell_type_names,
                                                         CELL_MIN = CELL_MIN_INSTANCE)
  if(class(spatialRNA.replicates) != 'list' ||
     any(!unlist(lapply(spatialRNA.replicates, function(x) class(x) == 'SpatialRNA'))))
    stop('create.RCTD.replicates: spatialRNA.replicates must be a list of SpatialRNA objects.')
  if(length(spatialRNA.replicates) <= 1)
    stop('create.RCTD.replicates: length(spatialRNA.replicates) <= 1. This object must be a list of at least two SpatialRNA objects.')
  if(is.null(group_ids))
    group_ids <- rep(1, length(spatialRNA.replicates))
  if(length(group_ids) != length(replicate_names))
    stop('create.RCTD.replicates: group_ids and replicate_names must both be the same length as the total number of replicates.')
  if(length(group_ids) != length(spatialRNA.replicates))
    stop('create.RCTD.replicates: group_ids must be the same length as the total number of replicates.')
  names(group_ids) <- replicate_names
  check_vector(group_ids, 'group_ids','create.RCTD.replicates', require_int = T)
  if(min(table(group_ids)) < 2)
    stop('create.RCTD.replicates: each group in group_ids must contain at least two replicates.')
  RCTD.reps <- list()
  for(i in 1:length(spatialRNA.replicates)) {
    message(paste('create.RCTD.replicates: creating RCTD for replicate',i))
    RCTD.reps[[i]] <- create.RCTD(spatialRNA.replicates[[i]], reference, max_cores = max_cores, test_mode = test_mode,
                                        gene_cutoff = gene_cutoff, fc_cutoff = fc_cutoff, gene_cutoff_reg = gene_cutoff_reg,
                                        fc_cutoff_reg = fc_cutoff_reg, UMI_min = UMI_min, UMI_max = UMI_max, UMI_min_sigma = UMI_min_sigma,
                                        class_df = class_df, CELL_MIN_INSTANCE = CELL_MIN_INSTANCE, cell_type_names = cell_type_names, MAX_MULTI_TYPES = MAX_MULTI_TYPES,
                                        cell_type_profiles = cell_type_info[[1]], keep_reference = F)
  }
  new("RCTD.replicates", RCTD.reps = RCTD.reps, group_ids = group_ids)
}

#' Runs the RCTD pipeline on a \code{\linkS4class{RCTD.replicates}} object
#'
#' For each \code{\linkS4class{SpatialRNA}} replicate in the \code{RCTD.replicates} object, runs the \code{\link{run.RCTD}}
#' function to assign cell types.
#'
#' If in doublet mode, fits at most two cell types per pixel. It classifies each pixel as 'singlet' or 'doublet' and searches for the cell types
#' on the pixel. If in full mode, can fit any number of cell types on each pixel. In multi mode, cell types are added using a greedy algorithm,
#' up to a fixed number.
#'
#' @param RCTD.replicates an \code{\linkS4class{RCTD.replicates}} object created using the \code{\link{create.RCTD.replicates}} function.
#' @param doublet_mode \code{character string}, either "doublet", "multi", or "full" on which mode to run RCTD. Please see above description.
#' @return an \code{\linkS4class{RCTD.replicates}} object containing the results of the RCTD algorithm. Please see \code{\linkS4class{RCTD.replicates}}
#' and \code{\linkS4class{RCTD}} documentation for more information on interpreting the content of this object.
#' @export
run.RCTD.replicates <- function(RCTD.replicates, doublet_mode = "doublet") {
  if(!(doublet_mode %in% c('doublet','multi','full')))
    stop(paste0("run.RCTD.replicates: doublet_mode=",doublet_mode, " is not a valid choice. Please set doublet_mode=doublet, multi, or full."))
  for(i in 1:length(RCTD.replicates@RCTD.reps)) {
    message(paste('run.RCTD.replicates: running RCTD for replicate',i))
    RCTD.replicates@RCTD.reps[[i]] <- run.RCTD(RCTD.replicates@RCTD.reps[[i]], doublet_mode = doublet_mode)
  }
  return(RCTD.replicates)
}

#' Runs CSIDE on a \code{\linkS4class{RCTD.replicates}} object
#'
#' Identifies cell type specific differential expression (DE) as a function of the explanatory variable
#' for each replicate. The design matrix contains an intercept column and a column of the explanatory variable. Uses maximum
#' likelihood estimation to estimate DE and standard errors for each gene and each cell type. Selects
#' genes with significant nonzero DE. Note: a minimum of three replicates are required for population mode.
#'
#' @param RCTD.replicates an \code{\linkS4class{RCTD.replicates}} object with annotated cell types e.g. from the \code{\link{run.RCTD.replicates}} function.
#' @param cell_types the cell types used for CSIDE. Each cell type must occur
#' at least `cell_type_threshold`, as aggregated by \code{\link{aggregate_cell_types}}
#' @param explanatory.variable.replicates (only used for de_mode = single) a list of the named numeric vectors representing for each replicate the explanatory variable used for explaining differential expression in CSIDE.
#' Names of the vectors are the \code{\linkS4class{SpatialRNA}} pixel names, and values should be standardized between 0 and 1.
#' @param X.replicates (only used for de_mode = general) a list for each replicate of matrices containing the covariates for running CSIDE. The rownames represent pixel names and
#' should be a subset of the pixels in the \code{\linkS4class{SpatialRNA}} object. The columns each represent a covariate for
#' explaining differential expression and need to be linearly independent.
#' @param cell_type_threshold (default 125) min occurence of number of cells for each cell type to be used, as aggregated by \code{\link{aggregate_cell_types}}
#' @param gene_threshold (default 5e-5) minimum average normalized expression required for selecting genes
#' @param doublet_mode (default TRUE) if TRUE, uses RCTD doublet mode weights. Otherwise, uses RCTD full mode weights
#' @param sigma_gene (default TRUE) if TRUE, fits gene specific overdispersion parameter. If FALSE, overdispersion parameter is same across all genes.
#' @param weight_threshold (default NULL) the threshold of total normalized weights across all cell types
#' in \code{cell_types} per pixel to be included in the model. Default 0.99 for doublet_mode or 0.95 for full_mode.
#' @param PRECISION.THRESHOLD (default 0.05) for checking for convergence, the maximum parameter change per algorithm step
#' @param cell_types_present cell types (a superset of `cell_types`) to be considered as occuring often enough
#' to consider for gene expression contamination during the step filtering out marker genes of other cell types.
#' @param fdr (default 0.01) false discovery rate for hypothesis testing
#' @param normalize_expr (default FALSE) if TRUE, constrains total gene expression to sum to 1 in each condition.
#' @param population_de whether population-level DE should be run (can also be run later using the \code{\link{CSIDE.population.inference}} function.)
#' @param replicate_index (default all replicates) integer list of replicate indices (subset of 1:N_replicates) to be run for CSIDE
#' @param de_mode (default 'single', otherwise 'nonparam' or 'general') if 'single', calls \code{\link{run.CSIDE.single}}.
#' If 'nonparam', calls \code{\link{run.CSIDE.nonparam}}. If 'general', calls \code{\link{run.CSIDE}}.
#' @param df (default 15) for de_mode = nonparam, the degrees of freedom, or number of basis functions to be used in the model.
#' @param log_fc_thresh (default 0.4) the natural log fold change cutoff for differential expression
#' @param test_error (default FALSE) if TRUE, first tests for error messages before running CSIDE.
#' If set to TRUE, this can be used to quickly evaluate if CSIDE will run without error.
#' @param test_genes_sig_individual (default FALSE) logical controlling whether on individual samples genes will be tested for significance.
#' @param params_to_test: (default 2 for test_mode = 'individual', all parameters for test_mode = 'categorical'). An integer vector of parameter
#' indices to test. For example c(1,4,5) would test only parameters corresponding to columns 1, 4, and 5 of the design matrix.
#' @param barcodes for de_mode = nonparam, the barcodes, or pixel names, of the \code{\linkS4class{SpatialRNA}} object to be used when fitting the model.
#' @return an \code{\linkS4class{RCTD.replicates}} object containing the results of the CSIDE algorithm. See \code{\linkS4class{RCTD.replicates}}
#' for documentation on the \code{population_de_results}, \code{population_sig_gene_list}, and \code{population_sig_gene_df} objects.
#' @export
run.CSIDE.replicates <- function(RCTD.replicates, cell_types, explanatory.variable.replicates = NULL, X.replicates = NULL, cell_type_threshold = 125,
                                 gene_threshold = 5e-5, doublet_mode = T, weight_threshold = NULL,
                                 sigma_gene = T, PRECISION.THRESHOLD = 0.05, cell_types_present = NULL,
                                 fdr = .01, population_de = T, replicate_index = NULL, normalize_expr = F, test_genes_sig_individual = F,
                                 de_mode = 'single', df = 15, barcodes = NULL, log_fc_thresh = 0.4, test_error = F,
                                 params_to_test = NULL, test_mode = 'individual') {
  if(!(de_mode %in% c('single','nonparam', 'general')))
    stop('run.CISDE.replicates: de_mode must be set to "single", "general", or "nonparam".')
  if(is.null(cell_types))
    stop('run.CSIDE.replicates: cell_types must not be null.')
  if(is.null(replicate_index))
    replicate_index <- 1:length(RCTD.replicates@RCTD.reps)
  if(any(!(replicate_index %in% 1:length(RCTD.replicates@RCTD.reps))))
    stop('run.CSIDE.replicates: replicate_index must be a subest of 1:N_replicates')
  if(test_error)
    warning('run.CSIDE.replicates: test_error is TRUE, so this run will just test C-SIDE for errors without running C-SIDE.')

  for(i in replicate_index) {
    if(test_error)
      message(paste('run.CSIDE.replicates: testing CSIDE for errors for replicate',i))
    else
      message(paste('run.CSIDE.replicates: running CSIDE for replicate',i))
    if(de_mode == 'single') {
      if(is.null(explanatory.variable.replicates))
        stop('run.CSIDE.replicates: if de_mode = single, explanatory.variable.replicates cannot be null.')
      if(class(explanatory.variable.replicates) != 'list')
        stop('run.CSIDE.replicates: explanatory.variable.replicates must be a list of explanatory variable vectors for each replicate.')
      if(length(RCTD.replicates@RCTD.reps) != length(explanatory.variable.replicates))
        stop('create.RCTD.replicates: length(explanatory.variable.replicates) is not equal to the number of RCTD replicates, as required.')
      RCTD.replicates@RCTD.reps[[i]] <- run.CSIDE.single(
        RCTD.replicates@RCTD.reps[[i]], explanatory.variable.replicates[[i]], cell_types = cell_types, cell_type_threshold = cell_type_threshold,
        gene_threshold = gene_threshold, doublet_mode = doublet_mode, weight_threshold = weight_threshold,
        sigma_gene = sigma_gene, PRECISION.THRESHOLD = PRECISION.THRESHOLD, test_genes_sig = test_genes_sig_individual,
        cell_types_present = cell_types_present, fdr = fdr, log_fc_thresh = log_fc_thresh, test_error = test_error)
    } else if(de_mode == 'nonparam') {
      RCTD.replicates@RCTD.reps[[i]] <- run.CSIDE.nonparam(
        RCTD.replicates@RCTD.reps[[i]], df = df, barcodes = barcodes, cell_types = cell_types, cell_type_threshold = cell_type_threshold,
        gene_threshold = gene_threshold, doublet_mode = doublet_mode, weight_threshold = weight_threshold, test_genes_sig = test_genes_sig_individual,
        sigma_gene = sigma_gene, PRECISION.THRESHOLD = PRECISION.THRESHOLD, cell_types_present = cell_types_present, fdr = fdr, test_error = test_error)
    } else {
      if(is.null(X.replicates))
        stop('run.CSIDE.replicates: if de_mode = single, X.replicates cannot be null.')
      if(class(X.replicates) != 'list')
          stop('run.CSIDE.replicates: X.replicates must be a list of design matrices for each replicate.')
      if(length(RCTD.replicates@RCTD.reps) != length(X.replicates))
        stop('run.CSIDE.replicates: length(X.replicates) is not equal to the number of RCTD replicates, as required.')
      X <- X.replicates[[i]]
      if(length(setdiff(rownames(X),rownames(RCTD_controls@RCTD.reps[[i]]@results$weights))) > 0)
        warning('run.CSIDE.replicates: some elements of rownames(X.replicates) do not appear in myRCTD object (myRCTD@results$weights) for this replicate, but they are required to be a subset.')
      RCTD.replicates@RCTD.reps[[i]] <- run.CSIDE(
        RCTD.replicates@RCTD.reps[[i]], X, rownames(X), cell_types = cell_types,
        cell_type_threshold = cell_type_threshold, gene_threshold = gene_threshold, doublet_mode = doublet_mode, test_genes_sig = test_genes_sig_individual,
        weight_threshold = weight_threshold, sigma_gene = sigma_gene, PRECISION.THRESHOLD = PRECISION.THRESHOLD,
        cell_types_present = cell_types_present, fdr = fdr, log_fc_thresh = log_fc_thresh, test_error = test_error,
        params_to_test = params_to_test)
    }
  }
  if(population_de)
    RCTD.replicates <- CSIDE.population.inference(RCTD.replicates, log_fc_thresh = log_fc_thresh)
  return(RCTD.replicates)
}

#' Creates an \code{\linkS4class{RCTD.replicates}} object across multiple \code{\linkS4class{RCTD}} objects
#'
#' @param RCTD.reps a list of multiple \code{\linkS4class{RCTD}} objects to merge into one \code{\linkS4class{RCTD.replicates}} object.
#' @param replicate_names a \code{character} vector of names for each replicate provided in \code{RCTD.reps}
#' @param group_ids (default constant across replicates) a named integer vector (length number of replicates) containing the group id for each replicate.
#' Names represent the replicate names, and replicates of the same group will be expected to be more similar than
#' replicates across groups
#' @return an \code{\linkS4class{RCTD.replicates}} object, containing each \code{\linkS4class{RCTD}} object in \code{RCTD.reps}
#' @export
merge.RCTD.objects <- function(RCTD.reps, replicate_names, group_ids = NULL) {
  if(class(RCTD.reps) != 'list' || any(!unlist(lapply(RCTD.reps, function(x) class(x) == 'RCTD'))))
    stop('merge.RCTD.objects: RCTD.reps must be a list of RCTD objects.')
  if(length(RCTD.reps) <= 1)
    stop('merge.RCTD.objects: length(RCTD.replicates) <= 1. This object must be a list of at least two RCTD objects.')
  if(is.null(group_ids))
    group_ids <- rep(1, length(RCTD.reps))
  if(length(group_ids) != length(replicate_names))
    stop('merge.RCTD.objects: group_ids and replicate_names must both be the same length as the total number of replicates.')
  if(length(group_ids) != length(RCTD.reps))
    stop('merge.RCTD.objects: group_ids must be the same length as the total number of replicates.')
  names(group_ids) <- replicate_names
  check_vector(group_ids, 'group_ids','create.RCTD.replicates', require_int = T)
  if(min(table(group_ids)) < 2)
    stop('create.RCTD.replicates: each group in group_ids must contain at least two replicates.')
  new("RCTD.replicates", RCTD.reps = RCTD.reps, group_ids = group_ids)
}

#' Runs population-level differential expression inference for a \code{\linkS4class{RCTD.replicates}} object
#'
#' First, CSIDE must have been run on all replicates using e.g. the \code{\link{run.CSIDE.replicates}} function.
#'
#' @param RCTD.replicates a \code{\linkS4class{RCTD.replicates}} object for which to perform population-level DE inference. Note, at least three
#' replicates must be provided.
#' @param use.groups (default FALSE) if TRUE, treats the replicates as having multiple groups (e.g. samples) according to the \code{group_ids} slot
#' @param MIN.CONV.REPLICATES (default 2) the minimum number of replicates (if not use.groups) for which a gene must converge
#' @param MIN.CONV.GROUPS (default 2) the minimum number of groups (if use.groups) for which a gene must converge
#' @param CT.PROP (default 0.5) minimum ratio of gene expression within cell type compared to other cell types
#' @param q_thresh (default 0.01) false discovery rate
#' @param log_fc_thresh (default 0.4) minimum natural log estimated DE threshold
#' @param params_to_test: (default 2 for test_mode = 'individual', all parameters for test_mode = 'categorical'). An integer vector of parameter
#' indices to test. Note, for population mode, only the first parameter is tested.
#' @param normalize_expr (default FALSE) if TRUE, constrains total gene expression to sum to 1 in each condition
#' @return an \code{\linkS4class{RCTD.replicates}} object containing the results of the CSIDE population-level algorithm. See \code{\linkS4class{RCTD.replicates}}
#' for documentation on the \code{population_de_results}, \code{population_sig_gene_list}, and \code{population_sig_gene_df} objects.
#' @export
CSIDE.population.inference <- function(RCTD.replicates, params_to_test = NULL, use.groups = FALSE, MIN.CONV.REPLICATES = 2,
                                        MIN.CONV.GROUPS = 2, CT.PROP = 0.5,
                                       q_thresh = 0.01, log_fc_thresh = 0.4,
                                       normalize_expr = F) {
  message(paste0('CSIDE.population.inference: running population DE inference with use.groups=', use.groups))
  MIN.REPS <- 3
  if(length(RCTD.replicates@RCTD.reps) < MIN.REPS)
    stop('CSIDE.population.inference: minimum of three replicates required for population mode.')
  RCTDde_list <- RCTD.replicates@RCTD.reps
  myRCTD <- RCTDde_list[[1]]
  if(is.null(params_to_test))
    params_to_test <- myRCTD@internal_vars_de$params_to_test[1]
  cell_types <- myRCTD@internal_vars_de$cell_types
  cell_types_present <- myRCTD@internal_vars_de$cell_types_present
  de_pop_all <- list()
  gene_final_all <- list()
  final_df <- list()
  for(i in 1:length(RCTDde_list)) {
    RCTDde_list[[i]] <- normalize_de_estimates(RCTDde_list[[i]], normalize_expr)
  }
  de_results_list <- lapply(RCTDde_list, function(x) x@de_results)
  for(cell_type in cell_types) {
    ct_pres <- sapply(RCTDde_list, function(x) cell_type %in% x@internal_vars_de$cell_types)
    if(sum(ct_pres) >= MIN.REPS) {
      res <- one_ct_genes(cell_type, RCTDde_list[ct_pres], de_results_list[ct_pres], NULL, cell_types_present, params_to_test,
                          plot_results = F, use.groups = use.groups,
                          group_ids = RCTD.replicates@group_ids, MIN.CONV.REPLICATES = MIN.CONV.REPLICATES,
                          MIN.CONV.GROUPS = MIN.CONV.GROUPS, CT.PROP = CT.PROP,
                          q_thresh = q_thresh, log_fc_thresh = log_fc_thresh, normalize_expr = normalize_expr)
      de_pop_all[[cell_type]] <- res$de_pop
      gene_final_all[[cell_type]] <- res$gene_final
      final_df[[cell_type]] <- res$final_df
    } else {
      warning(paste('CSIDE.population.inference: cell type', cell_type,
                     'was removed from population-level analysis because it was run on fewer than the minimum required three replicates.'))
    }
  }
  RCTD.replicates@population_de_results <- de_pop_all
  RCTD.replicates@population_sig_gene_list <- gene_final_all
  RCTD.replicates@population_sig_gene_df <- final_df
  return(RCTD.replicates)
}

#' Saves the CSIDE population-level differential expression results for a \code{\linkS4class{RCTD.replicates}} object
#'
#' First, CSIDE must have been run on all replicates at the population level using e.g. the \code{\link{run.CSIDE.replicates}} function.
#'
#' @param RCTD.replicates a \code{\linkS4class{RCTD.replicates}} object containing population-level DE inference results.
#' @param resultsdir a directory where to save the significant gene matrices for each cell type.
#' @export
save.CSIDE.replicates <- function(RCTD.replicates, resultsdir) {
  if(!dir.exists(resultsdir))
    dir.create(resultsdir)
  myRCTD <- RCTD.replicates@RCTD.reps[[1]]
  cell_types <- myRCTD@internal_vars_de$cell_types
  for(cell_type in cell_types) {
    write.csv(RCTD.replicates@population_de_results[[cell_type]],
              file.path(resultsdir,paste0(cell_type,'_cell_type_genes_all.csv')))
    write.csv(RCTD.replicates@population_sig_gene_df[[cell_type]],
              file.path(resultsdir,paste0(cell_type,'_cell_type_genes_sig.csv')))
  }
}
