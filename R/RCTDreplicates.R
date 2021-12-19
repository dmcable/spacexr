create.RCTD.replicates <- function(spatialRNA.replicates, reference, replicate_names, group_ids = NULL, max_cores = 4, test_mode = FALSE,
                                   gene_cutoff = 0.000125, fc_cutoff = 0.5, gene_cutoff_reg = 0.0002,
                                   fc_cutoff_reg = 0.75, UMI_min = 100, UMI_max = 20000000, UMI_min_sigma = 300,
                        class_df = NULL, CELL_MIN_INSTANCE = 25, cell_type_names = NULL, MAX_MULTI_TYPES = 4,
                        keep_reference = F) {
  if(is.null(cell_type_names))
    cell_type_names <- levels(reference@cell_types)
  cell_type_info <- list(info = process_cell_type_info(reference, cell_type_names = cell_type_names,
                                                         CELL_MIN = CELL_MIN_INSTANCE), renorm = NULL)
  if(class(spatialRNA.replicates) != 'list' ||
     any(!unlist(lapply(spatialRNA.replicates, function(x) class(x) == 'SpatialRNA'))))
    stop('create.RCTD.replicates: spatialRNA.replicates must be a list of SpatialRNA objects.')
  if(length(spatialRNA.replicates) <= 1)
    stop('create.RCTD.replicates: length(spatialRNA.replicates) <= 1. This object must be a list of at least two SpatialRNA objects.')
  if(is.null(group_ids))
    group_ids <- rep(1, length(spatialRNA.replicates))
  if(length(group_ids) != length(replicate_names))
    stop('create.RCTD.replicates: group_ids and replicate_names must both be the same length as the total number of replciates.')
  if(length(group_ids) != length(spatialRNA.replicates))
    stop('create.RCTD.replicates: group_ids must be the same length as the total number of replciates.')
  names(group_ids) <- replicate_names
  check_vector(group_ids, 'group_ids','create.RCTD.replicates', require_int = T)
  if(min(table(group_ids)) < 2)
    stop('create.RCTD.replicates: each group in group_ids must contain at least two replicates')
  RCTD.reps <- list()
  for(i in 1:length(spatialRNA.replicates)) {
    print(paste('create.RCTD.replicates: creating RCTD for replicate',i))
    RCTD.reps[[i]] <- create.RCTD(spatialRNA.replicates[[i]], reference, max_cores = max_cores, test_mode = test_mode,
                                        gene_cutoff = gene_cutoff, fc_cutoff = fc_cutoff, gene_cutoff_reg = gene_cutoff_reg,
                                        fc_cutoff_reg = fc_cutoff_reg, UMI_min = UMI_min, UMI_max = UMI_max, UMI_min_sigma = UMI_min_sigma,
                                        class_df = class_df, CELL_MIN_INSTANCE = CELL_MIN_INSTANCE, cell_type_names = cell_type_names, MAX_MULTI_TYPES = MAX_MULTI_TYPES,
                                        cell_type_info = cell_type_info, keep_reference = F)
  }
  new("RCTD.replicates", RCTD.reps = RCTD.reps, group_ids = group_ids)
}

run.RCTD.replicates <- function(RCTD.replicates, doublet_mode = "doublet") {
  if(!(doublet_mode %in% c('doublet','multi','full')))
    stop(paste0("run.RCTD.replicates: doublet_mode=",doublet_mode, " is not a valid choice. Please set doublet_mode=doublet, multi, or full."))
  for(i in 1:length(RCTD.replicates@RCTD.reps)) {
    print(paste('run.RCTD.replicates: running RCTD for replicate',i))
    RCTD.replicates@RCTD.reps[[i]] <- run.RCTD(RCTD.replicates@RCTD.reps[[i]], doublet_mode = doublet_mode)
  }
  return(RCTD.replicates)
}

run.RCTDE.replicates <- function(RCTD.replicates, explanatory.variable.replicates, cell_types, cell_type_threshold = 125,
                                 gene_threshold = 5e-5, doublet_mode = T, test_mode = 'individual', thresh_val = NULL,
                                 sigma_gene = T, PRECISION.THRESHOLD = 0.01, cell_types_present = NULL, fdr = .01, population_de = T) {
  if(is.null(cell_types))
    stop('run.RCTDE.replicates: cell_types must not be null.')
  if(class(explanatory.variable.replicates) != 'list')
    stop('run.RCTDE.replicates: explanatory.variable.replicates must be a list of explanatory variable vectors for each replicate.')
  if(length(RCTD.replicates@RCTD.reps) != length(explanatory.variable.replicates))
    stop('create.RCTD.replicates: length(explanatory.variable.replicates) is not equal to the number of RCTD replicates, as required.')
  for(i in 1:length(RCTD.replicates@RCTD.reps)) {
    print(paste('run.RCTDE.replicates: running RCTDE for replicate',i))
    RCTD.replicates@RCTD.reps[[i]] <- run.RCTDE.single(
      RCTD.replicates@RCTD.reps[[i]], explanatory.variable.replicates[[i]], cell_types = cell_types, cell_type_threshold = cell_type_threshold,
                                                    gene_threshold = gene_threshold, doublet_mode = doublet_mode, test_mode = test_mode, thresh_val = thresh_val,
                                                    sigma_gene = sigma_gene, PRECISION.THRESHOLD = PRECISION.THRESHOLD, cell_types_present = cell_types_present, fdr = fdr)

  }
  if(population_de)
    RCTD.replicates <- RCTDE.population.inference(RCTD.replicates)
  return(RCTD.replicates)
}

merge.RCTD.objects <- function(RCTD.reps, replicate_names, group_ids = NULL) {
  if(class(RCTD.reps) != 'list' || any(!unlist(lapply(RCTD.reps, function(x) class(x) == 'RCTD'))))
    stop('merge.RCTD.objects: RCTD.reps must be a list of RCTD objects.')
  if(length(RCTD.reps) <= 1)
    stop('merge.RCTD.objects: length(RCTD.replicates) <= 1. This object must be a list of at least two RCTD objects.')
  if(is.null(group_ids))
    group_ids <- rep(1, length(RCTD.reps))
  if(length(group_ids) != length(replicate_names))
    stop('merge.RCTD.objects: group_ids and replicate_names must both be the same length as the total number of replciates.')
  if(length(group_ids) != length(RCTD.reps))
    stop('merge.RCTD.objects: group_ids must be the same length as the total number of replciates.')
  names(group_ids) <- replicate_names
  check_vector(group_ids, 'group_ids','create.RCTD.replicates', require_int = T)
  if(min(table(group_ids)) < 2)
    stop('create.RCTD.replicates: each group in group_ids must contain at least two replicates')
  new("RCTD.replicates", RCTD.reps = RCTD.reps, group_ids = groups_ids)
}

RCTDE.population.inference <- function(RCTD.replicates, use.groups = FALSE) {
  print(paste0('RCTDE.population.inference: running population DE inference with use.groups=', use.groups))
  RCTDde_list <- RCTD.replicates@RCTD.reps
  de_results_list <- lapply(RCTDde_list, function(x) x@de_results)
  myRCTD <- RCTDde_list[[1]]
  cell_types <- myRCTD@internal_vars_de$cell_types
  cell_types_present <- myRCTD@internal_vars_de$cell_types_present
  de_pop_all <- list()
  gene_final_all <- list()
  final_df <- list()
  for(cell_type in cell_types) {
    res <- one_ct_genes(cell_type, RCTDde_list, de_results_list, NULL, cell_types_present,
                        plot_results = F, use.groups = use.groups, group_ids = RCTD.replicates@group_ids)
    de_pop_all[[cell_type]] <- res$de_pop
    gene_final_all[[cell_type]] <- res$gene_final
    final_df[[cell_type]] <- res$final_df
  }
  RCTD.replicates@population_de_results <- de_pop_all
  RCTD.replicates@population_sig_gene_list <- gene_final_all
  RCTD.replicates@population_sig_gene_df <- final_df
  return(RCTD.replicates)
}

save.RCTDE.replicates <- function(RCTD.replicates, resultsdir) {
  if(!dir.exists(resultsdir))
    dir.create(resultsdir)
  myRCTD <- RCTD.replicates@RCTD.reps[[1]]
  cell_types <- myRCTD@internal_vars_de$cell_types
  for(cell_type in cell_types) {
    write.csv(RCTD.replicates@population_sig_gene_df[[cell_type]],
              file.path(resultsdir,paste0(cell_type,'_cell_type_genes.csv')))
  }
}
