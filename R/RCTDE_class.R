check_designmatrix <- function(X, f_name, require_2d = FALSE) {
  tryCatch({
    X <- as(X,'matrix')
  }, error = function(e) {
    stop(paste0(f_name,': could not convert X to matrix using as(X,\'matrix\'). Please check that X is coercible to matrix, such as a matrix, dgCmatrix, or data.frame.'))
  })
  if(dim(X)[1] == 1) #check more than one gene
    stop(paste0(f_name,': the first dimension of X is 1, indicating only one pixel present. Please format X so that the first dimension is greater than 1.'))
  if(dim(X)[2] == 0 && require_2d)
    stop(paste0(f_name,': the second dimension of X is 0, as no covariates are present. Please format X so that the second dimension is at least 1.'))
  if(dim(X)[2] > 0 && !is.numeric(X[1,1]))
    stop(paste0(f_name,': elements of X are not numeric'))
  if(is.null(rownames(X)))
    stop(paste0(f_name,': rownames(X) is null. Please enter pixel names (i.e. barcodes) as rownames'))
  if(any(duplicated(rownames(X))))
    stop(paste0(f_name,': rownames(X) contain duplicated elements. Please ensure rownames are unique'))
  if(dim(X)[2] > 0) {
    if(Matrix::rankMatrix(X) < dim(X)[2])
      stop(paste0(f_name,': X is not full rank. Please ensure that columns are linearly independent.'))
    intercept <- !any(X[,1] != 1)
    if(intercept)
      start_ind <- 2
    else
      start_ind <- 1
    if(dim(X)[2] > start_ind) {
      for(i in start_ind:dim(X)[2]) {
        if(max(X[,i]) - min(X[,i]) < 1e-9)
          stop(paste0(f_name,': column, ',i, ' of X is a constant. Please ensure that contant term (the intercept) appears',
                      'as a vector of ones in the first column of X.'))
        X[,i] <- X[,i] - min(X[,i])
        X[,i] <- X[,i] / max(X[,i]) # standardize
      }
    }
  }
  return(X)
}

check_cell_type_specific <- function(cell_type_specific, D, f_name) {
  if(!is.atomic(cell_type_specific))
    stop(paste0(f_name,': cell_type_specific is not an atomic vector. Please format cell_type_specific as an atomic vector.'))
  if(!is.logical(cell_type_specific))
    stop(paste0(f_name,': cell_type_specific is not numeric'))
  if(length(cell_type_specific) != D)
    stop(paste0(f_name,': the length of nUMI is not currently equal to dim(X)[2], the number of covariates.'))
}

#' Constructs a design matrix for running RCTDE with a single explanatory variable
#'
#' The design matrix contains an intercept column and a column of the explanatory variable.
#'
#' @param myRCTD an \code{\linkS4class{RCTD}} object with annotated cell types e.g. from the \code{\link{run.RCTD}} function.
#' @param explanatory.variable a named numeric vector representing the explanatory variable used for explaining differential expression in RCTDE. Names of the variable
#' are the \code{\linkS4class{SpatialRNA}} pixel names, and values should be standardized between 0 and 1.
#' @return A design matrix containing the covariates for running RCTDE. The rownames represent pixel names and
#' are a subset of the pixels in the \code{\linkS4class{SpatialRNA}} object. The columns each represent a covariate for
#' explaining differential expression.
#' @export
build.designmatrix.single <- function(myRCTD, explanatory.variable) {
  check_vector(explanatory.variable, 'explanatory.variable', 'build.designmatrix.single')
  barcodes <- intersect(names(explanatory.variable), colnames(myRCTD@spatialRNA@counts))
  if(length(barcodes) <= 1)
    stop(paste0('build.designmatrix.single: ', length(barcodes),
                ' common barcode names found between explanatory.variable and myRCTD@spatialRNA. Please ensure that more common barcodes are found'))
  explanatory.variable <- explanatory.variable[barcodes]
  if(max(explanatory.variable) - min(explanatory.variable) < 1e-9)
    stop('build.designmatrix.single: values of explanatory.variable are constant. Please increase the range of this variable.')
  explanatory.variable <- explanatory.variable - min(explanatory.variable)
  explanatory.variable <- explanatory.variable / max(explanatory.variable) # standardize
  X <- cbind(1,explanatory.variable)
  rownames(X) <- barcodes;
  return(X)
}

#' Constructs a design matrix for running RCTDE nonparametrically
#'
#' The design matrix contains thin plate spline basis functions spanning the space of smooth functions.
#'
#' @param myRCTD an \code{\linkS4class{RCTD}} object with annotated cell types e.g. from the \code{\link{run.RCTD}} function.
#' @param df (default 15) the degrees of freedom, or number of basis functions to be used in the design matrix
#' @param barcodes (default NULL) the barcodes, or pixel names, of the \code{\linkS4class{SpatialRNA}} object to be used when creating the design matrix.
#' @return A design matrix containing the covariates for running RCTDE. The rownames represent pixel names and
#' are a subset of the pixels in the \code{\linkS4class{SpatialRNA}} object. The columns each represent a covariate for
#' explaining differential expression.
#' @export
build.designmatrix.nonparam <- function(myRCTD, barcodes = NULL, df = 15) {
  if(is.null(barcodes))
    barcodes <- colnames(myRCTD@spatialRNA@counts)
  else
    barcodes = intersect(barcodes, colnames(myRCTD@spatialRNA@counts))
  if(length(barcodes) <= 1)
    stop(paste0('build.designmatrix.nonparam: ', length(barcodes),
                ' common barcode names found between barcodes and myRCTD@spatialRNA. Please ensure that more common barcodes are found'))
  X2 <- get_spline_matrix(myRCTD@spatialRNA, df = df)
  rownames(X2) <- barcodes;
  return(X2)
}

#' Constructs a design matrix for running RCTDE across a set of regions
#'
#' The design matrix contains for each region a column of 0s and 1s representing membership in that region.
#'
#' @param myRCTD an \code{\linkS4class{RCTD}} object with annotated cell types e.g. from the \code{\link{run.RCTD}} function.
#' @param region_list a list of \code{character} vectors, where each vector contains pixel names, or barcodes, for a single region. These pixel names
#' should be a subset of the pixels in the \code{\linkS4class{SpatialRNA}} object
#' @return A design matrix containing the covariates for running RCTDE. The rownames represent pixel names and
#' are a subset of the pixels in the \code{\linkS4class{SpatialRNA}} object. The columns each represent a covariate for
#' explaining differential expression.
#' @export
build.designmatrix.regions <- function(myRCTD, region_list) {
  if(!is.list(region_list))
    stop('run.de.regions: error, region_list must be a list')
  n_regions <- length(region_list)
  if(n_regions < 3)
    stop('run.de.regions: length(region_list) <= 2. Must be at least 3 to continue.')
  for(i in 1:n_regions) {
    barcodes <- region_list[[i]]
    if(!is.character(barcodes) || !is.atomic(barcodes))
      stop('run.de.regions: error, region_list must be a list of atomic character vectors')
    shorter_barcodes <- intersect(barcodes, colnames(myRCTD@spatialRNA@counts))
    if(length(barcodes) > length(shorter_barcodes))
      warning('run.de.regions: some barcodes from region_list not found in RCTD object')
    region_list[[i]] <- shorter_barcodes
    if(length(shorter_barcodes) < 2)
      stop('run.de.regions: error, region_list must be a list of atomic character vectors of length at least 2.')
  }
  barcodes <- Reduce(union, region_list)
  X2 <- matrix(0, nrow = length(barcodes), ncol = n_regions)
  rownames(X2) <- barcodes;
  for(i in 1:n_regions)
    X2[region_list[[i]],i] <- 1
  return(X2)
}

