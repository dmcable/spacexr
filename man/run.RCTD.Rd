% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spacexr.R
\name{run.RCTD}
\alias{run.RCTD}
\title{Runs the RCTD pipeline on a \code{\linkS4class{RCTD}} object}
\usage{
run.RCTD(RCTD, doublet_mode = "doublet")
}
\arguments{
\item{RCTD}{an \code{\linkS4class{RCTD}} object created using the \code{\link{create.RCTD}} function.}
}
\value{
an \code{\linkS4class{RCTD}} object containing the results of the RCTD algorithm. Please see \code{\linkS4class{RCTD}}
documentation for more information on interpreting the content of the RCTD object.
}
\description{
Equivalent to sequentially running the functions \code{\link{fitBulk}}, \code{\link{choose_sigma_c}}, and \code{\link{fitPixels}}
}
\details{
If in doublet mode, fits at most two cell types per pixel. It classifies each pixel as 'singlet' or 'doublet' and searches for the cell types
on the pixel. If in full mode, can fit any number of cell types on each pixel. In multi mode, cell types are added using a greedy algorithm,
up to a fixed number.
}
