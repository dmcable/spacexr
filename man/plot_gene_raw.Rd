% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CSIDE_plots.R
\name{plot_gene_raw}
\alias{plot_gene_raw}
\title{Makes a spatial plot of continuous gene expression for a particular gene}
\usage{
plot_gene_raw(myRCTD, gene, cell_type, ymax = 10)
}
\arguments{
\item{myRCTD}{\code{\linkS4class{RCTD}} object containing \code{de_results}, after running CSIDE}

\item{gene}{gene to be plotted}

\item{cell_type}{cell_type to be plotted (only single cell type pixels)}

\item{ymax}{(default 10) maximum expression (in counts per 500) for color scale}
}
\value{
gene expression plot
}
\description{
Units counts per 500
}
