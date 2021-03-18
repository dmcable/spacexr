
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Robust Cell Type Decomposition (RCTD)

<!-- badges: start -->

<!-- badges: end -->

Welcome to RCTD, an R package for assigning cell types to spatial
transcriptomics data. RCTD inputs a spatial transcriptomics dataset,
which consists of a set of *pixels*, which are spatial locations that
measure RNA counts across many genes. RCTD additionally uses a single
cell RNA-seq (scRNA-seq) dataset, which is labeled for cell types. RCTD
learns cell type profiles from the scRNA-seq dataset, and uses these to
label the spatial transcriptomics pixels as cell types. RCTD has been
tested across a variety of spatial transcriptomics technologies
including imaging-based (e.g. MERFISH) and sequencing-based
(e.g. Slide-seq, Visium). Notably, RCTD allows for individual pixels to
be cell type *mixtures*; that is, they can potentially source RNA from
multiple cell types. That said, RCTD can still handle the case where
there is only one cell per pixel. RCTD identifies the cell types on each
pixel, and estimates the proportion of each of these cell types.
Additionally, RCTD has a platform effect normalization step, which
normalizes the scRNA-seq cell type profiles to match the *platform
effects* of the spatial transcriptomics dataset. A platform effect is
the tendency of a sequencing technology to capture individual genes at
different rates than another sequencing technology.

Code for generating the figures of our paper, Robust decomposition of
cell type mixtures in spatial transcriptomics, is located
[here](https://github.com/dmcable/RCTD/tree/dev/AnalysisPaper). Our
*Nature Biotechnology* paper can be found
[here](https://www.nature.com/articles/s41587-021-00830-w).

## News and Updates

March, 18th, 2021: Our RCTD paper has been published in *Nature
Biotechnology*
[here](https://www.nature.com/articles/s41587-021-00830-w). Also, we
have just released a new version of RCTD (version 1.2.0) with a
different method to input data (please see `SpatialRNA` and `Reference`
constructor functions). Please report any bugs associated with this new
update.

## Installation

You can install the current version of RCTD from
[GitHub](https://github.com/dmcable/RCTD) with:

``` r
# install.packages("devtools")
devtools::install_github("dmcable/RCTD", build_vignettes = TRUE)
```

To subsequently view the vignette (recommended for learning how to use
RCTD):

``` r
browseVignettes('RCTD')
```

## Quick Guide to Getting Started with RCTD

In this section, we aim to explain how to use RCTD as quickly as
possible on your data:

1.  Open the ‘spatial-transcriptomics.Rmd’ vignette for a complete
    explanation of the RCTD workflow. Expected output of the vignette is
    provided
    [here](https://raw.githack.com/dmcable/RCTD/dev/vignettes/spatial-transcriptomics.html).
    If you have any questions about how to run RCTD, please first make
    sure you run this Vignette on your computer and make sure you
    understand how it works.

2.  As described in the ‘Data Preprocessing’ step of the vignette,
    convert your spatial transcriptomics data to a `SpatialRNA` object
    (called `puck` here) and your scRNA-seq reference to a `Reference`
    object (called `reference` here). In order to create these objects,
    you need to load your coordinate matrices, counts matrices, etc into
    R. Type `?Reference` or `?SpatialRNA` into R to learn more about
    these constructor functions.

3.  Run RCTD. You can optionally set `test_mode` to `TRUE` in
    `create.RCTD` to quickly test RCTD, but you should set it to `FALSE`
    for the official run.

<!-- end list -->

``` r
myRCTD <- create.RCTD(puck, reference, max_cores = 8, test_mode = FALSE) # here puck is the SpatialRNA object, and reference is the Reference object.
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
```

4.  Observe RCTD results. RCTD results are stored in the `@results`
    field. Of particular interest is `@results$weights`, a data frame of
    cell type weights for each pixel (for full mode \[i.e. doublet\_mode
    = ‘full’\]). Alternatively This section will generate various plots
    which can be found in `resultsdir`. The results of
    ‘doublet\_mode=“doublet”’ are stored in `@results$results_df`
    and `@results$weights_doublet`, the weights of each cell type. More
    specifically, the `results_df` object contains one column per pixel
    (barcodes as rownames). Important columns are:

<!-- end list -->

  - `spot_class`, a factor variable representing RCTD’s classification
    in doublet mode: “singlet” (1 cell type on pixel),
    “doublet\_certain” (2 cell types on pixel), “doublet\_uncertain”
    (2 cell types on pixel, but only confident of 1), “reject” (no
    prediction given for pixel).
  - Next, the `first_type` column gives the first cell type predicted on
    the bead (for all spot\_class conditions except “reject”).
  - The `second_type column` gives the second cell type predicted on the
    bead for doublet spot\_class conditions (not a confident prediction
    for “doublet\_uncertain”).

For some example of summary plots, follow the ‘RCTD results’ section of
the ‘spatial-transcriptomics’ vignette.

### Dependencies

  - R version \>= 3.5.0.
  - R packages: readr, pals, ggplot2, Matrix, doParallel, foreach,
    quadprog, tibble, dplyr, reshape2.

For optimal performance, we recommend at least 4 GB of RAM, and multiple
cores may be used to speed up runtime.

Installation time: Less than five minutes, after installing dependent
packages. RCTD takes up approximately 145 MB of space due to
pre-computed data tables that significantly improve performance.

Runtime: The example dataset provided (Vignette) can be run in less than
2 minutes on a normal desktop computer. Approximately 20 minutes on a
Slide-seq cerebellum dataset (approximately 3,000 genes and 11,000
pixels) on a laptop computer with 8 cores. Using less cores will lead to
longer runtime.

Operating systems (version 1.2.0 RCTD) tested on:

  - macOS Mojave 10.14.6
  - GNU/Linux (GNU coreutils) 8.22

### License

RCTD is [licensed](https://github.com/dmcable/RCTD/blob/master/LICENSE)
under the GNU General Public License v3.0.
