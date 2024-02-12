
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spacexr (Spatial-eXpression-R): Robust Cell Type Decomposition (RCTD) and Cell type-Specific Inference of Differential Expression (C-SIDE)

<!-- badges: start -->
<!-- badges: end -->

Welcome to *spacexr*, an R package for learning cell types and cell
type-specific differential expression in spatial transcriptomics data.

Robust Cell Type Decomposition (RCTD) inputs a spatial transcriptomics
dataset, which consists of a set of *pixels*, which are spatial
locations that measure RNA counts across many genes. RCTD additionally
uses a single cell RNA-seq (scRNA-seq) dataset, which is labeled for
cell types. RCTD learns cell type profiles from the scRNA-seq dataset,
and uses these to label the spatial transcriptomics pixels as cell
types. RCTD has been tested across a variety of spatial transcriptomics
technologies including imaging-based (e.g. MERFISH) and sequencing-based
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

Cell type-Specific Inference of Differential Expression (C-SIDE) learns
cell type-specific differential expression on spatial transcriptomics
dataset. C-SIDE inputs one or more user-defined *covariates*, which are
biologically-relevant axes along which differential expression is
hypothesized. These variables can be generally defined, but some choices
for covariates include spatial position, cellular microenvironment, or
cell-to-cell interactions. C-SIDE then identifies, for each cell type,
genes that are significantly differentially expressed as a function of
the covariates. Similar to RCTD, C-SIDE can operate on spatial data
containing single cells or cell type mixtures. Differential expression
can be learned at the population level across one or multiple biological
replicates or samples.

Code for generating the figures of our RCTD paper, Robust decomposition
of cell type mixtures in spatial transcriptomics, is located
[here](https://github.com/dmcable/spacexr/tree/master/AnalysisPaper).
Our *Nature Biotechnology* paper can be found
[here](https://www.nature.com/articles/s41587-021-00830-w).

Code for generating the figures of our C-SIDE paper, Cell type-specific
differential expression inference in spatial transcriptomics, is located
[here](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE).
Our C-SIDE paper is available in *Nature Methods*
[here](https://www.nature.com/articles/s41592-022-01575-3).

## News and Updates

February 1st, 2023: Version 2.2.0 released. This is a significant
expansion of the functionality of C-SIDE. For multi-sample analysis,
using meta-regression, we have added the ability to add in covariates,
such as time, age, case vs. control, etc, and to test along these
covariates. An example can be found in [Population-level RCTD and
C-SIDE](https://raw.githack.com/dmcable/spacexr/master/vignettes/replicates.html)
Vignette under the “Population inference: meta regression” subheader.

Additionally, linear interpolation for the log-likelihood has been
replaced with cubic spline interpolation. The result is increased
numerical stability, convergence, and speed, especially for C-SIDE.

December 13th, 2022: Version 2.1.0 released. This is a substantial
update to the behind-the-scenes implementation. Linear interpolation for
the log-likelihood has been replaced with cubic spline interpolation.
The result is increased numerical stability, convergence, and speed,
especially for C-SIDE.

September 1st, 2022: Our C-SIDE paper has been published in *Nature
Methods* [here](https://www.nature.com/articles/s41592-022-01575-3). We
have also written a research briefing about this article in *Nature
Methods* [here](https://www.nature.com/articles/s41592-022-01576-2).

December 22nd, 2021: We are renaming this package (formerly RCTD) as
*spacexr* (Spatial eXpression R package). We are also releasing
*spacexr* 2.0, now featuring cell type-specific differential expression.
The new algorithm, called C-SIDE, is introduced in our new paper which
will be available soon. <!--[here](BIORXIV LINK). --> We are also
introducing a feature where RCTD and C-SIDE can be run in batch across
multiple experimental replicates.

March 18th, 2021: Our RCTD paper has been published in *Nature
Biotechnology*
[here](https://www.nature.com/articles/s41587-021-00830-w). Also, we
have just released a new version of RCTD (version 1.2.0) with a
different method to input data (please see `SpatialRNA` and `Reference`
constructor functions). Please report any bugs associated with this new
update.

## Installation

You can install the current version of *spacexr* from
[GitHub](https://github.com/dmcable/spacexr) with:

``` r
# install.packages("devtools")
options(timeout = 600000000) ### set this to avoid timeout error
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
```

If you would like to build vignettes (it will take some time), modify
the above by setting `build_vignettes = TRUE`.

## Vignettes and Documentation

A complete guide to spacexr vignettes can be found
[here](https://github.com/dmcable/spacexr/tree/master/vignettes), which
includes diverse applications of RCTD and C-SIDE on Slide-seq, MERFISH,
and Visium datasets. These vignettes also include multiple applications
of C-SIDE to differential expression problems including spatial
position, cell-to-cell interactions, and joint analysis of multiple
experimental samples/replicates.

The *spacexr* manual can be found
[here](https://github.com/dmcable/spacexr/blob/master/spacexr_manual_2.2.1.pdf).

Additional detailed recommended reading (documentation, tutorials, and
tips) can be found
[here](https://github.com/dmcable/spacexr/tree/master/documentation).

## Quick Guide to Getting Started with RCTD

In this section, we aim to explain how to use RCTD as quickly as
possible on your data:

1.  Open the ‘spatial-transcriptomics.Rmd’ vignette for a complete
    explanation of the RCTD workflow. Expected output of the vignette is
    provided
    [here](https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html).
    If you have any questions about how to run RCTD, please first make
    sure you run this Vignette on your computer and make sure you
    understand how it works. For other modes of RCTD, be sure to check
    out the rest of our
    [vignettes](https://github.com/dmcable/spacexr/tree/master/vignettes).

2.  As described in the ‘Data Preprocessing’ step of the vignette,
    convert your spatial transcriptomics data to a `SpatialRNA` object
    (called `puck` here) and your scRNA-seq reference to a `Reference`
    object (called `reference` here). In order to create these objects,
    you need to load your coordinate matrices, counts matrices, etc
    into R. Type `?Reference` or `?SpatialRNA` into R to learn more
    about these constructor functions.

3.  Run RCTD. You can optionally set `test_mode` to `TRUE` in
    `create.RCTD` to quickly test RCTD, but you should set it to `FALSE`
    for the official run.

``` r
myRCTD <- create.RCTD(puck, reference, max_cores = 8, test_mode = FALSE) # here puck is the SpatialRNA object, and reference is the Reference object.
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
```

4.  Observe RCTD results. RCTD results are stored in the `@results`
    field. Of particular interest is `@results$weights`, a data frame of
    cell type weights for each pixel (for full mode \[i.e. doublet_mode
    = ‘full’\]). Alternatively This section will generate various plots
    which can be found in `resultsdir`. The results of
    ‘doublet_mode=“doublet”’ are stored in `@results$results_df` and
    `@results$weights_doublet`, the weights of each cell type. More
    specifically, the `results_df` object contains one column per pixel
    (barcodes as rownames). Important columns are:

-   `spot_class`, a factor variable representing RCTD’s classification
    in doublet mode: “singlet” (1 cell type on pixel), “doublet_certain”
    (2 cell types on pixel), “doublet_uncertain” (2 cell types on pixel,
    but only confident of 1), “reject” (no prediction given for pixel).
-   Next, the `first_type` column gives the first cell type predicted on
    the bead (for all spot_class conditions except “reject”).
-   The `second_type column` gives the second cell type predicted on the
    bead for doublet spot_class conditions (not a confident prediction
    for “doublet_uncertain”).

For some example of summary plots, follow the ‘RCTD results’ section of
the ‘spatial-transcriptomics’ vignette.

## Quick Guide to Getting Started with C-SIDE

Here, we discuss how to run C-SIDE on your data and detect cell
type-specific differential expression:

1.  Open the ‘differential-expression.Rmd’ vignette for a complete
    explanation of the C-SIDE workflow. Expected output of the vignette
    is provided
    [here](https://raw.githack.com/dmcable/spacexr/master/vignettes/differential-expression.html).
    If you have any questions about how to run C-SIDE, please first make
    sure you run this Vignette on your computer and make sure you
    understand how it works. For other applications of C-SIDE, be sure
    to check out the rest of our
    [vignettes](https://github.com/dmcable/spacexr/tree/master/vignettes).

2.  Assign cell types to your spatial transcriptomics dataset. Since
    C-SIDE detects cell type-specific differential expression, it needs
    to first identify cell types. It is recommended to use the internal
    RCTD procedure to identify cell types, as described above. However,
    cell types can also be imported from another source using the
    `import_weights` function.

3.  As shown in the vignettes, the most important step is to define the
    covariates along which to detect differential expression. Each
    covariate should be a numeric vector, representing the value of the
    covariate at each pixel. Please standardize each covariate between 0
    and 1, and the names of each covariate variable should match the
    names, or barcodes, of the spatial transcriptomics pixels.

4.  Run C-SIDE. If a single covariate is present, the simplest way to
    run C-SIDE is to use the `run.CSIDE.single` function as in the
    ‘differential-expression.Rmd’ vignette. Here you will pass in your
    covariate and the RCTD object containing previously-computed cell
    type assignments. If multiple covariates are present, consider our
    other functions for running: `run.CSIDE.regions` (multiple regions),
    `run.CSIDE.nonparametric` (nonparametric modelling), or `run.CSIDE`
    (inputs a general design matrix).

5.  Observe C-SIDE results. We recommend making plots representing DE
    results using the `make_all_de_plots` function. Also, DE results can
    be directly examined in the `@de_results` object. Please see the
    vignettes for more examples.

6.  Multiple replicates. If multiple experimental replicates are
    available, we recommend running RCTD and C-SIDE in batch mode
    (operates on a `RCTD.replicates` object) and using C-SIDE to do
    population-level statistical inference. This procedure is detailed
    in the [Population-level RCTD and C-SIDE
    vignette](https://raw.githack.com/dmcable/spacexr/master/vignettes/replicates.html).

### Dependencies

-   R version \>= 3.5.0.
-   R packages: readr, pals, ggplot2, Matrix, parallel, doParallel,
    foreach, quadprog, tibble, dplyr, reshape2, knitr, rmarkdown,
    fields, and mgcv.

For optimal performance, we recommend at least 4 GB of RAM, and multiple
cores may be used for RCTD and C-SIDE to speed up runtime.

Installation time: Less than five minutes, after installing dependent
packages. RCTD takes up approximately 145 MB of space due to
pre-computed data tables that substantially improve performance.

Runtime: The example datasets provided can be run in less than 5 minutes
on a normal desktop computer (both RCTD and C-SIDE vignettes). RCTD runs
in approximately 20 minutes on a Slide-seq cerebellum dataset
(approximately 3,000 genes and 11,000 pixels) on a laptop computer with
4 cores. Using less cores will lead to longer runtime. C-SIDE (excluding
cell type identification step) was found to run in approximately 15
minutes (4 cores) on differential expression between two regions for the
Slide-seq cerebellum dataset (approximately 2,776 pixels, 4,812 genes,
and 5 cell types used).

Operating systems (version 2.0 spacexr) tested on:

-   macOS Big Sur 11.6
-   GNU/Linux (GNU coreutils) 8.22 (version 1.0 spacexr tested)

### Citations

If you use this work for cell type estimation, please cite:

Cable, Dylan M., et al. “Robust decomposition of cell type mixtures in
spatial transcriptomics.” *Nature Biotechnology* 40.4 (2022): 517-526.

If you use this work for differential expression, please cite:

Cable, Dylan M., et al. “Cell type-specific inference of differential
expression in spatial transcriptomics.” *Nature methods* 19.9 (2022):
1076-1087.

### License

spacexr is
[licensed](https://github.com/dmcable/spacexr/blob/master/LICENSE) under
the GNU General Public License v3.0.
