
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
label the spatial transcriptomics pixels as cell types. Notably, RCTD
allows for individual pixels to be cell type *mixtures*; that is, they
can potentially source RNA from multiple cell types. RCTD identifies the
cell types on each pixel, and estimates the proportion of each of these
cell types. Additionally, RCTD has a platform effect normalization step,
which normalizes the scRNA-seq cell type profiles to match the *platform
effects* of the spatial transcriptomics dataset. A platform effect is
the tendency of a sequencing technology to capture individual genes at
different rates.

## Installation

You can install the current version of RCTD from
[GitHub](https://github.com/dmcable/RCTD) with:

``` r
# install.packages("devtools")
devtools::install_github("dmcable/RCTD")
```

### Downloading RCTD Folder

We reccommend downloading the ‘RCTD\_base’ folder, and adapting it to be
your working directory for RCTD. ‘RCTD\_base’ can be downloaded as
follows:

``` bash
#!/bin/bash
wget https://raw.githubusercontent.com/dmcable/RCTD/dev/RCTD_base.zip 
unzip RCTD_base.zip
```

For the rest of this document, we assume that you are working in
‘RCTD\_base’ as a base directory. This folder contains several
important subfolders/files:

  - conf: Configuration files required for running RCTD.
      - This folder must be placed in the directory where you are
        running RCTD.
      - Data locations are entered into the ‘conf/dataset.yml’ file, and
        RCTD parameters in ‘conf/default.yml’.
      - For selecting a different RCTD parameter configuration file, you
        can change the `config_mode` field in ‘conf/dataset.yml’ to
        point to e.g. ‘conf/default.yml’, ‘conf/test.yml’, or another
        file.
      - Set `n_puck_folds` for number of folds to split the dataset
        (RCTD runs on batches/folds of the dataset).
  - data: A recommended folder to contain the data. Contains two
    subfolders ‘Reference’ for containing the scRNA-seq data, and
    ‘SpatialRNA’ for containing the spatial transcriptomics data.
    Examples are included with the ‘Vignette’ subfolders. RCTD will
    place its results in the ‘SpatialRNA’ subdirectory.
  - R\_scripts: R scripts for running RCTD
      - See the following ‘Workflow’ and ‘Recommended Guidelines’
        sections for explanation of how to use these scripts.
  - bash\_scripts: bash scripts for submitting RCTD jobs to the computer
    cluster.
      - These sample scripts were configured for the Broad Institute’s
        PBS cluster, but they are not likely to directly work without
        modification on other clusters.
      - See the section ‘Running on a Cluster’ for more details.
  - ‘spatial-transcriptomics.Rmd’ a vignette for RCTD, which you can run
    after downloading the ‘RCTD\_base’ folder
      - This vignette explains most of the features of RCTD.

## Workflow

The basic workflow for RCTD consists of the following steps:

1.  Data Preprocessing. Representing the spatial transcriptomics data as
    a `SpatialRNA` object, and the scRNA-seq reference as a `Seurat`
    object.
      - Save these objects as ‘RDS’ files, as shown in the ‘Data
        Preprocessing’ section of the ‘spatial-transcriptomics’
        vignette.
2.  Setup. Edit the configuration files (e.g. ‘conf/dataset.yml’ and
    ‘conf/default.yml’) to point to the data files and e.g. determine
    parameters for selecting differentially expressed genes.
      - Follow the ‘Setup’ section of the ‘spatial-transcriptomics’
        vignette.
3.  Platform Effect Normalization. Estimates platform effects between
    scRNA-seq dataset and spatial transcriptomics dataset. Uses this to
    normalize cell type profiles.
      - Run the ‘R\_scripts/fitBulk.R’ script or, equivalently, follow
        the ‘Platform Effect Normalization’ section of the
        ‘spatial-transcriptomics’ vignette.
4.  Hyperparameter optimization (choosing sigma). Handles overdispersion
    by determining the maximum likelihood variance for RCTD’s lognormal
    random effects. Precomputes likelihood function.
      - Run the ‘R\_scripts/chooseSigma.R’ script or, equivalently,
        follow the ‘Hyperparameter optimization (choosing sigma)’
        section of the ‘spatial-transcriptomics’ vignette.
5.  Robust Cell Type Decomposition. Assigns cell types to each spatial
    transcriptomics pixel, and estimates the cell type proportions.
      - Run the ‘R\_scripts/callDoublets.R’ for each data fold. An
        example of how this works is outlined in the ‘Robust Cell Type
        Decomposition’ section of the ‘spatial-transcriptomics’
        vignette.
6.  Collecting RCTD results. Obtain results objects and make summary
    plots.
      - Run the ‘R\_scripts/gather\_results.R’ script or, equivalently,
        follow the ‘Collecting RCTD results’ section of the
        ‘spatial-transcriptomics’ vignette.

### Recommended Guidelines

We recommend running steps 1-2 in the R console (e.g. in RStudio), as
demonstrated in the ‘spatial-transcriptomics’ vignette. After setting
the configuration files, we recommend that steps 3-5 are run from
command line as follows (where there are, for example, three data
folds):

``` bash
#!/bin/bash
Rscript R_scripts/fitBulk.R
Rscript R_scripts/chooseSigma.R
Rscript R_scripts/callDoublets.R 1 # first fold of data
Rscript R_scripts/callDoublets.R 2 # second fold of data
Rscript R_scripts/callDoublets.R 3 # third fold of data
Rscript R_scripts/gather_results.R
```

Step 6 (gathering RCTD results) can be optionally run from command line
as above or in the R console.

### Running on a Cluster

One (recommended) option is to automate this workflow on a computer
cluster for steps 3-5. Each of steps 3-5 must be run sequentially, but
each fold may be run in parallel (as a job array) for step 5. We provide
an example of scripts used to automate RCTD on the Broad Institute’s PBS
cluster in the ‘bash\_scripts’ folder. Note that other clusters may have
different syntax for submitting jobs. If ‘pipeline\_sample.sh’ is edited
to have the correct number of folds, and ‘pipeline\_sample.sh’,
‘doub\_sample.sh’, and ‘pre\_sample.sh’ are edited with paths to the
RCTD folder, then one can, for example, run RCTD by the following
command:

``` bash
#!/bin/bash
bash ./bash_scripts/pipeline_sample.sh
```
