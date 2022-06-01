
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Cell type-specific inference of differential expression for spatial transcriptomics

<!-- badges: start -->
<!-- badges: end -->

Here, we will explain how the analysis occurred for our paper ‘Cell
type-specific inference of differential expression in spatial
transcriptomics’, which introduces and validates the Cell type-Specific
Inference of Differential Expression (C-SIDE) method for detecting cell
type-specific differential expression in spatial transcriptomics. You
may access C-SIDE within our open-source R package
[here](https://github.com/dmcable/spacexr).

### Obtaining Data

The data generated and/or used in this study may be accessed at the
[Broad Institute’s Single Cell
Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP1663).
This repository contains both the Slide-seq datasets used in this study,
and the single-cell RNA-sequencing references.

### DE Simulation

Using the helper functions in
[de_simulation_helper.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/helper_functions/de_simulation_helper.R),
code and results of C-SIDE on simulated data is shown in Figure 2 below.

### Cerebellum Slide-seq

The script
[cer_reps.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/Preprocessing_and_RCTD/cer_reps.R)
pre-processes the cerebellum Slide-seq data including cropping ROIs and
creating `SpatialRNA` objects. This script also runs RCTD on our three
cerebellum replicates.

The scripts
[run_de_cer_08.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/run_CSIDE/run_de_cer_08.R),
[run_de_cer_09.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/run_CSIDE/run_de_cer_09.R),
and
[run_de_cer_11.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/run_CSIDE/run_de_cer_11.R)
run C-SIDE on the three Slide-seq replicates. Merging multiple
replicates uses functions from
[merge_de_helper.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/helper_functions/merge_de_helper.R).
The script
[hcr.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/image_analysis/hcr.R)
quantitatively processes the HCR data. Results are shown in Figure 3
below.

### Testes Slide-seq

The script
[run_RCTD_testes.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/Preprocessing_and_RCTD/run_RCTD_testes.R)
loads in the testes data and runs RCTD. Next, C-SIDE is run using the
script
[run_de_testes.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/run_CSIDE/run_de_testes.R).
For downstream testes analysis (results in Figure 4 - Testes below),
helper functions are used from
[testes_helper.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/helper_functions/testes_helper.R).

### J20 Hippocampus Slide-seq

The file
[preprocess_j20.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/j20/preprocess_j20.R)
preprocesses the Slide-seq data including cropping ROI. Next, RCTD is
run using the scripts
[run_RCTD_j20_21.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/j20/run_RCTD_j20_21.R),
[run_RCTD_j20_22.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/j20/run_RCTD_j20_22.R),
[run_RCTD_j20_23.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/j20/run_RCTD_j20_23.R),
and
[run_RCTD_j20_24.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/j20/run_RCTD_j20_24.R).

In the following analysis, helper functions from
[alzheimers_helper.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/helper_functions/alzheimers_helper.R)
are used.

In order to create the C-SIDE predictive variable, we must next align
the amyloid plaque antibody stains to the Slide-seq data. To do this,
[align_plaque_density.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/j20/align_plaque_density.R)
aligns the images (manually) using the DAPI channel to the Slide-seq
data. Next, the script
[created_pd.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/j20/create_pd.R)
creates and saves smoothed versions of the plaque images. Finally,
[pd_to_exvar.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/j20/pd_to_exvar.R)
calculates the plaque density predictive variable at each Slide-seq
pixel.

After generating the predictive variable, we run C-SIDE on the four
samples using the scripts
[run_CSIDE_j20_21.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/j20/run_CSIDE_j20_21.R),
[run_CSIDE_j20_22.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/j20/run_CSIDE_j20_22.R),
[run_CSIDE_j20_23.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/j20/run_CSIDE_j20_23.R),
and
[run_CSIDE_j20_24.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/j20/run_CSIDE_j20_24.R).
Overall results are shown in Figure 4 - J20 below, which uses helper
functinos from
[merge_de_helper.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/helper_functions/merge_de_helper.R).

### MERFISH Hypothalamus Slide-seq

The script
[run_de_merfish.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/run_CSIDE/run_de_merfish.R)
runs C-SIDE (including RCTD) on the MERFISH hypothalamus data for both
linear and quadratic models. Results and other analysis are shown in
Figure 4 - MERFISH below.

### KP Tumor Slide-seq

The script
[overlay.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/image_analysis/overlay.R)
generates a plot of the KP tumor H&E annotations.

The script
[run_de_tumor.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/run_CSIDE/run_de_tumor.R)
runs C-SIDE (including RCTD) on the Slide-seq KP tumor data for the
parametric case (DE as a function of immune cell density). Results and
other analysis are shown in Figure 5 - parametric below.

The script
[run_de_nonparametric.R](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/run_CSIDE/run_de_nonparametric.R)
runs C-SIDE nonparametrically on the tumor data. Results and other
analysis are shown in Figure 5 - nonparametric below.

### Running C-SIDE

For each dataset, RCTD and C-SIDE were run according to the instructions
for the [spacexr package](https://github.com/dmcable/spacexr).

### Generating Main Figures

We provide R Markdown files that were used to create the main figures:

-   [Validating C-SIDE on simulated
    data](https://raw.githack.com/dmcable/spacexr/master/AnalysisCSIDE/Figures/figure2.html)
    (Figure 2)
-   [C-SIDE on Slide-seq
    cerebellum](https://raw.githack.com/dmcable/spacexr/master/AnalysisCSIDE/Figures/figure3.html)
    (Figure 3)
-   [C-SIDE on Slide-seq J20
    hippocampus](https://raw.githack.com/dmcable/spacexr/master/AnalysisCSIDE/Figures/figure4_j20.html)
    (Figure 4)
-   [C-SIDE on Slide-seq
    testes](https://raw.githack.com/dmcable/spacexr/master/AnalysisCSIDE/Figures/figure4_testes.html)
    (Figure 4)
-   [C-SIDE on MERFISH
    hypothalamus](https://raw.githack.com/dmcable/spacexr/master/AnalysisCSIDE/Figures/figure4_merfish.html)
    (Figure 4)
-   [C-SIDE on Visium lymph
    node](https://raw.githack.com/dmcable/spacexr/master/AnalysisCSIDE/Figures/figure4_visium.html)
    (Figure 5)
-   [Parametric C-SIDE on Slide-seq KP
    tumor](https://raw.githack.com/dmcable/spacexr/master/AnalysisCSIDE/Figures/figure5_parametric.html)
    (Figure 6)
-   [Nonparametric C-SIDE on Slide-seq KP
    tumor](https://raw.githack.com/dmcable/spacexr/master/AnalysisCSIDE/Figures/figure5_nonparametric.html)
    (Figure 6)

### Supplementary Figures

R Markdown code for generating supplemental figures can be found at:

-   [Supplementary figures part
    1](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/Figures/supp1.Rmd)
-   [Supplementary figures part
    2](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/Figures/supp2.Rmd)
-   [Supplementary figures part
    3](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/Figures/supp3.Rmd)
-   [Supplementary figures part
    4](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/Figures/supp4.Rmd)
-   [Supplementary figures part
    5](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/Figures/supp5.Rmd)
-   [Supplementary figures part
    6](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/Figures/supp6.Rmd)
-   [Supplementary figures part
    7](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/Figures/supp7.Rmd)
-   [Supplementary figures part
    8](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/Figures/supp8.Rmd)

### C-SIDE Results

A list of significant genes found by C-SIDE on all datasets (for each
cell type) can be found at [C-SIDE
results](https://github.com/dmcable/spacexr/tree/master/AnalysisCSIDE/paper_results)
