
<!-- README.md is generated from README.Rmd. Please edit that file -->

# INSERT TITLE OF THE PAPER

<!-- badges: start -->

<!-- badges: end -->

Here, we will explain how the analysis occured for our paper ‘Cell
type-specific differential expression for spatial transcriptomics’,
which introduces and validates the Robust Cell Type
Differential-Expression (RCTDE) method for detecting cell type-specific
differential expression in spatial transcriptomics. You may access RCTDE
within our open-source R package
[here](https://github.com/dmcable/RCTD).

### DE Simulation

Using the helper functions in
[de\_simulation\_helper.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/helper_functions/de_simulation_helper.R),
code and results of RCTDE on simulated data is shown in Figure 2 below.

### Cerebellum Slide-seq

The script
[cer\_reps.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/Preprocessing_and_RCTD/cer_reps.R)
pre-processes the cerebellum Slide-seq data including cropping ROIs and
creating `SpatialRNA` objects. This script also runs RCTD on our three
cerebellum replicates.

The scripts
[run\_de\_cer\_08.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/run_RCTDE/run_de_cer_08.R),
[run\_de\_cer\_09.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/run_RCTDE/run_de_cer_09.R),
and
[run\_de\_cer\_11.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/run_RCTDE/run_de_cer_11.R)
run RCTDE on the three Slide-seq replicates. Merging multiple replicates
uses functions from
[merge\_de\_helper.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/helper_functions/merge_de_helper.R).
The script
[hcr.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/image_analysis/hcr.R)
quantitatively processes the HCR data. Results are shown in Figure 3
below.

### Testes Slide-seq

The script
[run\_RCTD\_testes.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/Preprocessing_and_RCTD/run_RCTD_testes.R)
loads in the testes data and runs RCTD. Next, RCTDE is run using the
script
[run\_de\_testes.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/run_RCTD/run_de_testes.R).
For downstream testes analysis (results in Figure 4 - Testes below),
helper functions are used from
[testes\_helper.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/helper_functions/testes_helper.R).

### J20 Hippocampus Slide-seq

The file
[preprocess\_j20.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/j20/preprocess_j20.R)
preprocesses the Slide-seq data including cropping ROI. Next, RCTD is
run using the scripts
[run\_RCTD\_j20\_21.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/j20/run_RCTD_j20_21.R),
[run\_RCTD\_j20\_22.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/j20/run_RCTD_j20_22.R),
[run\_RCTD\_j20\_23.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/j20/run_RCTD_j20_23.R),
and
[run\_RCTD\_j20\_24.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/j20/run_RCTD_j20_24.R).

In the following analysis, helper functions from
[alzheimers\_helper.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/helper_functions/alzheimers_helper.R)
are used.

In order to create the RCTDE predictive variable, we must next align the
amyloid plaque antibody stains to the Slide-seq data. To do this,
[align\_plaque\_density.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/j20/align_plaque_density.R)
aligns the images (manually) using the DAPI channel to the Slide-seq
data. Next, the script
[created\_pd.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/j20/create_pd.R)
creates and saves smoothed versions of the plaque images. Finally,
[pd\_to\_exvar.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/j20/pd_to_exvar.R)
calculates the plaque density predictive variable at each Slide-seq
pixel.

After generating the predictive variable, we run RCTDE on the four
samples using the scripts
[run\_RCTD\_j20\_21.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/j20/run_RCTD_j20_21.R),
[run\_RCTD\_j20\_22.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/j20/run_RCTD_j20_22.R),
[run\_RCTD\_j20\_23.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/j20/run_RCTD_j20_23.R),
and
[run\_RCTD\_j20\_24.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/j20/run_RCTD_j20_24.R).
Overall results are shown in Figure 4 - J20 below, which uses helper
functinos from
[merge\_de\_helper.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/helper_functions/merge_de_helper.R).

### MERFISH Hypothalamus Slide-seq

The script
[run\_de\_merfish.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/run_RCTD/run_de_merfish.R)
runs RCTDE (including RCTD) on the MERFISH hypothalamus data for both
linear and quadratic models. Results and other analysis are shown in
Figure 4 - MERFISH below.

### KP Tumor Slide-seq

The script
[overlay.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/image_analysis/overlay.R)
generates a plot of the KP tumor H\&E annotations.

The script
[run\_de\_tumor.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/run_RCTD/run_de_tumor.R)
runs RCTDE (including RCTD) on the Slide-seq KP tumor data for the
parametric case (DE as a function of immune cell density). Results and
other analysis are shown in Figure 5 - parametric below.

The script
[run\_de\_nonparametric.R](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/run_RCTD/run_de_nonparametric.R)
runs RCTDE nonparametrically on the tumor data. Results and other
analysis are shown in Figure 5 - nonparametric below.

### Running RCTDE

For each dataset, RCTD and RCTDE were run according to the instructions
for the [RCTD package](https://github.com/dmcable/RCTD).

### Generating Main Figures

We provide R Markdown files that were used to create the main figures:

  - [Validating RCTDE on simulated
    data](https://raw.githack.com/dmcable/RCTD/master/AnalysisRCTDE/Figures/figure2.html)
    (Figure 2)
  - [RCTDE on Slide-seq
    cerebellum](https://raw.githack.com/dmcable/RCTD/master/AnalysisRCTDE/Figures/figure3.html)
    (Figure 3)
  - [RCTDE on Slide-seq J20
    hippocampus](https://raw.githack.com/dmcable/RCTD/master/AnalysisRCTDE/Figures/figure4_j20.html)
    (Figure 4)
  - [RCTDE on Slide-seq
    testes](https://raw.githack.com/dmcable/RCTD/master/AnalysisRCTDE/Figures/figure4_testes.html)
    (Figure 4)
  - [RCTDE on MERFISH
    hypothalamus](https://raw.githack.com/dmcable/RCTD/master/AnalysisRCTDE/Figures/figure4_merfish.html)
    (Figure 4)
  - [Parametric RCTDE on Slide-seq KP
    tumor](https://raw.githack.com/dmcable/RCTD/master/AnalysisRCTDE/Figures/figure5_parametric.html)
    (Figure 5)
  - [Nonparametric RCTDE on Slide-seq KP
    tumor](https://raw.githack.com/dmcable/RCTD/master/AnalysisRCTDE/Figures/figure5_nonparametric.html)
    (Figure 5)

### Supplementary Figures

R Markdown code for generating supplemental figures can be found at:

  - [Supplementary figures
    part 1](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/Figures/supp1.Rmd)
  - [Supplementary figures
    part 2](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/Figures/supp2.Rmd)
  - [Supplementary figures
    part 3](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/Figures/supp3.Rmd)
  - [Supplementary figures
    part 4](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/Figures/supp4.Rmd)
  - [Supplementary figures
    part 5](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/Figures/supp5.Rmd)
  - [Supplementary figures
    part 6](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/Figures/supp6.Rmd)
  - [Supplementary figures
    part 7](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/Figures/supp7.Rmd)

### RCTDE Results

A list of significant genes found by RCTDE on all datasets (for each
cell type) can be found at [RCTDE
results](https://github.com/dmcable/RCTD/tree/master/AnalysisRCTDE/paper_results)
