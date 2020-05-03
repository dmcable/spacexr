library(RCTD)
library(readr)
slideseqdir = "Data/Slideseq/NewCerPuck_190926_08"
puck = readRDS(file.path(slideseqdir, "puckCropped.RDS"))
outdir = "NMFreg/DataCer"
nmf_weights <- as.data.frame(read_csv(file.path(outdir,"bead_deconv_df_norm.csv")))
rownames(nmf_weights) <- colnames(puck@counts)
nmf_weights <- nmf_weights[puck@nUMI >= 100,]
iv <- init_RCTD(load_info_renorm = T)
nmf_norm_weights <- nmf_weights[,1:iv$cell_type_info[[3]]]
resultsdir = "Data/Slideseq/NewCerPuck_190926_08/NMFResults/"
dir.create(resultsdir)
puck = iv$puck
thresh_nmf <- as.data.frame(read_csv(file.path(outdir,"thresh_certainty.csv")))
rownames(thresh_nmf) <- iv$cell_type_info[[2]]
plot_weights_nmf(iv$cell_type_info, puck, resultsdir, nmf_norm_weights, thresh_nmf)
plot_weights_unthreshold(iv$cell_type_info, puck, resultsdir, nmf_norm_weights)
plot_cond_occur_nmf(cell_type_info, resultsdir, nmf_norm_weights, thresh_nmf)
plot_occur_unthreshold(cell_type_info, resultsdir, nmf_norm_weights)
save(iv, nmf_norm_weights, thresh_nmf, file = "Plotting/Results/NMFreg")
