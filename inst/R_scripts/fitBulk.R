### Platform Effect Normalization

# The first step in running RCTD is the Platform Effect Normalization step.
# This is accomplished simply by feeding `iv` into the `fitBulk` function.
# Platform Effect Normalization simultaneously estimates bulk cell type proportion in the
# SpatialRNA dataset, and the platform effects (changes in gene expression from scRNA-seq to
# spatial transcriptomics). It uses the platform effects to renormalize the scRNA-seq
# cell type profiles, which is returned as `bulkResults$cell_type_info_renorm`.
# You can examine the predicted bulk cell type proportion (`bulkResults$proportions`).
# If this drastically does not fit your expectations of the cell type proportions,
# something may have gone wrong.

library(RCTD)
library(Matrix)
iv <- init_RCTD(gene_list_reg = F, get_proportions = F, load_info=F) #initial variables
bulkResults <- fitBulk(iv)

