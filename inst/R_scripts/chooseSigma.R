### Hyperparameter optimization: choosing sigma

# After running platform effect estimation, we can determine `sigma_c`,
# a hyperparameter that represents the variance of the random effects
# (note, this is also called $\sigma_\varepsilon$ in the paper). Note, that we reload in
# `init_RCTD` to include the normalized `cell_type_info` after platform effect normalization
# (it is now loading from 'MetaData/cell_type_info_renorm.RDS' in the SpatialRNA directory).
# The trace files for optimizing sigma can be found in the `resultsdir` directory.
# In fact, this results directory contains all the main results of RCTD.
# After selecting sigma, we calculate the quadrature, `Q_mat`, which is saved
# in `resultsdir` as `Q_mat.RDS`.

library(RCTD)
library(Matrix)
#puck_file = paste0("SplitPuck/puck",fold_index,".RDS"),
iv <- init_RCTD(MIN_OBS=0, load_info_renorm = T) #initial variables
resultsdir <- paste0(iv$SpatialRNAdir,"/results")
sigma <- choose_sigma_c(iv,resultsdir)
Q_mat <- get_Q_mat(iv, sigma)
saveRDS(Q_mat, file.path(resultsdir,'Q_mat.RDS'))
