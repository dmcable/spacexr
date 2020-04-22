library(RCTD)
library(Matrix)
iv <- init_RCTD(gene_list_reg = F, get_proportions = F, load_info=F) #initial variables
bulkResults <- fitBulk(iv)

