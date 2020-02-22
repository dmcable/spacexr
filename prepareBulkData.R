#this script takes the single cell reference and puck and prepares the bulk data files so that you can run fitBulkWeights.py
#working directory must be the Cell Demixing folder.
library(RCTD)
library(Matrix)
print("prepareBulkData: begin")
iv <- init_RCTD(gene_list_reg = F, get_proportions = F) #initial variables
prepareBulkData(iv$bulkdir, iv$cell_type_info[[1]], iv$puck, iv$gene_list)
print("prepareBulkData: end")
