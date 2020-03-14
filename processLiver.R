library(readr)
dge <- read_csv('Data/Reference/Liver/GSM3714747_chow1_filtered_gene_bc_matrices.csv')
dgeMat <- as(as(dge,"matrix"),"dgCMatrix")
nUMI = colSums(dge[,-1])
