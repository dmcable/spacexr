library(RCTD)
library(tidyr)
iv <- init_RCTD(gene_list_reg = F, get_proportions = T)
proportions = iv$proportions
cell_type_info_unnorm <- iv$cell_type_info

puck = iv$puck
#iv <- init_RCTD(load_info_renorm = T) #initial variables


bulk_vec = rowSums(puck@counts);
total_UMI <- sum(puck@nUMI)
gene_means_unnorm_pred <- as.matrix(cell_type_info_unnorm[[1]][iv$gene_list,]) %*% proportions
pred_platform_effect = log(bulk_vec[iv$gene_list] / total_UMI,2) -log(gene_means_unnorm_pred,2)
hist(pred_platform_effect, breaks =50,xlab = "log2(Platform Effect)", main = "Measured Platform effects")
N_up = 250; N_down = 300;
write.csv(names(head(pred_platform_effect[order(pred_platform_effect),],N_down)),file.path(resultsdir,'down_genes.csv'))
write.csv(names(tail(pred_platform_effect[order(pred_platform_effect),],N_up)),file.path(resultsdir,'up_genes.csv'))
#do the gene ontology
up_go <- read.csv(file.path(resultsdir,'GO/up_all.csv'))
rownames(up_go) <- as.character(up_go$GO.cellular.component.complete)
down_go <- read.csv(file.path(resultsdir,'GO/down_all.csv'))
rownames(down_go) <- as.character(down_go$GO.cellular.component.complete)
down_go <- down_go[rownames(up_go),]
gene_list <- rownames(down_go)[(up_go$upload_1..FDR. < 0.05 | down_go$upload_1..FDR. < 0.05)]
up_scores <- unlist(lapply(as.character(up_go[gene_list,]$upload_1..fold.Enrichment.), parse_number))
down_scores <- unlist(lapply(as.character(down_go[gene_list,]$upload_1..fold.Enrichment.), parse_number))
log_fc <- log(up_scores/down_scores)
min_genes = 5; FC_thresh = 0.25
frequent <- rownames(up_go)[up_go$upload_1..238. > min_genes | down_go$upload_1..292. > min_genes]
write.csv(log_fc[log_fc > FC_thresh & names(log_fc) %in% frequent],file.path(resultsdir,'GO/up_results.csv'))
write.csv(log_fc[log_fc < -FC_thresh & names(log_fc) %in% frequent],file.path(resultsdir,'GO/down_results.csv'))
if(F) {
  names(up_scores) <- gene_list; names(down_scores) <- gene_list
  common <- intersect(rownames(up_go), rownames(down_go))
  up_genes <- setdiff(rownames(up_go), rownames(down_go))
  up_genes <- up_genes[up_go[up_genes,]$upload_1..fold.Enrichment. >= 1.5]
  down_genes <- setdiff(rownames(down_go), rownames(up_go))
  down_genes <- down_genes[down_go[down_genes,]$upload_1..fold.Enrichment. >= 1.5]
  log_fc <- log(up_go[common, ]$upload_1..fold.Enrichment. / down_go[common, ]$upload_1..fold.Enrichment.)
  names(log_fc) <- common
  down_genes = c(down_genes,common[log_fc < -FC_thresh])
  up_genes = c(up_genes,common[log_fc > FC_thresh])
}
