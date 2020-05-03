#used to find genes within astrocytes dependent on environment

ast_gene_list <- rownames(iv$cell_type_info[[1]][marker_data$cell_type == "Astrocyte",])
ast_gene_list <- intersect(ast_gene_list, rownames(puck_d@counts))
ast_gene_list <- ast_gene_list[iv$cell_type_info[[1]][ast_gene_list,"Astrocyte"]/apply(iv$cell_type_info[[1]][ast_gene_list,cell_type_list_not],1,max) >= 10]
non_reject <- results_df$spot_class != "reject"
loc_df <- cbind(results_df[non_reject, "first_type"], puck@coords[non_reject, c('x','y')])
doub_ind <- results_df$spot_class == "doublet_certain"
loc_df2 <- cbind(results_df[doub_ind, "second_type"], puck@coords[doub_ind, c('x','y')])
colnames(loc_df) = c('cell_type','x','y'); colnames(loc_df2) = c('cell_type','x','y')
loc_df <- rbind(loc_df, loc_df2)
conversion <- .65
D_vec_microns = 40
D_vec = D_vec_microns / conversion
D = D_vec
library(fields)
library(tidyr)
mat <- rdist(loc_df[,c('x','y')])
mat_ind = mat <= D
xcoord <- floor((which(mat_ind)-1)/dim(loc_df)[1]) + 1
ycoord <- (which(mat_ind)-1) %% dim(loc_df)[1] + 1
small_df <- data.frame(xcoord, ycoord, mat[mat_ind])
colnames(small_df) = c('ind1', 'ind2', 'dist')
small_df$type1 <- loc_df[small_df$ind1,"cell_type"]
small_df$type2 <- loc_df[small_df$ind2,"cell_type"]
neighbors <- as.data.frame((small_df %>% group_by(ind1, type2) %>% summarise(count = n())) %>% pivot_wider(id_cols = ind1, names_from = type2, values_from=count))
neighbors[is.na(neighbors)] <- 0
my_map <- c(which(results_df[non_reject,]$spot_class == "doublet_certain"),tail((1:dim(loc_df)[1]),sum(results_df$spot_class == "doublet_certain")),which(results_df[non_reject,]$spot_class == "singlet"))
my_neigh <- neighbors[my_map,]
rownames(my_neigh) <- colnames(puck_d@counts)
my_neigh$ind1 <- NULL
hist(my_neigh[puck_d@cell_labels=="Astrocyte","Astrocyte"]/rowSums(my_neigh[puck_d@cell_labels=="Astrocyte",]))
ast_neigh <- my_neigh[puck_d@cell_labels=="Astrocyte",]
ast_prop <- sweep(ast_neigh, 1, rowSums(ast_neigh),'/')
cutoff <- 0.25
max_type <- apply(ast_prop,1, function(x) colnames(ast_prop)[which.max(x[2:17])+1])
for(cell_type in cell_type_names) {
  print(cell_type)
  print(sum(ast_prop[,cell_type] > 0.25 & max_type == cell_type))
}

plot_puck_continuous(puck_d,rownames(ast_prop)[ast_prop[,"Astrocyte"] > 0.8], 1,title="Astrocyte")
cell_type_list <- c("CA1","CA3","Astrocyte","Oligodendrocyte","Denate","Interneuron")
cell_type="Oligodendrocyte"
gene_mean_mat <- Matrix(0,nrow = length(ast_gene_list), ncol = length(cell_type_list))
gene_sd_mat <- Matrix(0,nrow = length(ast_gene_list), ncol = length(cell_type_list))
rownames(gene_mean_mat) = ast_gene_list; colnames(gene_mean_mat) = cell_type_list
rownames(gene_sd_mat) = ast_gene_list; colnames(gene_sd_mat) = cell_type_list
for(cell_type in cell_type_list) {
  if(cell_type == "Astrocyte")
    my_ind <- rownames(ast_prop)[ast_prop[,"Astrocyte"] > 0.8]
  else
    my_ind <- rownames(ast_prop)[ast_prop[,cell_type] > 0.25 & max_type == cell_type]
  gene_mean_mat[,cell_type] <- rowMeans(puck_d@counts[ast_gene_list,my_ind])
  gene_sd_mat[,cell_type] <- apply(puck_d@counts[ast_gene_list,my_ind],1,sd)/sqrt(length(my_ind))
  #plot_puck_continuous(puck_d,my_ind, puck_d@counts[gene,],title=cell_type,ylimit = c(0,1e-2))
  #plot_puck_continuous(puck_d,my_ind, puck_d@counts[gene,]*puck_d@nUMI,title=cell_type,ylimit = c(0,1))
}
Z = abs(gene_mean_mat - gene_mean_mat[,"Astrocyte"])/(sqrt(gene_sd_mat^2 + gene_sd_mat[,"Astrocyte"]^2))
log_fc <- log(apply(gene_mean_mat,1,max)/apply(gene_mean_mat,1,mean),2)
log_fc <- log_fc[order(-log_fc)]
save(log_fc, gene_mean_mat, ast_gene_list, cell_type_list, ast_prop, puck_d, max_type, file = 'Plotting/Results/Astrocytes.RData')
save(log_fc, gene_mean_mat, gene_sd_mat, ast_gene_list, cell_type_list, ast_prop, puck_d, max_type, res_df, pos_genes, neg_genes, file = 'Plotting/Results/AstrocytesFull.RData')
