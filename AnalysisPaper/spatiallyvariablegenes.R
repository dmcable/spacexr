#find globally spatially variable genes

#spatial genes
hippo <- CreateSeuratObject(puck@counts, project = "Slideseq",assay = "RNA",min.cells = 0,min.features = 0,names.field = 1,names.delim = "_",meta.data = NULL)
hippo <- NormalizeData(hippo, verbose = FALSE)
hippo <- FindVariableFeatures(hippo, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
hippo <- ScaleData(hippo)
hippo[['image']] <- new(Class='SlideSeq',assay='RNA',coordinates = puck@coords[colnames(puck@counts),])
resultsdir = "Data/Slideseq/NewCerPuck_190926_08/SeuratResults/"
DefaultAssay(hippo) <- "RNA"
hippo <- FindSpatiallyVariableFeatures(hippo, assay = "RNA", slot = "scale.data", features = VariableFeatures(hippo)[1:1000],
                                       selection.method = "moransi", x.cuts = 100, y.cuts = 100)
spatial_genes <- t(iv$cell_type_info[[1]][head(SpatiallyVariableFeatures(hippo, selection.method = "moransi"),20),])
SpatialFeaturePlot(hippo, features = head(SpatiallyVariableFeatures(hippo, selection.method = "moransi"),
                                          20), ncol = 10, alpha = c(0.1, 1), max.cutoff = "q95")
library(DescTools)
#start here for plots
spatial_gini <- apply(spatial_genes,2,Gini)
all_gini <- apply(iv$cell_type_info[[1]][gene_list,],1,Gini)
hist(all_gini)
hist((rank(all_gini) / length(all_gini))[colnames(spatial_genes)], xlim=c(0,1))
means_df[colnames(spatial_genes),]$log_fc
fc_df <- data.frame(means_df[int_genes_CA3,],"Local")
CA3_spatial <- colnames(spatial_genes)[!is.na(means_df[colnames(spatial_genes),"log_fc"])]
fc_df2 <- data.frame(means_df[CA3_spatial,],"Global")
colnames(fc_df)[length(colnames(fc_df))] <- 'class'; colnames(fc_df2)[length(colnames(fc_df))] <- 'class'
fc_df <- rbind(fc_df,fc_df2)
p <- ggplot(fc_df, aes(x = class, y = log_fc)) +
  geom_jitter(alpha = 0.25) +
  geom_boxplot(fill = NA, outlier.alpha = 0) +
  theme_classic()
p


#find spatially variable genes within CA3
#filter for high expression
gene_list = intersect(rownames(cell_type_info[[1]]), rownames(puck_d@counts))
CA1_genes <- names(which(iv$cell_type_info[[1]][gene_list,]$CA1*2 > apply(iv$cell_type_info[[1]][gene_list,!(iv$cell_type_info[[2]] %in% c('CA1','CA3','Denate','Neuron.Slc17a6', "Entorihinal", "Neurogenesis"))],1,max)))
CA3_genes <- names(which(iv$cell_type_info[[1]][gene_list,]$CA3*2 > apply(iv$cell_type_info[[1]][gene_list,!(iv$cell_type_info[[2]] %in% c('CA1','CA3','Denate','Neuron.Slc17a6', "Entorihinal", "Neurogenesis"))],1,max)))
CA3_genes <- CA3_genes[rowMeans(sweep(puck_d@counts[CA3_genes,],2,puck_d@nUMI,'/')) > 2.5e-5]
Dentate_genes <- names(which(iv$cell_type_info[[1]][gene_list,]$Denate*2 > apply(iv$cell_type_info[[1]][gene_list,!(iv$cell_type_info[[2]] %in% c('CA1','CA3','Denate','Neuron.Slc17a6', "Entorihinal"))],1,max)))
get_int_genes <- function(gene_df,cell_type, gene_list_ct, puck_d, cell_barc) {
  pvals <- numeric(length(gene_list_ct)); names(pvals) = gene_list_ct
  cvs <- numeric(length(gene_list_ct)); names(cvs) = gene_list_ct
  toti = 0
  for(gene in gene_list_ct) {
    print(toti)
    toti = toti + 1
    gene_df$nUMI <- puck_d@nUMI[cell_barc]
    Tr = 100
    f_dist <- numeric(Tr)
    F_keep = 0
    big_count = 0
    for(i in 1:Tr) {
      gene_df$gene <- puck_d@counts[gene,cell_barc]/ puck_d@nUMI[cell_barc]
      if(i > 1)
        gene_df$gene <- sample(gene_df$gene)
      n = dim(gene_df)[1]
      fit <- loess(gene ~ x*y, span = 0.8, data = gene_df, degree = 1)
      #with(dat, plot(x, y, col = getcol(f), pch = 16, cex = 0.5, main = i))
      p <- fit$enp+1
      ss_total <-   sum((gene_df$gene- mean(gene_df$gene))^2)
      rss <- sum(fit$residuals^2)
      ss_reg <- ss_total - rss
      ms_reg <- ss_reg / (p-1)
      ms_res <- rss/(n-p)
      fstat <- ms_reg/ms_res
      if(i == 1) {
        gene_df$fitted <- fitted(fit)
        F_keep = fstat
        cvs[gene] <- sd(gene_df$fitted)/mean(gene_df$fitted)
        if(cvs[gene] < 0.5)
          break
      } else {
        if(fstat > F_keep)
          big_count = big_count + 1
      }
      if(big_count > 3)
        break
      #print(c(fstat, 1  - pf(fstat, p-1, n-p)))
      f_dist[i] <- fstat
    }
    pvals[gene] = (big_count + 1) / Tr
  }
  means_df <- data.frame(pvals,cvs)
  return(means_df)
}
means_df_CA3 <- get_int_genes(gene_df_CA3,"CA3", CA3_genes, puck_d, cell_barc_CA3)
int_genes_CA3 <- rownames(means_df_CA3[means_df_CA3$cvs >= 0.5 & means_df_CA3$pvals < 0.05,])
save(means_df_CA3,globallyvar, int_genes_CA3, puck_d, gene_df,file = "AnalysisPaper/Plotting/Results/2dloess.RData")
