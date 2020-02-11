library(RCTD)
library(Matrix)
config <- config::get()
refdir <- file.path("Data/Reference",config$reffolder)
dir.create(file.path(refdir, "results"), showWarnings = FALSE) #folder to save results
slideseqdir <- file.path("Data/Slideseq",config$slideseqfolder)
resultsdir = file.path(slideseqdir,"results")
if(!dir.exists(resultsdir))
  dir.create(resultsdir)
bulkdir <- paste(slideseqdir,"results/Bulk",sep="/")
if(!dir.exists(bulkdir))
  dir.create(bulkdir)
puck = read.slideseq(slideseqdir, count_file = config$puckfile)
reference <- readRDS(paste(refdir,config$reffile,sep="/"))
cell_type_info <- get_cell_type_info(reference@assays$RNA@counts, reference@meta.data$liger_ident_coarse, reference@meta.data$nUMI)
cell_type_means = cell_type_info[[1]]; cell_type_names = cell_type_info[[2]]
n_cell_types = cell_type_info[[3]]

proportions <- read.csv(file.path(bulkdir,"weights.csv"))$Weight
names(proportions) = cell_type_names
print("callBeads: estimated bulk composition: ")
print(proportions)
#now we switch to a different gene list
gene_list = get_gene_list(cell_type_means, puck, cutoff_val = config$gene_cutoff_reg)
print(paste("callBeads: number of genes used for regression:", length(gene_list)))

X = read.csv(file.path(bulkdir,"X_bulk.csv"))
rownames(X) = X$X
X$X <- NULL
b = read.csv(file.path(bulkdir,"b_bulk.csv"))
rownames(b) = b$X
b$X <- NULL
diff = abs(log(b) - log(as.matrix(X) %*% as.matrix(proportions)))
rownames(diff)[diff > 3]
marker_data = get_marker_data(cell_type_names, reference, cell_type_means, gene_list, marker_provided = TRUE)
my_type = marker_data[marker_data$cell_type == "CA3",]
cell_type_means[rownames(my_type),"CA1"] /cell_type_means[rownames(my_type),"CA3"]
sum(apply(cell_type_means,1, max) > .0003)

#get the marker genes

total_gene_list = get_de_genes(cell_type_means, puck, fc_thresh = .5, expr_thresh = .0001)
cell_type_counts = table(reference@meta.data$liger_ident_coarse)

print(length(total_gene_list))
gene = "Fibcd1"
print(cell_type_means[gene,])
print(sum(total_gene_list == gene))

#debugging CA3
proportions = smooth_proportions(proportions, constrain = F)
cell_type_info_renorm = cell_type_info
cell_type_info_renorm[[1]] = get_norm_ref(puck, cell_type_info[[1]], gene_list, proportions)
cell_type_means_renorm = cell_type_info_renorm[[1]]
marker_data = get_marker_data(cell_type_names, reference, cell_type_means, gene_list, marker_provided = TRUE)
gene_means = rowMeans(cell_type_means[rownames(marker_data),])
search_x = 3500
search_y = 4000
search_dist = 100
search_res = puck@coords[abs(puck@coords$x - search_x) < search_dist & abs(puck@coords$y - search_y) < search_dist,]
search_res = rownames(search_res)
box_vec = rowSums(puck@counts[,search_res])
nUMI = sum(puck@nUMI[search_res])
print(get_marker_score_type(marker_data, box_vec,nUMI,"Neurogenesis",gene_means))
print(get_marker_score_type(marker_data, box_vec,nUMI,"CA3",gene_means))
decompose(cell_type_means_renorm, gene_list_good, nUMI, box_vec[gene_list_good], verbose=T)
ind = 6
decompose(cell_type_means_renorm, gene_list_good, puck@nUMI[search_res[ind]], puck@counts[gene_list_good,search_res[ind]], verbose=T)
print(get_marker_score_type(marker_data, puck@counts[,ind],puck@nUMI[ind],"Neurogenesis",gene_means))
print(get_marker_score_type(marker_data, puck@counts[,ind],puck@nUMI[ind],"CA3",gene_means))

wrong_gene = "Gm2000"
#debug getNormRef
bulk_vec = rowSums(puck@counts)
weight_avg = rowSums(sweep(cell_type_means[gene_list,],2,proportions,'*'))
target_means = bulk_vec[gene_list]/sum(bulk_vec)
cell_type_means_renorm = sweep(cell_type_means[gene_list,],1,weight_avg / target_means,'/')

gene_list_good = gene_list[bulk_vec[gene_list] > 10]

marker_data_rn = get_marker_data(cell_type_names, reference, cell_type_means_renorm, gene_list_good, marker_provided = TRUE)
gene_means_rn = rowMeans(cell_type_means_renorm[rownames(marker_data_rn),])
print(get_marker_score_type(marker_data_rn, box_vec,nUMI,"Neurogenesis",gene_means_rn))
print(get_marker_score_type(marker_data_rn, box_vec,nUMI,"CA3",gene_means_rn))
dev_df = plot_deviance_hist(box_vec,"CA3", cell_type_means_renorm, gene_list, nUMI, to_disp = 10, verbose = F)
marker_data[rownames(dev_df),]


#bayesian decomposiiotn
bead1 <- bead_mix(reference, gene_list, 1000, 1000, "CA3", "Denate")
test_ref = reference
UMI1 = 2000; UMI2 = 2000
type1 = "CA3"; type2 = "Denate"
nUMI = test_ref@meta.data$nUMI
firstInd = sample(intersect(which(test_ref@meta.data$liger_ident_coarse == type1) , which(nUMI > 1000)),1)
secondInd = sample(intersect(which(test_ref@meta.data$liger_ident_coarse == type2) , which(nUMI > 1000)),1)
cell1 = sub_sample_cell(gene_list, test_ref@assays$RNA@counts, firstInd, UMI1)
cell2 = sub_sample_cell(gene_list, test_ref@assays$RNA@counts, secondInd, UMI2)
bead = data.matrix((as.matrix(cell1 + cell2)))
res <- decompose_sparse(cell_type_means, gene_list, 2000, bead, "CA3", "Denate")
p_1 = res[1]
p_2 = res[2]
int_genes = tail(rownames(bead)[order(bead)],40)
N_genes = length(int_genes)
expect_1 = vector(mode="numeric",length = N_genes)
expect_2 = vector(mode="numeric",length = N_genes)
variance = vector(mode="numeric",length = N_genes)
sq_err_model = 0
sq_err_null = 0
sq_err_complete = 0
for(ind in 1:N_genes) {
  gene = int_genes[ind]
  denom = p_1 * cell_type_means_renorm[gene,type1] + p_2 * cell_type_means_renorm[gene,type2]
  posterior_1 = p_1 * cell_type_means_renorm[gene,type1] / denom
  expect_1[[ind]] = posterior_1 * bead[gene,]
  expect_2[[ind]] = bead[gene,] - posterior_1 * bead[gene,]
  variance[[ind]] = posterior_1 * bead[gene,] * (1 - posterior_1)
  k = cell1[gene]; n = bead[gene,]
  #log_l_model = log_l_model + log(choose(n,k)) + k*log(posterior_1) + (n-k)*log(1-posterior_1)
  #log_l_null = log_l_null + log(choose(n,k)) + k*log(0.5) + (n-k)*log(0.5)
  #log_l_complete = log_l_complete + log(choose(n,k)) + k*log(k/n + epsilon) + (n-k)*log(1-k/n + epsilon)
  sq_err_model = sq_err_model + (posterior_1 - k/n)^2
  sq_err_null = sq_err_null + (0.5 - k/n)^2
  sq_err_complete = sq_err_complete + (k/n)*(1 - k/n)/n
}
pseudo_R2 = 1 - (sq_err_model - sq_err_complete) / (sq_err_null - sq_err_complete)
#pseudo_R2 = 1 - (log_l_model - log_l_complete) / (log_l_null - log_l_complete)
st_dev = sqrt(variance) * const_stderr
const_stderr = 1.96
cell_1 = as.vector(cell1[int_genes])
Total_UMI = bead[int_genes,]
plot_df <- data.frame(int_genes, expect_1,cell_1,st_dev, Total_UMI)
plot_df <- plot_df[order(expect_1),]
plot_df$int_genes <- factor(plot_df$int_genes, levels= plot_df$int_genes)
p<- ggplot2::ggplot(plot_df, ggplot2::aes(x=int_genes, y=cell_1, colour = "cell_1")) +
  ggplot2::geom_point(size=3)+
  ggplot2::geom_point(size=3,ggplot2::aes(y=Total_UMI,colour = "all_UMI")) +
  ggplot2::geom_point(size=3, ggplot2::aes(y=expect_1,colour = "expect_1")) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin=expect_1-st_dev, ymax=expect_1+st_dev,colour = "expect_1"), width=.05,
                         position=ggplot2::position_dodge(0.05))
#pseudo R2
log_l_model = Total_UMI*log(posterior_1_vec) -
print(p)
if(save.file) {
  out_dir = paste(mydir,"Output/CompositionRecovery",sep="/")
  saveRDS(plot_df,file = paste(out_dir,"weight_recovery_df.RDS",sep="/"))
}
