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
puck = readRDS(file.path(resultsdir,"puck_labeled.RDS"))
reference <- readRDS(paste(refdir,config$reffile,sep="/"))
cell_type_info <- get_cell_type_info(reference@assays$RNA@counts, reference@meta.data$liger_ident_coarse, reference@meta.data$nUMI)
cell_type_means = cell_type_info[[1]]; cell_type_names = cell_type_info[[2]]
n_cell_types = cell_type_info[[3]]

proportions <- read.csv(file.path(bulkdir,"weights.csv"))$Weight
names(proportions) = cell_type_names
print("callBeads: estimated bulk composition: ")
print(proportions)
#now we switch to a different gene list
gene_list = get_de_genes(cell_type_means, puck, fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg)
print(paste("callBeads: number of genes used for regression:", length(gene_list)))
cell_type_means_renorm <- get_norm_ref(puck, cell_type_means, gene_list, proportions)
marker_data = get_marker_data(cell_type_names, reference, cell_type_means_renorm, gene_list, marker_provided = TRUE)
cell_type = "Neurogenesis"
my_type = marker_data[marker_data$cell_type == cell_type,]
my_type$nPuck = rowSums(puck@counts[rownames(my_type),])
#plots a continuous value over the puck
plot_puck_continuous <- function(puck, barcodes, plot_val) {
  my_pal = pals::kovesi.rainbow(20)
  my_table = puck@coords[barcodes,]
  my_table$value = plot_val[barcodes]
  plot <- ggplot2::ggplot(my_table, ggplot2::aes(x=x, y=y)) + ggplot2::geom_point(ggplot2::aes(size = 1.5, shape=19,color=value)) +
    scale_colour_gradientn(colors = my_pal) + ggplot2::scale_shape_identity() + ggplot2::theme_bw() + ggplot2::scale_size_identity()
  plot
  #pdf(file.path(results_dir,"all_cell_types.pdf"))
  #invisible(print(plot))
  #dev.off()
}

plot_puck_wrapper <- function(puck, plot_val, cell_type = NULL, minUMI = 0) {
  my_cond = puck@nUMI > minUMI
  if(!is.null(cell_type))
    my_cond = my_cond & (puck@cell_labels == cell_type)
  plot_puck_continuous(puck, names(which(my_cond)), plot_val)
}

plot_puck_gene <- function(puck, gene, cell_type = NULL, minUMI = 0) {
  gene_vals = puck@counts[gene,]
  plot_puck_wrapper(puck, gene_vals, cell_type, minUMI)
}

plot_puck_wrapper(puck, log(puck@nUMI,2), cell_type, minUMI = 100)
plot_puck_gene(puck, "Marcksl1", cell_type = cell_type, minUMI =100)
proportions_sm = smooth_proportions(proportions, constrain = F, smoothing_par = 1e-6)
cell_type_means_renorm <- get_norm_ref(puck, cell_type_means, gene_list, proportions_sm)
cell_type_means_regood <- get_norm_ref(puck, cell_type_means, gene_list, proportions)
type_ind = which(puck@cell_labels == cell_type)
bulk_vec = rowSums(puck@counts[,type_ind])
total_UMI = sum(puck@nUMI[type_ind])
decompose(cell_type_means_soft, gene_list, total_UMI, bulk_vec, verbose=F)
decompose(cell_type_means_regood, gene_list, total_UMI, bulk_vec, verbose=F)
decompose(cell_type_means_renorm, gene_list, puck@nUMI[type_ind[4]], puck@counts[,type_ind[4]], verbose = F)
res = decompose(cell_type_means_renorm, gene_list, puck@nUMI[type_ind[ind]], puck@counts[,type_ind[ind]], verbose = F)
gene_means = rowMeans(cell_type_means[rownames(marker_data),])
print(get_marker_score_type(marker_data, bulk_vec,total_UMI,cell_type,gene_means))
print(get_marker_score_type(marker_data, puck@counts[,type_ind[2]],puck@nUMI[type_ind[2]],cell_type,gene_means))
for (cell_type_oth in cell_type_names)
  print(paste(cell_type_oth,(get_marker_score_type(marker_data, bulk_vec,total_UMI,cell_type_oth,gene_means))))

for (cell_type_oth in cell_type_names)
  print(paste(cell_type_oth,(get_marker_score_type(marker_data, puck@counts[,type_ind[ind]],puck@nUMI[type_ind[ind]],cell_type_oth,gene_means))))
ind = 1004
get_marker_score_type(marker_data, puck@counts[,type_ind[ind]],puck@nUMI[type_ind[ind]],cell_type,gene_means, score_mode = F)
curr_vec = puck@counts[,type_ind[ind]]
for(gene in rownames(my_type)) {
  if(curr_vec[gene] > 0) {
    val = curr_vec[gene]
    fake_vec = curr_vec
    fake_vec[gene] = 0
    res = decompose(cell_type_means_renorm, gene_list, puck@nUMI[type_ind[ind]], fake_vec, verbose = F)
    print(paste(gene, val, res["Neurogenesis"]))
  }
}
gene1 = "Tuba1a"
gene2 = "Rps14"
for (ind in 1000:1200) {
  curr_vec = puck@counts[,type_ind[ind]]
  if(curr_vec[gene1] > 0 || curr_vec[gene2] > 0 ) {
    reso = decompose(cell_type_means_renorm, gene_list, puck@nUMI[type_ind[ind]], curr_vec, verbose = F)
    val1 = curr_vec[gene1]; val2 = curr_vec[gene2]
    fake_vec = curr_vec
    fake_vec[gene1] = 0; fake_vec[gene2] = 0
    res = decompose(cell_type_means_renorm, gene_list, puck@nUMI[type_ind[ind]] - val1 - val2, fake_vec, verbose = F)
    print(paste(val1, val2, reso["Neurogenesis"], res["Neurogenesis"]))
  } else {
    print("None found")
  }
}

for (ind in 1000:1200) {
  curr_vec = puck@counts[,type_ind[ind]]
  res = decompose(cell_type_means_renorm, gene_list, puck@nUMI[type_ind[ind]], curr_vec, constrain=F)
  print(res)
}
for(gene in rownames(my_type)) {
  if(bulk_vec[gene] > 0) {
    val = bulk_vec[gene]
    fake_vec = bulk_vec
    fake_vec[gene] = 0
    res = decompose(cell_type_means_renorm, gene_list, total_UMI - val, fake_vec, verbose=F)
    print(paste(gene, val, res["Neurogenesis"]))
  }
}
complete_vec = rowSums(puck@counts)
gene_prop = bulk_vec / complete_vec
load_df = data.frame(gene_prop, bulk_vec)
sorted_df = load_df[order(gene_prop),]
tail(sorted_df[sorted_df$bulk_vec > 20, ])
tail(gene_prop[order(gene_prop)])
fake_vec["Rps14"] = 0
fake_vec["Tuba1a"] = 0
fake_vec["2410006H16Rik"] = 0

#soften
epsilon = 1e-9
alpha = 0.6
cell_type_means_soft = exp(log(cell_type_means_renorm + epsilon)*alpha)
ratio = (as.matrix(cell_type_means_soft) %*% proportions) / (as.matrix(cell_type_means_renorm) %*% proportions)
cell_type_means_soft = sweep(cell_type_means_soft, 1, ratio,'/')

curr_vec = puck@counts[,type_ind[ind]]
cUMI = puck@nUMI[type_ind[ind]]
pred = as.matrix(cell_type_means_renorm) %*% weights * cUMI
canames = get_marker_score_type(marker_data, puck@counts[,type_ind[ind]],puck@nUMI[type_ind[ind]],"CA3",gene_means, score_mode = F)
canames = tail(names(canames))
ngnames = get_marker_score_type(marker_data, puck@counts[,type_ind[ind]],puck@nUMI[type_ind[ind]],"Neurogenesis",gene_means, score_mode = F)
ngnames = tail(names(ngnames))
sml = sample(gene_list,1000)
decompose(cell_type_means_renorm, sml, puck@nUMI[type_ind[ind]], puck@counts[sml,type_ind[ind]], verbose = F)

pres_genes = names(which(curr_vec > 0))
zer_genes = names(which(curr_vec == 0))
for (cell_type_o in cell_type_names)
  print(paste(cell_type_o,sum(log(cell_type_means_renorm[pres_genes, cell_type_o]*cUMI))))

others = (1 - pred[zer_genes,]/cUMI)^cUMI
head(others[order(others)],20)
pred_comp = cell_type_means_renorm[gene_list,"Mural"]*cUMI/2 + cell_type_means_renorm[gene_list,"Endothelial_Stalk"]*cUMI /2

cell_type_o = cell_type

names(pred_comp) = gene_list
others = (1 - pred_comp[zer_genes]/cUMI)^cUMI
head(others[order(others)],20)
for (cell_type_o in cell_type_names) {
  print(cell_type_o)
  big_genes =  gene_list[which(cell_type_means_renorm[gene_list,cell_type_o] > 0.5/(cUMI))]
  #print(cell_type_means_renorm[big_genes,cell_type_o]*cUMI)
  #print(curr_vec[big_genes])
  print(sum(curr_vec[big_genes]) / sum(cell_type_means_renorm[big_genes,cell_type_o]*cUMI))
}

#test low sample regime
nUMI = 100
bead = bead_singlet(reference, gene_list, nUMI, cell_type)
decompose(cell_type_means, gene_list, nUMI, bead, constrain=F)
pred = cell_type_means_renorm[gene_list,cell_type]*nUMI
names(pred) = gene_list
pres_genes = gene_list[(which(bead > 0))]
pred[pres_genes]
bead[pres_genes,]
X = cell_type_means_renorm[gene_list,] * nUMI
b = bead[gene_list]
b[b > 1] = 1
write.csv(as.matrix(X),file.path(bulkdir,"X_bead.csv"))
write.csv(as.matrix(b),file.path(bulkdir,"b_bead.csv"))



decompose(cell_type_means, gene_list, nUMI, bead, constrain=F)

#generate the fake 1d data
p = .1
G = 1000
N = round(G*p)
bead = c(rep(1,N),rep(0,G-N))
dumb = c(rep(0,N),rep(1,G-N))/(G-N)*N
X = data.frame(rep(p,G),dumb)
decompose(X,1:G,N,bead,constrain=F)
decompose(X,1:G,N,bead,constrain=T)
pres_genes = which(bead > 0)
zer_genes = which(bead == 0)
q = .12
9*(-q) + 1*(log(1 - exp(-q)))
