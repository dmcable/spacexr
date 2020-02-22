library(RCTD)
library(Matrix)
config <- config::get()
iv <- init_RCTD(load_info_renorm = T) #initial variables
cell_type_names = iv$cell_type_info[[2]]
results <- list()
for (fold_index in 1:config$n_puck_folds) {
  results <- append(results, readRDS(paste0(iv$slideseqdir,"/SplitPuckResults/results",fold_index,".RDS")))
}
metadir <- file.path(iv$slideseqdir,"MetaData")
meta_data <- readRDS(file.path(metadir,"meta_data.RDS"))
meta_df <- meta_data$meta_df
barcodes <- rownames(meta_df)
N <- meta_data$N_samples
weights = Matrix(0, nrow = N, ncol = iv$cell_type_info[[3]])
weights_doublet = Matrix(0, nrow = N, ncol = 2)
rownames(weights) = barcodes; rownames(weights_doublet) = barcodes
colnames(weights) = iv$cell_type_info[[2]]; colnames(weights_doublet) = c('first_type', 'second_type')
empty_cell_types = factor(character(N),levels = cell_type_names)
results_df <- data.frame(spot_class = factor(character(N),levels=c("reject", "singlet", "doublet_certain", "doublet_uncertain")),
                         first_type = empty_cell_types, second_type = empty_cell_types,
                         other_type = empty_cell_types, min_score = numeric(N), second_score = numeric(N),
                         conv_all = logical(N), conv_doublet = logical(N))
for(i in 1:N) {
  weights_doublet[i,] = results[[i]]$doublet_weights
  weights[i,] = results[[i]]$all_weights
  results_df[i, "spot_class"] = results[[i]]$spot_class
  results_df[i, "first_type"] = results[[i]]$first_type
  results_df[i, "second_type"] = results[[i]]$second_type
  results_df[i, "other_type"] = results[[i]]$other_type
  results_df[i, "min_score"] = results[[i]]$min_score
  results_df[i, "second_score"] = results[[i]]$second_score
  results_df[i, "conv_all"] = results[[i]]$conv_all
  results_df[i, "conv_doublet"] = results[[i]]$conv_doublet
}

class_df = data.frame(cell_type_names, row.names = cell_type_names)
colnames(class_df)[1] = "class"
class_df["Bergmann","class"] = "Astrocytes"
class_df["Fibroblast","class"] = "Endothelial"
class_df["MLI2","class"] = "MLI1"
class_df["Polydendrocytes","class"] = "Oligodendrocytes"
library(dplyr)
library(ggplot2)
common_cell_types = c("Astrocytes", "Bergmann", "Endothelial", "Fibroblast", "Golgi", "Granule", "MLI1", "MLI2", "Oligodendrocytes", "Polydendrocytes", "Purkinje", "UBCs")
resultsdir = file.path(iv$slideseqdir, "results")
n_levels = choose(meta_data$n_cell_types,2)
plots <- vector(mode = "list", length = n_levels)
index = 1
use_class = T
toss_reject = T
square_results = Matrix(0, nrow = meta_data$n_cell_types, ncol = meta_data$n_cell_types)
rownames(square_results) = common_cell_types[1:meta_data$n_cell_types]
colnames(square_results) = common_cell_types[1:meta_data$n_cell_types]
for (i in (1:(meta_data$n_cell_types-1))) {
  for (j in ((i + 1):meta_data$n_cell_types)) {
    type1 = common_cell_types[i]
    type2 = common_cell_types[j]
    curr_barcodes = meta_df$first_type == type1 & meta_df$second_type == type2
    if(toss_reject)
      curr_barcodes = curr_barcodes & (results_df$spot_class != "reject")
    if(use_class) {
      first_type_class = class_df[results_df[curr_barcodes,]$first_type,"class"]
      second_type_class = class_df[results_df[curr_barcodes,]$second_type,"class"]
      type1_class = class_df[type1,"class"]; type2_class = class_df[type2,"class"]
      first_type_pres <- first_type_class == type1_class | second_type_class == type1_class
      second_type_pres <- first_type_class == type2_class | second_type_class == type2_class
    } else {
      first_type_pres <- results_df[curr_barcodes,]$first_type == type1 | results_df[curr_barcodes,]$second_type == type1
      second_type_pres <- results_df[curr_barcodes,]$first_type == type2 | results_df[curr_barcodes,]$second_type == type2
    }
    plot_df <- aggregate(first_type_pres, list(meta_df[curr_barcodes, "first_UMI"]), mean)
    colnames(plot_df) = c("nUMI", "first_type")
    rownames(plot_df) = plot_df$nUMI
    plot_df$second_type = aggregate(second_type_pres, list(meta_df[curr_barcodes, "first_UMI"]), mean)$x
    plots[[index]] = ggplot2::ggplot() +
      ggplot2::geom_line(data = plot_df, aes(x = nUMI, y = first_type), color = "blue") +
      ggplot2::geom_line(data = plot_df, aes(x = nUMI, y = second_type), color = "red") +
      ggplot2::xlab('nUMI') +
      ggplot2::ylab('percent present') + labs(title = paste(type1,type2))
    square_results[type1,type2] = colMeans(plot_df[c("250","500","750"),])["first_type"]
    square_results[type2, type1] = colMeans(plot_df[c("250","500","750"),])["second_type"]
    index = index + 1
  }
}
if(use_class) {
  pdf(file.path(resultsdir,"doublet_classification_class.pdf"))
} else {
  pdf(file.path(resultsdir,"doublet_classification.pdf"))
}
invisible(lapply(plots, print))
dev.off()


type1 = "Astrocytes"
type2 = "Bergmann"


source('Plotting/figure_utils.R')
table(results_df[curr_barcodes,]$spot_class)
aggregate(meta_df$first_UMI, meta_df$second_type)

#next make the confusion matrix
true_types = unlist(list(meta_df[meta_df$first_UMI == 0,"second_type"], meta_df[meta_df$first_UMI == 1000,"first_type"]))
pred_types = unlist(list(results_df[meta_df$first_UMI == 0, "first_type"], results_df[meta_df$first_UMI == 1000, "first_type"]))
conf_mat <- caret::confusionMatrix(pred_types,factor(true_types,levels = iv$cell_type_info[[2]]))
plot_heat_map(conf_mat$table)

singlets = meta_df$first_UMI == 0 | meta_df$first_UMI == 1000
hist(weights_doublet[singlets,"first_type"],breaks = 30)
p1 <- hist(weights_doublet[singlets,"first_type"], breaks = 20)
p2 <- hist(weights_doublet[!singlets,"first_type"], breaks = 20)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T)
do.call(pmax, weights[singlets,])
apply(weights[singlets,],1,max)

p1 <- hist(apply(weights[singlets,],1,max), breaks = 20)
p2 <- hist(apply(weights[!singlets,],1,max), breaks = 20)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T)


