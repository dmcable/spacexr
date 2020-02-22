#this script converts a DGE and coords matrix to a puck object saved as RDS
library(RCTD)
library(Matrix)
print("compressPuck: begin")
config <- config::get()
slideseqdir <- file.path("Data/Slideseq",config$slideseqfolder)
puck = read.slideseq(slideseqdir, count_file = config$puckfile)
split_puck(puck, slideseqdir, config$n_puck_folds)
saveRDS(puck, file.path(slideseqdir, config$puckrds))
split_puck(puck, slideseqdir, config$n_puck_folds)
# to read: puck = readRDS(file.path(slideseqdir, config$puckrds))


