library(RCTD)
library(Matrix)
config <- config::get()
slideseqdir <- file.path("Data/Slideseq",config$slideseqfolder)
counts <- readr::read_csv(file = paste(slideseqdir,config$puckfile,sep="/"))
N_keep = 1000
restr_counts <- counts[,1:N_keep]
readr::write_csv(restr_counts, paste(slideseqdir,"smallDGE.csv",sep="/"))
