id <- "CGOV1T.bam"
id.rds <- paste0(id, ".rds")
library(GenomicRanges)
library(svovarian)
library(svclasses)
library(svrearrange)
library(svfusions)
build <- "hg19"
wd <- getwd()
odir <- "~/Dropbox/OvarianCellLines"
setwd(odir)
rfile <- file.path("rearrangements/3blat_unmapped", id.rds)
rlist <- readRDS(rfile)
r <- rlist[[2]]
saveRDS(r, file="../inst/extdata/cgov1t_complex_rearrangement.rds")
