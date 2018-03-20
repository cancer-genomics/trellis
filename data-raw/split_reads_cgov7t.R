library(svbams)
library(svfilters.hg19)
library(svcnvs)
library(svalignments)
library(svrearrange)
library(svfusions)
data(rear_cgov7t)
## no bam file available for cgov7t
x <- readRDS("fusionList.641eb81.rds")
sr.file <- file.path("~/Dropbox/OvarianCellLines/structuralvar",
                     "data/alignments/5parsed_unmapped",
                     "CGOV7T.bam.rds")
rid <- unique(x$rearrangement.id)
split_reads <- readRDS(sr.file)
split_reads <- split_reads[rid]
saveRDS(split_reads, file="../inst/extdata/cgov7t_split_reads.rds")
