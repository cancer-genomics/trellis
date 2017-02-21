library(svfilters.hg19)
library(svbams)
library(Rsamtools)
library(svalignments)
data(germline_filters)
seq.info <- seqinfo(bins1kb)
extdata <- system.file("extdata", package="svbams")
bview <- BamViews(bamPaths=file.path(extdata, "cgov44t_revised.bam"))
seqinfo(bamRanges(bview)) <- seq.info

##
## Assumes improperly paired reads for the entire bam file have already been
## extracted and saved to disk
##
iparams <- improperAlignmentParams(mapqFilter=30)
irp <- getImproperAlignmentPairs(bview, iparams)
saveRDS(irp, file="../inst/extdata/improper_cgov44t.bam.rds")
##aview <- AlignmentViews2(bview, path="./")
##save(aview, file="../data/aview.rda")
