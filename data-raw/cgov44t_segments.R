gr <- readRDS("~/Dropbox/OvarianCellLines/structuralvar/data/segment/0cbs/CGOV44T.bam.rds")
segs <- keepSeqlevels(gr, c("chr5", "chr8", "chr15"), pruning.mode="coarse")
##gr <- readRDS("~/Dropbox/OvarianCellLines/structuralvar/data/segment/0cbs/CGOV44T.bam.rds")
##gr <- keepSeqlevels(gr, c("chr5", "chr8"), pruning.mode="coarse")
saveRDS(segs, file="../inst/extdata/cgov44t_segments.rds")

## for unit test
library(Rsamtools)
library(svfilters.hg19)
library(svpreprocess)
data(germline_filters, package="svfilters.hg19")
extdata <- system.file("extdata", package="svbams")
##id <- "CGOV44T.bam"
id <- "cgov44t_revised.bam"
id.rds <- paste0(id, ".rds")
bamfile <- file.path(extdata, id)
bview <- Rsamtools::BamViews(bamPaths=bamfile)
data(segments, package="svcnvs")
data(deletion, package="svcnvs")
gr <- variant(deletion)
gr <- expandGRanges(gr, 10000)
segs <- keepSeqlevels(segments, "chr15", pruning.mode="coarse")
segs <- segs[overlapsAny(segs, gr)]
seqlevels(segs) <- seqlevels(gr)
seqinfo(segs) <- seqinfo(gr)
saveRDS(segs, file="../tests/testthat/segs.4adcc78.rds")
