gr <- readRDS("~/Dropbox/OvarianCellLines/structuralvar/data/segment/0cbs/CGOV44T.bam.rds")
gr <- keepSeqlevels(gr, c("chr5", "chr8"), pruning.mode="coarse")
saveRDS(gr, file="../inst/extdata/cgov44t_segments.rds")
