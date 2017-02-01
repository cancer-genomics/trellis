## region segmentation results from ovarian cell lines
path <- file.path("~/Dropbox/OvarianCellLines",
                  "structuralvar/data/segment")
segs <- readRDS(file.path(path, "grl_hg19.rds"))
segs <- segs[["CGOV44T.bam"]]
segments <- segs
save(segments, file="../data/segments.rda")
