library(trellis)
data(rlist)
set.seed(123)
extdata <- system.file("extdata", package="svbams")
bamfile <- file.path(extdata, "cgov44t_revised.bam")
tags <- getSequenceOfReads(rlist, bamfile,
                           MAX=25L, build = "hg19")
save(tags, file="../data/tags.rda", compression_level=9)
