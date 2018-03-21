library(svovarian)
path <- projectOvarian()["rearrangements/3blat_unmapped"]
r <- readRDS(file.path(path, "CGOV7T.bam.rds"))
r <- r[chromosome(linkedBins(r)) == "chr11"]
rear_cgov7t <- r[chromosome(linkedBins(r)$linked.to) == "chr11"]
save(rear_cgov7t, file="~/Software/svpackages/svfusions/data/rear_cgov7t.rda")
