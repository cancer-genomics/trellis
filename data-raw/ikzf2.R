library(svfilters.hg19)
library(svrearrange)
fusion.tx <- c(transcripts[transcripts$gene_name=="ERBB4"],
               transcripts[transcripts$gene_name=="IKZF2"])
fusion.tx <- reduce(fusion.tx)
fusion.tx$gene <- c("ERBB4", "IKZF2")



rlist <- readRDS("~/Dropbox/OvarianCellLines/rearrangements/3blat_unmapped/CGOV11T_1.bam.rds")
lb <- linkedBins(rlist)
index1 <- findOverlaps(fusion.tx, lb)
index2 <- findOverlaps(fusion.tx, lb$linked.to)
index <- unique(c(subjectHits(index1), subjectHits(index2)))
rears <- rlist[index]
r <- rears[[2]]
saveRDS(r, file="../inst/extdata/ikzf2_rear.rds")
