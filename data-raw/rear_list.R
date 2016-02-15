extdata <- system.file("extdata", package="StructuralVariants")
rlist <- readRDS(file.path(extdata, "cgovt_filtered_rearrangements.rds"))
cl <- attributes(rlist)$class
attributes(cl)$package <- "svclasses"
attributes(rlist)$class <- cl
rear_list <- rlist
save(rear_list, file="~/Software/svpackages/svclasses/data/rear_list.rda")

library(TxDb.Hsapiens.UCSC.hg19.refGene)
library(BSgenome.Hsapiens.UCSC.hg19)
txdb <- TxDb.Hsapiens.UCSC.hg19.refGene
genome <- BSgenome.Hsapiens.UCSC.hg19
tx <- transcripts(txdb)
cds.all <- cdsBy(txdb, "tx", use.names=TRUE)
r <- rear_list[["18557-18736"]]
rear_cds <- svfusions:::getCDS(r, tx, cds.all)
save(rear_cds, file="~/Software/svpackages/svclasses/data/rear_cds.rda")

clipped <- svfusions::clip(transcripts)
fuse(clipped)
