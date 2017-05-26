##extdata <- system.file("extdata", package="StructuralVariants")
updateAttributes <- function(x){
  cl <- attributes(x)$class
  attributes(cl)$package <- "svclasses"
  attributes(x)$class <- cl
  x
}
path <- "~/Software/StructuralVariants/inst/extdata"
fname <- file.path(path, "cgovt_filtered_rearrangements.rds")
rlist <- readRDS(fname)
rlist <- updateAttributes(rlist)
for(i in seq_along(rlist)){
  r <- rlist[[i]]
  r <- updateAttributes(r)
  r <- updateObject(r)
  rlist[[i]] <- r
}
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
