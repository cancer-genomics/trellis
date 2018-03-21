wd <- "/dcl01/scharpf1/data/rscharpf/dorothy/ovarian/rearrangements"
setwd(wd)
id <- "t_PGDX3310T_WGS_Ex.bam.rds.txt"
library(svclasses)
library(GenomicRanges)
library(svcnvs)
library(svfilters.hg18)
svpacks(); load_all("svalignments"); setwd(wd)
library(ggplot2)
suppressPackageStartupMessages(library(gridExtra))
library(scales)
##library(svrearrange)
svpacks(); load_all("svrearrange"); setwd(wd)
library(TxDb.Hsapiens.UCSC.hg18.refGene)
seqlengths(TxDb.Hsapiens.UCSC.hg18.refGene)
suppressPackageStartupMessages(library(svrearrange))
setwd("germline_filtered_split_read_supported_rearrangements/split_reads")
rlist <- readRDS(id)
rlist
if (length(rlist) == 0)  next()
##
## - temporary fix
##
for(j in seq_along(rlist)){
  r <- rlist[[j]]
  irp <- improper(r)
  si <- seqinfo(irp)
  genome(si) <- "hg18"
  seqlengths(si) <- seqlengths(TxDb.Hsapiens.UCSC.hg18.refGene)[seqlevels(si)]
  seqinfo(irp) <- si
  lb <- linkedBins(r)
  seqinfo(lb) <- si
  r@improper <- irp
  r@linkedBins <- lb
  rlist[[j]] <- r
}
##arrange each rearrangement in the rlist object by two orientations of possible rearranged pair
## no invalid split reads at this point
n.sr <- elementNROWS(splitReads(rlist))
if(all(n.sr < 1)) next()
rlist <- rlist[ n.sr > 0 ]
build <- "hg18"
near.coding <- seqJunctionNearTx(rlist, build)
rlist <- rlist[ near.coding ]
svpacks(); setwd("svrearrange/inst/extdata")
saveRDS(rlist, file="rlist_endometrioid_project.rds")
n.sr1 <- elementNROWS(splitReads(rlist))
if(length(rlist) < 1) next()
rlist2 <-  fiveTo3List(rlist, build)
## One of the rearrangements has both ++ and -- split reads
## - subsetting the split reads on the most common type results in a rearrangement object that is no longer supported by reads on both sides of the sequence junction
n.sr2 <- elementNROWS(splitReads(rlist2))

