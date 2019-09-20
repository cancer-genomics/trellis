library(tidyverse)
library(ggplot2)
library(grid)
library(trellis)
library(eschelman.data)
library(GenomicRanges)
library(scales)
wd <- getwd()
setwd("~/Dropbox/Labs/Eschelman/eschelman.data/data-derived")
data(bamviews)
ids <- paste0(colnames(bamviews), ".rds")
i <- 1
extdir <- "../inst/extdata"
segdir <- file.path(extdir, "segments")
bindir <- file.path(extdir, "bins")
deldir <- file.path(extdir, "deletions")
rdir   <- file.path(extdir, "rearrangements_blat")
patid  <- substr(ids, 1, 5)
select <- dplyr::select
fnames <- tibble(patient=patid,

                 labid  =ids) %>%
  group_by(patient) %>%
  summarize(samples=paste(labid, collapse=","),
            cancer=paste0(unique(patient), "C.rds"),
            matched=paste0(unique(patient), "N.rds"),
            n=n()) %>%
  filter(n == 2) %>%
  select(patient, cancer, matched)
fnames <- fnames[i, ]
pid <- paste0(fnames$patient[[1]], "C")
figdir <- file.path("../../figures/rearrangements",
                    pid)
if(!dir.exists(figdir)) dir.create(figdir)
rfileN <- file.path(rdir, fnames$matched)
rfileC <- file.path(rdir, fnames$cancer)

rlistC <- readRDS(rfileC)
rlistN <- readRDS(rfileN)

## Any sequence junctions in cancer are also in the normal sample
is_overlap <- overlapsAny(linkedBins(rlistC), linkedBins(rlistN),
                           maxgap=500) &
  overlapsAny(linkedTo(rlistC), linkedTo(rlistN))
if(any(is_overlap)){
  ## some of these would be thrown out anyway if no split reads were found
  tmpC <- rlistC[ is_overlap ]
  hits <- findOverlaps(linkedBins(tmpC), linkedBins(rlistN),
                       maxgap=500)
  tmpN <- rlistN[ subjectHits(hits) ]
}

rlistC <- rlistC[ !is_overlap ]
##
## this requirement eliminates many candidates
##
at_least_one <- lengths(splitReads(rlistC)) >= 1
rlistC <- rlistC[ at_least_one ]
rlist2 <- fiveTo3List(rlistC, build="hg19")
rlist3 <- rlist2[ is_valid_splits(rlist2, maxgap=50) ]
fig.list <- vector("list", length(rlist3))
figName <- function(r){
  type <- modalRearrangement(r)
  nm <- trellis:::conciseGRangeSummary(linkedBins(r), type)
  nm
}
fignames <- file.path(figdir, paste0(sapply(rlist3, figName), ".pdf"))
j <- 7

rear <- rlist3[[j]]
chr6.12_rear <- rear
setwd(wd) ## back to trellis directory
saveRDS(chr6.12_rear, "../inst/extdata/chr6-12_rear.rds")
