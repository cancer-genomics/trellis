context("Regression tests for filtering candidate rearrangements")

getBamView <- function(){
  seq.info <- seqinfo(bins1kb)
  extdata <- system.file("extdata", package="svbams")
  bview <- BamViews(bamPaths=file.path(extdata, "cgov44t_revised.bam"))
  seqinfo(bamRanges(bview)) <- seq.info
  bins <- keepSeqlevels(bins1kb, paste0("chr", c(1:22, "X")),
                        pruning.mode="coarse")
  bamRanges(bview) <- bins
  bview
}

improper.path <- function(){
  extdata <- system.file("extdata", package="svbams")
  file.path(extdata, "improper_cgov44t.bam.rds")
}

test_that("filterRearrangementList", {
  library(svfilters.hg19)
  library(svbams)
  if(FALSE){
    library(Rsamtools)
    data(germline_filters)
    bview <- getBamView()
    rlist <- readRDS("findCandidates.fe9b1f6.rds")
    rparams <- RearrangementParams()
    sl <- seqlevels(bamRanges(bview))
    gf <- reduceGenomeFilters(germline_filters, sl)
    rlist2 <- filterRearrangementList2(rlist,
                                       params=rparams,
                                       germline_filters=gf)
    if(FALSE){
      saveRDS(rlist2, file="filterRearrangementList.37be188.rds")
    }
  }
  path <- system.file("extdata", package="svbams")
  expected <- readRDS(file.path(path, "filterRearrangementList.37be188.rds"))
##  x <- expected
##  attributes(class(x))$package <- "trellis"
##  x1 <- x[[1]]
##  x2 <- x[[2]]
##  attributes(class(x1))$package <- "trellis"
##  attributes(class(x2))$package <- "trellis"
##  x[[1]] <- x1
##  x[[2]] <- x2
##  saveRDS(x, file=file.path(path, "filterRearrangementList.37be188.rds"))
  ##expect_identical(rlist.37be188, rlist2)
  rlist <- readRDS(file.path(path, "findCandidates.fe9b1f6.rds"))
  sl <- paste0("chr", c(1:22, "X"))
  gf <- reduceGenomeFilters(germline_filters, sl)
  rf <- rFilters(germline=gf)
  expect_identical(names(rf),
                   c("germline", "rear", "deletions", "amplicons"))
  expect_identical(length(rf$amplicons), 0L)
  rparams <- RearrangementParams()
  rdat <- rearrangementData(rlist=rlist,
                            read_pairs=list(),
                            filters=rf)
  expect_is(rdat, "list")
  rlist3 <- filterRear(rdat, params=rparams)
  expect_equivalent(rlist3, expected)
})

## this should be tested in svalignments
##test_that("getSequenceOfReads", {
##  path <- system.file("extdata", package="svbams")
##  rlist <- readRDS(file.path(path, "filterRearrangementList.37be188.rds"))
##  extdata <- system.file("extdata", package="svbams")
##  bam.file <- file.path(extdata, "cgov44t_revised.bam")
##  tags <- getSequenceOfReads(rlist, bam.file, MAX=25L)
##  if(FALSE){
##    saveRDS(tags, file="tags.34bd277.rds")
##  }
##  path <- system.file("extdata", package="svbams")
##  expected <- readRDS(file.path(path, "tags.34bd277.rds"))
##  tags <- tags[order(rownames(tags)), ]
##  expected <- expected[order(rownames(expected)), ]
##  expect_identical(tags, expected)
##})
