context("Regression tests for candidateRearrangements")


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
  extdata <- system.file("extdata", package="svalignments")
  file.path(extdata, "improper_cgov44t.bam.rds")
}

cgov44t_preprocess<- function(){
  extdata <- system.file("extdata", package="svbams")
  id <- "cgov44t_revised.bam"
  bamfile <- file.path(extdata, id)

  cnvpath <- system.file("extdata", package="svcnvs")
  gr <- readRDS(file.path(cnvpath, "cgov44t_segments.rds"))
  segs <- keepSeqlevels(gr, "chr15", pruning.mode="coarse")

  irp.file <- file.path(extdata, "cgov44t_improper.rds")
  irp <- readRDS(irp.file)
  ddir <- system.file("extdata", package="svpreprocess",
                      mustWork=TRUE)
  lr <- readRDS(file.path(ddir, "preprocessed_coverage.rds"))/1000
  seqlevels(bins1kb, pruning.mode="coarse") <- paste0("chr", c(1:22, "X"))
  bins1kb$log_ratio <- lr
  del.gr <- segs[segs$seg.mean < hemizygousThr(DeletionParam())]
  proper.del <- properReadPairs(bamfile,
                                gr=reduce(del.gr,
                                          min.gapwidth=2000))
  rps <- list(improper=irp, proper_del=proper.del)
  pdat <- svcnvs::preprocessData(bam.file=bamfile,
                                 genome=genome(segs)[[1]],
                                 segments=segs,
                                 read_pairs=rps,
                                 bins=bins1kb)
}

test_that("findCandidates", {
  library(svfilters.hg19)
  library(svbams)
  library(svalignments)
  pdat <- cgov44t_preprocess()
  rparam <- RearrangementParams()
  rlist <- findCandidates2(pdat, rparam)
  if(FALSE){
    saveRDS(rlist, file="findCandidates.fe9b1f6.rds")
  }
  path <- system.file("extdata", package = "svrearrange")
  rlist.fe9b1f6 <- readRDS(file.path(path, "findCandidates.fe9b1f6.rds"))
  expect_equivalent(rlist.fe9b1f6, rlist)
})

test_that("seqJunctionsInferredByPairedTags", {
  library(svfilters.hg19)
  pdat <- cgov44t_preprocess()
  rparam <- RearrangementParams()
  candidates <- seqJunctionsInferredByPairedTags2(pdat, rparam)
  if(FALSE){
    saveRDS(candidates, file="seqJunctionsInferredByPairedTags.fe9b1f6.rds")
  }
  path <- system.file("extdata", package = "svrearrange")
  cand.fe9b1f6 <- readRDS(file.path(path, "seqJunctionsInferredByPairedTags.fe9b1f6.rds"))
  expect_equivalent(candidates, cand.fe9b1f6)
})

test_that("RearrangementList", {
  path <- system.file("extdata", package = "svrearrange")
  candidates <- readRDS(file.path(path, "seqJunctionsInferredByPairedTags.fe9b1f6.rds"))
  cand.list <- RearrangementList(candidates)
  if(FALSE){
    saveRDS(cand.list, file="RearrangementList.fe9b1f6.rds")
  }
  path <- system.file("extdata", package = "svrearrange")
  expected.fe9b1f6 <- readRDS(file.path(path, "RearrangementList.fe9b1f6.rds"))
  expect_identical(cand.list, expected.fe9b1f6)
})

test_that("type_each", {
  path <- system.file("extdata", package = "svrearrange")
  cand.list <- readRDS(file.path(path, "RearrangementList.fe9b1f6.rds"))
  result <- type_each(cand.list)
  if(FALSE){
    saveRDS(result, file="type_each.fe9b1f6.rds")
  }
  path <- system.file("extdata", package = "svrearrange")
  result.fe9b1f6 <- readRDS(file.path(path, "type_each.fe9b1f6.rds"))
  expect_identical(result, result.fe9b1f6)
})

test_that("modalRearrangement", {
  path <- system.file("extdata", package = "svrearrange")
  candidate.list <- readRDS(file.path(path, "type_each.fe9b1f6.rds"))
  colData(candidate.list)$modal_rearrangement <- modalRearrangement(candidate.list)
  colData(candidate.list)$percent_rearrangement <- percentRearrangement(candidate.list)
  path <- system.file("extdata", package = "svrearrange")
  rlist.fe9b1f6 <- readRDS(file.path(path, "findCandidates.fe9b1f6.rds"))
  expect_identical(candidate.list, rlist.fe9b1f6)
})
