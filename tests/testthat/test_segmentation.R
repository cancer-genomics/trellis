context("Segmentation")

test_that("duplicate segments", {
  library(svcnvs)
  library(svbams)
  library(svfilters.hg19)
  library(Rsamtools)
  library(svpreprocess)
  path <- system.file("extdata", package="svbams", mustWork=TRUE)
  data(bins1kb)
  data(germline_filters, package="svfilters.hg19")
  bins <- keepSeqlevels(bins1kb, "chr3", pruning.mode="coarse")
  bins <- bins[seq(1, length(bins), by=5)]

  ##Not running this line causes the duplication behavior
  ##bins <- subsetByOverlaps(bins, GRanges("chr3", IRanges(59600000, 61000000)))
  bam.file <- file.path(path, "cgov10t.bam")
  bview <- BamViews(bamRanges=bins, bamPaths=bam.file)
  bins$cnt <- binCounts(bview)
  bins$std_cnt <- binNormalize(bins)
  set.seed(123)
  gc.adj <- binGCCorrect(bins)
  gc.adj <- gc.adj - 0.6 ## why?
  bins$log_ratio <- gc.adj
  seg.params <- SegmentParam()
  bins$adjusted <- bins$log_ratio
  ##trace(segmentBins, browser)
  g <- segmentBins(bins, seg.params)
  expect_true(!any(duplicated(g)))
})
