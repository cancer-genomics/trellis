test_that("binGCCorrect", {
  library(Rsamtools)
  library(svbams)
  library(svfilters.hg19)
  data(bins1kb)
  extdir <- system.file("extdata", package="svbams", mustWork=TRUE)
  bamfile <- file.path(extdir, "cgov10t.bam")
  bins <- keepSeqlevels(bins1kb, "chr3", pruning.mode="coarse")
  bins <- subsetByOverlaps(bins, GRanges("chr3", IRanges(59600000, 61000000)))
  bviews <- BamViews(bamRanges=bins, bamPaths=bamfile)
  bins$cnt <- binCounts(bviews)
  std_cnt <- binNormalize(bins)
  bins$std_cnt <- std_cnt
  correct <- binGCCorrect(bins)
  expect_that(length(correct), equals(1378))
  width <- width(bins)
  expect_that(width[1], equals(1000))
  expect_that(correct, is_a("numeric"))
}
)


