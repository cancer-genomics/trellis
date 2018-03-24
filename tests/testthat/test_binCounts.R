test_that("binnedCounts", {
  library(svbams)
  library(Rsamtools)
  library(svfilters.hg19)
  data(bins1kb)
  extdir <- system.file("extdata", package="svbams", mustWork=TRUE)
  bamfile <- file.path(extdir, "cgov10t.bam")
  ## Restrict to chr10 to make this faster
  bins <- keepSeqlevels(bins1kb, "chr3", pruning.mode="coarse")
  region <- GRanges("chr3", IRanges(59600000, 61000000), seqinfo=seqinfo(bins))
  bins <- subsetByOverlaps(bins, region)
  bviews <- BamViews(bamRanges=bins, bamPaths=bamfile)
  bins$cnt <- binnedCounts(bviews)
  ## Select a 'random' bin
  set.seed(123)
  select <- sample(seq_along(bins), 10)
  gr <- bins[select]

  flags <- scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=FALSE)
  scan_param <- ScanBamParam(flag=flags, which=gr)
  answer <- countBam(bamfile, param=scan_param)$records
  result <- bins$cnt[select]
  expect_identical(result, answer)
  ##
  ## binNormalize
  ##
  bins$std_cnt <- binNormalize(bins)
  bins$gc_adj <- binGCCorrect(bins)
})
