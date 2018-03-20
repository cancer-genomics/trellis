context("Merge adjacent amplicons")

test_that("joinNearGRanges", {
  library(GenomicRanges)
  gr <- GRanges(rep("chr1", 5),
                IRanges(c(1001, 3001, 5001, 7001, 9001), width=1000),
                seg.mean=c(0, 2.4, 2.5, 2.35, 0))
  names(gr) <-letters[1:5]
  params <- ampliconParams(MIN_SEGMEAN_DIFF=0.2)
  merged <- joinNearGRanges(gr, params)
  expected <- gr
  expect_identical(expected, merged) ## only merges if amplicon is at least 2kb
  params$MIN_WIDTH <- 500
  merged <- joinNearGRanges(gr, params)
  expected <- GRanges(rep("chr1", 3), IRanges(c(1001, 3001, 9001), width=c(1000, 5000, 1000)))
  expected$seg.mean <- c(0, 2.42, 0)
  names(expected) <- c("a", "", "e")
  expect_identical(merged, expected)
})
