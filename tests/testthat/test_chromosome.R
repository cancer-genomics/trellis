context("chromosome")
test_that("chromosome", {
  gr <- GRanges(c("chr1", "chrX"), IRanges(1:2, 3:4))
  chr <- chromosome(gr)
  expect_identical(chr, c("chr1", "chrX"))

  library(SummarizedExperiment)
  se <- SummarizedExperiment()
  rowRanges(se) <- gr
  expect_identical(chromosome(se), c("chr1", "chrX"))
})
