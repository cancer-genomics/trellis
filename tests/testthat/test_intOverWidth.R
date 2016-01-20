context("Intersection divided by width")
test_that("intOverWidth", {
  autosomes <- function() paste0("chr", 1:22)
  del <- GRanges(c("chr1", "chr2"), IRanges(c(1, 1), c(10, 10)))
  filters <- GRanges(c("chr1", "chr2"), IRanges(c(5, 8), c(20, 10000)))
  test <- intOverWidth(del, filters)
  expect_identical(test, c(0.6, 0.3))

  filters <- GRanges(c("chr1", "chrX"), IRanges(c(5, 8), c(20, 10000)))
  seqlevels(filters, force=TRUE) <- autosomes()
  seqlevels(del, force=TRUE) <- autosomes()
  test <- intOverWidth(del, filters)
  expect_identical(test, c(0.6, 0))

  del <- GRanges(c("chr1", "chr1"), IRanges(c(1, 101), c(10, 200)))
  filters <- GRanges("chr1", IRanges(1, 1000))
  test <- intOverWidth(del, filters)
  expect_identical(test, c(1, 1))

  cnv <- GRanges("chr1", IRanges(1, 5))
  filters <- GRanges("chr1", IRanges(c(1,5), c(3, 7)))
  fraction <- intOverWidth(cnv, filters)
  expect_equal(as.numeric(fraction), 0.8)

  filters <- GRanges("chr1", IRanges(1, 5))
  cnv <- GRanges("chr1", IRanges(c(1,5), c(3, 7)))
  fraction <- intOverWidth(cnv, filters)
  expect_equal(as.numeric(fraction), c(1, 1/3))

  ## cnv not spanned by any filter should have proportion 0
  filters <- GRanges("chr1", IRanges(1, 5))
  cnv <- GRanges(c("chr1", "chr1", "chr2"),
                 IRanges(c(1,5, 50), c(3, 7, 96)))
  fraction <- intOverWidth(cnv, filters)
  expect_equal(as.numeric(fraction), c(1, 1/3, 0))

  ## filters can be overlapping
  cnv <- GRanges(c("chr1", "chr1", "chr2"),
                 IRanges(c(1,5, 50), c(3, 7, 96)))
  filters <- GRanges("chr1", IRanges(c(1,3), c(5, 8)))
  expect_equal(as.numeric(intOverWidth(cnv, filters)), c(1, 1, 0))
})
