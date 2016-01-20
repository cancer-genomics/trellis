test_that("expandGRanges", {
  library(GenomeInfoDb)
  g_tmp <- GRanges("chr1", IRanges(50, 75))
  seqlengths(g_tmp) <- 5e6
  ##trace(expandGRanges, browser)
  test <- expandGRanges(g_tmp, 2*width(g_tmp))
  ans <- GRanges("chr1", IRanges(1, 124))
  seqlengths(ans) <- seqlengths(test)
  expect_identical(ans, test)
})
