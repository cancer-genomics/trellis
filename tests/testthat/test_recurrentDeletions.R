context("recurrent deletions")
testthat("recurrentDeletions", {
  genes <- GRanges("chr1", IRanges(5, 10))
  genes$gene_name <- "a"
  gr1 <- GRanges(rep("chr1", 2), IRanges(c(4, 8), c(6, 10)), id=rep("id1", 2))
  gr2 <- GRanges("chr1", IRanges(3, 9), id="id2")
  grl <- GRangesList(id1=gr1, id2=gr2)
  ans <- recurrentDeletions(genes, grl, maxgap=5)
  ## should only count first sample once
  expect_identical(ans$freq, 2L)

  ## gaps allowed
  genes <- GRanges("chr1", IRanges(50e3, 55e3))
  genes$gene_name <- "a"
  gr1 <- GRanges(rep("chr1", 2), IRanges(c(48e3, 60.1e3), c(49e3, 65e3)), id=rep("id1", 2))
  gr2 <- GRanges("chr1", IRanges(50e3, 50e3), id="id2")
  grl <- GRangesList(id1=gr1, id2=gr2)
  ans <- recurrentDeletions(genes, grl, maxgap=5000)
  expect_identical(ans$freq, 2L)

  ## genes that are not recurrent (count of 1) are not returned
  ans <- recurrentDeletions(genes, grl, maxgap=0)
  expect_identical(ans$freq, integer())
})
