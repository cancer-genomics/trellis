context("PreprocessViews:")
test_that("PreprocessViews", {
  library(Rsamtools)
  extdir <- system.file("extdata", package="Rsamtools", mustWork=TRUE)
  bamfile <- list.files(extdir, pattern="ex1.bam$", full.names=TRUE)
  bR <- GRanges(c("seq1", "seq2", "seq2"),
                IRanges(c(1000L, 100L, 1000L),
                        c(2000L, 1000L, 2000L)))
  bview <- BamViews(bamPaths=bamfile, bamRanges=bR)
  ##
  ## Use class that directly extends BamViews
  ##
  colnames(bview) <- "ex1"
  tree <- DataPaths(tempdir(), "Test", dryrun=FALSE)
  all(file.exists(tree))

  pv <- PreprocessViews2(bview)
})
