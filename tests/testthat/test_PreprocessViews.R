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
  expect_true(file.exists(tree[["top"]]))
##
##  pview <- countExperiment2(bview, tree)
##  expect_identical(indexRanges(pview), 1:3)
##  pview2 <- pview[1, ]
##  expect_identical(nrow(pview2), 1L)
##  expect_identical(indexRanges(pview2), 1L)
##  pview2 <- pview[2:3, ]
##  expect_identical(nrow(pview2), 2L)
##  expect_identical(indexRanges(pview2), 2:3)
##  expect_identical(getScale(pview2), 1)
##  cnts <- assays(pview)
##
##  expect_is(seqinfo(pview), "Seqinfo")
##  expect_is(seqlevels(pview), "character")
##  test5 <- keepSeqlevels(pview, "seq1")
##  expect_identical(seqlevels(test5), "seq1")
##  expect_identical(indexRanges(test5), 1L)
})
