context("AlignmentViews")
test_that("AlignmentViews", {
  library(Rsamtools)
  require(TestBams)
  extdata <- system.file("extdata", package="TestBams")
  bam.file <- list.files(extdata, pattern="\\.bam$", full.name=TRUE)
  bv <- BamViews(bam.file)
  dp <- DataPaths(tempdir())
  aviews <- AlignmentViews2(bv, dp)
  expect_identical(length(improperPaths(aviews)),
                   ncol(bv))
  expect_true(validObject(aviews[, 1]))
})
