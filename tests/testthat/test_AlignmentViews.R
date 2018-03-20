context("AlignmentViews")
test_that("AlignmentViews", {
  av <-  AlignmentViews2()
  expect_true(validObject(av))
  expect_identical(rdsId(av), character())

  library(Rsamtools)
  require(svbams)
  extdata <- system.file("extdata", package="svbams")
  bam.file <- list.files(extdata, pattern="cgov10t\\.bam$", full.name=TRUE)
  bv <- BamViews(bam.file)
  dp <- tempfile()
  aviews <- AlignmentViews2(bv, dp)
  ## allow non-existing files since this can create headaches
  ## expect_error(validObject(aviews))
  tmp.file <- tempfile()
  aviews <- AlignmentViews2(bv, tmp.file)
  expect_identical(length(improperPaths(aviews)),
                   ncol(bv))
  expect_true(validObject(aviews[, 1]))
})
