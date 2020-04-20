context("reviseEachJunction")

## this is not a regression
test_that("reviseEachJunction2", {
  path <- system.file("extdata", package = "svbams")
  sv <- readRDS(file.path(path, "reviseEachJunction.ca5194d.rds"))
  param <- DeletionParam()
  skip("Same GAlignments slots in class definition but not in object: elementType error")
  sv2 <- hemizygousBorders(sv, param)
  expect_identical(colnames(mcols(variant(sv2))),
                   colnames(mcols(variant(sv))))
})
