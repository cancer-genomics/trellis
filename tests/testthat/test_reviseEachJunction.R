context("reviseEachJunction")

## this is not a regression
test_that("reviseEachJunction2", {
  sv <- readRDS("reviseEachJunction.ca5194d.rds")
  sv2 <- hemizygousBorders(sv, param)
  expect_identical(colnames(mcols(variant(sv2))),
                   colnames(mcols(variant(sv))))
})
