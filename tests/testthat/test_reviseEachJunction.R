context("reviseEachJunction")

## this is not a regression
test_that("reviseEachJunction2", {
  sv.bak <- sv <- readRDS("reviseEachJunction.ca5194d.rds")
  sv2 <- hemizygousBorders(sv, param)

  sv3 <- hemizygousBorders2(sv, param)
  expect_identical(colnames(mcols(variant(sv2))),
                   colnames(mcols(variant(sv.bak))))

  expect_identical(granges(variant(sv2)),
                   granges(variant(sv3)))
  expect_identical(colnames(mcols(variant(sv3))),
                   colnames(mcols(variant(sv.bak))))
})
